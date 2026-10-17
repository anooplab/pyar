#!/usr/bin/env python3
"""Run the frozen MS0 matrix in a caller-selected scratch directory."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import time


ROOT = Path(__file__).resolve().parent
FIXTURES = ROOT / "fixtures"
SOLUTES = {
    "methane": {"atoms": 5, "charge": 0, "multiplicity": 1},
    "ammonia": {"atoms": 4, "charge": 0, "multiplicity": 1},
    "formamide": {"atoms": 6, "charge": 0, "multiplicity": 1},
    "acetate": {"atoms": 7, "charge": -1, "multiplicity": 1},
}


def run_case(workdir: Path, solute: str, count: int, methods: tuple[str, ...]) -> None:
    case = workdir / f"{solute}-water-n{count}"
    if (case / "run.json").exists():
        raise FileExistsError(f"Run manifest already exists; choose a fresh workdir: {case / 'run.json'}")
    case.mkdir(parents=True, exist_ok=True)
    inputs = case / "inputs"
    inputs.mkdir(exist_ok=True)
    for name in (solute, "water"):
        source = FIXTURES / f"{name}.xyz"
        target = inputs / f"{name}.xyz"
        if target.exists() and target.read_bytes() != source.read_bytes():
            raise FileExistsError(f"Refusing to overwrite different benchmark input: {target}")
        if not target.exists():
            shutil.copyfile(source, target)
    solute_state = SOLUTES[solute]
    record = {
        "solute": solute, "solvent": "water", "count": count,
        "expected_atoms": solute_state["atoms"] + 3 * count,
        "charge": solute_state["charge"], "multiplicity": solute_state["multiplicity"], "methods": {},
    }
    base = ["--count", str(count), "--orientations", "8",
            "--maximum-number-of-seeds", "4", "--surface-points-per-atom", "32"]
    for method in methods:
        out = case / method
        if out.exists() and any(out.iterdir()):
            record["methods"][method] = {"status": "preserved_existing_output"}
            continue
        out.mkdir(exist_ok=True)
        log = out / "run.log"
        if method in {"pyar-geometry", "pyar-xtb"}:
            cmd = ["pyar", "microsolvate", str(inputs / f"{solute}.xyz"),
                   str(inputs / "water.xyz"), *base, "--charge", str(solute_state["charge"]), "0",
                   "--multiplicity", str(solute_state["multiplicity"]), "1", "--output", str(out / "run")]
            if method == "pyar-xtb":
                cmd += ["--backend", "xtb", "--nprocs", "1"]
        elif method == "crest-qcg-grow":
            cmd = ["crest", str(inputs / f"{solute}.xyz"), "-qcg", str(inputs / "water.xyz"),
                   "-nsolv", str(count), "-grow", "-nofix", "-T", "1", "-gfn2", "-samerand"]
            if solute_state["charge"]:
                cmd += ["-chrg", str(solute_state["charge"])]
            if solute_state["multiplicity"] != 1:
                cmd += ["-uhf", str(solute_state["multiplicity"] - 1)]
        elif method == "crest-qcg-ensemble":
            cmd = ["crest", str(inputs / f"{solute}.xyz"), "-qcg", str(inputs / "water.xyz"),
                   "-nsolv", str(count), "-ensemble", "-md", "-mdlen", "0.2", "-mddump", "100",
                   "-nofix", "-T", "1", "-gfn2", "-samerand", "-keepdir"]
            if solute_state["charge"]:
                cmd += ["-chrg", str(solute_state["charge"])]
            if solute_state["multiplicity"] != 1:
                cmd += ["-uhf", str(solute_state["multiplicity"] - 1)]
        elif method == "crest-conformers":
            cluster = case / "crest-qcg-grow" / "grow" / "cluster.xyz"
            if not cluster.is_file():
                record["methods"][method] = {
                    "status": "not_run", "reason": "requires completed crest-qcg-grow cluster.xyz",
                }
                continue
            shutil.copyfile(cluster, inputs / "crest-qcg-cluster.xyz")
            cmd = ["crest", str(inputs / "crest-qcg-cluster.xyz"), "-gfn2", "-T", "1"]
        elif method == "orca-solvator":
            xyz = (inputs / f"{solute}.xyz").read_text().splitlines()
            coords = "\n".join(xyz[2:])
            inp = out / "solvator.inp"
            inp.write_text(
                "! GFN2-XTB ALPB(WATER)\n%pal nprocs 1 end\n"
                f"%SOLVATOR\n  NSOLV {count}\n  FIXSOLUTE FALSE\n  VACUUMSEARCH TRUE\nEND\n"
                f"* XYZ {solute_state['charge']} {solute_state['multiplicity']}\n{coords}\n*\n"
            )
            orca = os.environ.get("ORCA", "orca")
            cmd = [orca, inp.name]
        else:
            raise ValueError(f"Unknown method: {method}")
        executable = shutil.which(cmd[0])
        if not executable:
            record["methods"][method] = {"status": "not_run", "reason": f"{cmd[0]} not found"}
            continue
        env = os.environ.copy()
        env["OMP_NUM_THREADS"] = "1"
        if method == "orca-solvator":
            env.setdefault("XTBEXE", shutil.which("xtb") or "xtb")
        started = time.monotonic()
        try:
            with log.open("w") as handle:
                result = subprocess.run(cmd, cwd=out, env=env, text=True,
                                        stdout=handle, stderr=subprocess.STDOUT,
                                        timeout=3600, check=False)
            record["methods"][method] = {
                "status": "completed" if result.returncode == 0 else "failed",
                "exit_code": result.returncode,
                "duration_seconds": round(time.monotonic() - started, 3),
                "command": cmd,
                "log": str(log),
            }
        except subprocess.TimeoutExpired:
            record["methods"][method] = {
                "status": "failed", "reason": "3600 s timeout",
                "duration_seconds": 3600, "log": str(log),
            }
    (case / "run.json").write_text(json.dumps(record, indent=2) + "\n")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workdir", required=True, type=Path, help="empty or previously populated scratch root")
    parser.add_argument("--counts", nargs="+", type=int, default=[4])
    parser.add_argument("--solutes", nargs="+", choices=tuple(SOLUTES), default=list(SOLUTES))
    parser.add_argument("--methods", nargs="+", choices=(
        "pyar-geometry", "pyar-xtb", "crest-qcg-grow", "crest-qcg-ensemble", "crest-conformers", "orca-solvator",
    ),
                        default=["pyar-geometry", "pyar-xtb", "crest-qcg-grow", "crest-conformers",
                                 "crest-qcg-ensemble", "orca-solvator"])
    args = parser.parse_args(argv)
    if any(count < 1 for count in args.counts):
        parser.error("water counts must be positive")
    args.workdir = args.workdir.expanduser().resolve()
    args.workdir.mkdir(parents=True, exist_ok=True)
    for solute in args.solutes:
        for count in args.counts:
            run_case(args.workdir, solute, count, tuple(args.methods))
    failures = []
    for manifest in args.workdir.glob("*/run.json"):
        record = json.loads(manifest.read_text())
        for method, result in record["methods"].items():
            if result.get("status") == "failed":
                failures.append(f"{manifest.parent.name}/{method}")
    print(f"failed_methods={len(failures)}" + (f": {', '.join(failures)}" if failures else ""))
    print(args.workdir)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
