#!/usr/bin/env python3
"""Common-vacuum GFN2-xTB single-point rescoring for benchmark XYZ pools."""

from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path
import re
import shutil
import subprocess


ENERGY = re.compile(r"TOTAL ENERGY\s+(-?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)\s+Eh")
ROOT = Path(__file__).resolve().parent
SOLUTE_STATES = {"methane": (0, 1), "ammonia": (0, 1), "formamide": (0, 1), "acetate": (-1, 1)}


def read_xyz(path: Path) -> list[str]:
    lines = path.read_text().splitlines()
    try:
        count = int(lines[0].strip())
    except (IndexError, ValueError) as exc:
        raise ValueError(f"Invalid XYZ atom count: {path}") from exc
    if count < 1 or len(lines) < count + 2:
        raise ValueError(f"Truncated XYZ: {path}")
    symbols = [line.split()[0] for line in lines[2:count + 2] if len(line.split()) >= 4]
    if len(symbols) != count:
        raise ValueError(f"Expected {count} atom records: {path}")
    return symbols


def read_xyz_frames(path: Path) -> list[tuple[list[str], str]]:
    """Read every frame in an XYZ trajectory/ensemble without losing its order."""
    lines = path.read_text().splitlines()
    frames = []
    index = 0
    while index < len(lines):
        if not lines[index].strip():
            index += 1
            continue
        try:
            count = int(lines[index].strip())
        except ValueError as exc:
            raise ValueError(f"Invalid XYZ frame atom count at line {index + 1}: {path}") from exc
        end = index + count + 2
        if count < 1 or end > len(lines):
            raise ValueError(f"Truncated XYZ frame at line {index + 1}: {path}")
        frame = lines[index:end]
        symbols = [line.split()[0] for line in frame[2:] if len(line.split()) >= 4]
        if len(symbols) != count:
            raise ValueError(f"Expected {count} atom records in frame at line {index + 1}: {path}")
        frames.append((symbols, "\n".join(frame) + "\n"))
        index = end
    if not frames:
        raise ValueError(f"XYZ contains no coordinate frames: {path}")
    return frames


def count_water_fragments(symbols: list[str], solute_atom_count: int) -> int:
    """Validate the documented solute-prefix and O,H,H solvent fragment map."""
    solvent = symbols[solute_atom_count:]
    if not solvent:
        return 0
    if len(solvent) % 3:
        raise ValueError("Candidate atom count does not contain whole water fragments")
    fragments = [solvent[index:index + 3] for index in range(0, len(solvent), 3)]
    if any(fragment != ["O", "H", "H"] for fragment in fragments):
        raise ValueError("Candidate solvent atoms do not preserve O,H,H fragment order")
    return len(fragments)


def match_solute_prefix(symbols: list[str]) -> str | None:
    """Identify a benchmark solute by its fixed atom-order prefix."""
    return next((name for name in SOLUTE_STATES
                 if symbols[:len(read_xyz(ROOT / "fixtures" / f"{name}.xyz"))]
                 == read_xyz(ROOT / "fixtures" / f"{name}.xyz")), None)


def energy_for(xyz: Path, executable: str, scratch: Path, *, charge: int = 0,
               multiplicity: int = 1, geometry_text: str | None = None) -> float:
    scratch.mkdir(parents=True, exist_ok=False)
    local_xyz = scratch / "candidate.xyz"
    if geometry_text is None:
        shutil.copyfile(xyz, local_xyz)
    else:
        local_xyz.write_text(geometry_text)
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"
    result = subprocess.run([executable, "candidate.xyz", "--gfn", "2", "--sp", "--chrg", str(charge),
                             "--uhf", str(multiplicity - 1)],
                            cwd=scratch, env=env, text=True, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, timeout=600, check=False)
    (scratch / "xtb.log").write_text(result.stdout)
    if result.returncode:
        raise RuntimeError(f"xTB failed for {xyz}: exit {result.returncode}; see {scratch / 'xtb.log'}")
    match = ENERGY.search(result.stdout)
    if not match:
        raise RuntimeError(f"No final xTB energy found for {xyz}; see {scratch / 'xtb.log'}")
    return float(match.group(1))


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workdir", type=Path)
    parser.add_argument("--output", type=Path, help="CSV destination (default: <workdir>/ms0-rescored.csv)")
    parser.add_argument("--xtb", default="xtb")
    args = parser.parse_args(argv)
    executable = shutil.which(args.xtb)
    if not executable:
        parser.error(f"xTB executable not found: {args.xtb}")
    workdir = args.workdir.resolve()
    output = (args.output if args.output and args.output.is_absolute()
              else workdir / (args.output or Path("ms0-rescored.csv")))
    candidates = []
    for candidate in workdir.rglob("*.xyz"):
        if "inputs" in candidate.parts or "rescoring" in candidate.parts:
            continue
        if "selected" in candidate.parts and candidate.name.startswith("selected_"):
            candidates.append(candidate)
        elif candidate.name.endswith(".solvator.xyz"):
            candidates.append(candidate)
        elif candidate.name in {"crest_conformers.xyz", "crest_ensemble.xyz"}:
            candidates.append(candidate)
        elif candidate.name == "cluster.xyz" and candidate.parent.name == "grow":
            candidates.append(candidate)
    candidates.sort()
    rows = []
    for xyz in candidates:
        try:
            frames = read_xyz_frames(xyz)
            metadata = next((parent / "run.json" for parent in xyz.parents
                             if (parent / "run.json").is_file()), None)
            if metadata:
                case = json.loads(metadata.read_text())
                solute_name = case["solute"]
                charge, multiplicity = int(case["charge"]), int(case["multiplicity"])
            else:
                solute_name = None
                charge, multiplicity = 0, 1
            solute_atoms = len(read_xyz(ROOT / "fixtures" / f"{solute_name}.xyz")) if solute_name else 0
            for frame_index, (symbols, geometry_text) in enumerate(frames, start=1):
                if solute_name is None:
                    solute_name = match_solute_prefix(symbols)
                    if solute_name is None:
                        raise ValueError("Candidate atom prefix does not match a benchmark solute fixture")
                    charge, multiplicity = SOLUTE_STATES[solute_name]
                    solute_atoms = len(read_xyz(ROOT / "fixtures" / f"{solute_name}.xyz"))
                scratch = workdir / "rescoring" / xyz.relative_to(workdir).with_suffix("") / f"frame_{frame_index:04d}"
                energy = energy_for(xyz, executable, scratch, charge=charge, multiplicity=multiplicity,
                                    geometry_text=geometry_text)
                water_count = count_water_fragments(symbols, solute_atoms)
                label = str(xyz.relative_to(workdir))
                if len(frames) > 1:
                    label += f"#frame={frame_index}"
                rows.append({"path": label, "atoms": len(symbols), "water_count": water_count,
                             "charge": charge, "energy_hartree": f"{energy:.12f}", "status": "ok"})
        except (ValueError, RuntimeError, subprocess.TimeoutExpired) as exc:
            rows.append({"path": str(xyz.relative_to(workdir)), "atoms": "", "water_count": "", "charge": "",
                         "energy_hartree": "", "status": f"failed: {exc}"})
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["path", "atoms", "water_count", "charge", "energy_hartree", "status"])
        writer.writeheader()
        writer.writerows(rows)
    print(f"rescored={sum(row['status'] == 'ok' for row in rows)} failed={sum(row['status'] != 'ok' for row in rows)} output={output}")
    return 1 if not rows or any(row["status"] != "ok" for row in rows) else 0


if __name__ == "__main__":
    raise SystemExit(main())
