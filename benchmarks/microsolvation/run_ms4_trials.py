#!/usr/bin/env python3
"""Run reproducible independent GFN2-xTB PyAR cluster-growth trials for MS4 readiness."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parent
PROJECT = ROOT.parents[1]
sys.path.insert(0, str(PROJECT))

from pyar.core.molecule import Molecule  # noqa: E402
from pyar.structure_comparison import GraphRMSDComparator  # noqa: E402

SOLUTES = {"methane": (0, 1), "ammonia": (0, 1), "acetate": (-1, 1)}
TRIAL_SEEDS = (1103, 2207, 3301)
HARTREE_TO_KCAL_MOL = 627.509474


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rotate_xyz(source: Path, destination: Path, seed: int) -> None:
    lines = source.read_text().splitlines()
    count = int(lines[0])
    symbols, coordinates = [], []
    for line in lines[2:2 + count]:
        fields = line.split()
        symbols.append(fields[0])
        coordinates.append([float(value) for value in fields[1:4]])
    xyz = np.asarray(coordinates, dtype=float)
    center = np.average(xyz, axis=0, weights=np.asarray([15.999, 1.008, 1.008]))
    rng = np.random.default_rng(seed)
    axis = rng.normal(size=3)
    axis /= np.linalg.norm(axis)
    angle = float(rng.uniform(0.0, 2.0 * math.pi))
    cross = np.array([[0.0, -axis[2], axis[1]],
                      [axis[2], 0.0, -axis[0]],
                      [-axis[1], axis[0], 0.0]])
    rotation = (np.eye(3) + math.sin(angle) * cross
                + (1.0 - math.cos(angle)) * (cross @ cross))
    rotated = (xyz - center) @ rotation.T + center
    destination.write_text(
        f"{count}\nwater input orientation seed={seed}; rigid rotation\n"
        + "".join(f"{symbol:2s} {x: .10f} {y: .10f} {z: .10f}\n"
                 for symbol, (x, y, z) in zip(symbols, rotated))
    )


def write_xtb_counter(bin_dir: Path, executable: str) -> Path:
    wrapper = bin_dir / "xtb"
    wrapper.write_text(
        "#!/usr/bin/env bash\n"
        "set -u\n"
        "log=${XTB_CALL_LOG:?XTB_CALL_LOG is required}\n"
        "blocked=${XTB_BLOCKED_LOG:?XTB_BLOCKED_LOG is required}\n"
        "limit=${XTB_CALL_BUDGET:-0}\n"
        "count=$(wc -l < \"$log\")\n"
        "if (( limit > 0 && count >= limit )); then printf '%s\\n' \"$*\" >> \"$blocked\"; exit 124; fi\n"
        "printf '%s\\n' \"$*\" >> \"$log\"\n"
        f"exec {executable!r} \"$@\"\n"
    )
    wrapper.chmod(0o755)
    return wrapper


def read_xyz(path: Path):
    lines = path.read_text().splitlines()
    count = int(lines[0])
    records = [line.split() for line in lines[2:2 + count]]
    if len(records) != count or any(len(row) < 4 for row in records):
        raise ValueError(f"malformed XYZ: {path}")
    return [row[0] for row in records]


def collect_trial(run_dir: Path, solute: str, trial: int, call_log: Path, blocked_log: Path):
    state_path = run_dir / "pyar-run" / "state.json"
    if not state_path.is_file():
        return {"solute": solute, "trial": trial, "status": "failed/no-state",
                "xtb_evaluations": sum(1 for _ in call_log.open()) if call_log.exists() else 0,
                "blocked_backend_calls": sum(1 for _ in blocked_log.open()) if blocked_log.exists() else 0,
                "final_candidates": []}
    state = json.loads(state_path.read_text())
    requested_charge, requested_multiplicity = SOLUTES[solute]
    candidates = []
    for rank, ref in enumerate(state.get("final_seeds", [])):
        path = run_dir / "pyar-run" / ref["path"]
        symbols = read_xyz(path)
        fixture_symbols = read_xyz(ROOT / "fixtures" / f"{solute}.xyz")
        waters = symbols[len(fixture_symbols):]
        if (len(symbols) != len(fixture_symbols) + 12 or symbols[:len(fixture_symbols)] != fixture_symbols
                or len(waters) != 12 or any(waters[i:i + 3] != ["O", "H", "H"]
                                             for i in range(0, len(waters), 3))):
            raise ValueError(f"final candidate identity/composition mismatch: {path}")
        energy_record = ref.get("energy_record") or {}
        if (energy_record.get("unit") != "hartree"
                or abs(float(energy_record.get("value")) - float(ref["energy"])) > 1e-12):
            raise ValueError(f"final candidate lacks a consistent physical Hartree energy record: {path}")
        candidates.append({
            "rank": rank, "path": path.relative_to(run_dir).as_posix(),
            "sha256": sha256(path), "energy_hartree": float(ref["energy"]),
            "charge": requested_charge, "multiplicity": requested_multiplicity,
            "atoms": len(symbols), "waters": 4,
        })
    return {
        "solute": solute, "trial": trial, "status": state["status"],
        "charge": requested_charge, "multiplicity": requested_multiplicity,
        "composition": {"waters": 4, "atoms": len(read_xyz(ROOT / "fixtures" / f"{solute}.xyz")) + 12},
        "xtb_evaluations": sum(1 for _ in call_log.open()) if call_log.exists() else 0,
        "blocked_backend_calls": sum(1 for _ in blocked_log.open()) if blocked_log.exists() else 0,
        "completed_steps": len(state.get("completed_steps", [])),
        "step_costs": [{"step": record["step"],
                        "generated": record["selection"].get("generated_count"),
                        "prescreened": record["selection"].get("prescreened_count"),
                        "optimized": record["selection"].get("optimized_count"),
                        "optimization_failures": len(record["selection"].get("optimization_failures", []))}
                       for record in state.get("completed_steps", [])],
        "final_candidates": candidates,
    }


def cluster_basins(records, *, threshold_angstrom=0.5, near_minimum_kcal_mol=1.0):
    candidates = []
    for trial in records:
        for item in trial["final_candidates"]:
            path = trial["run_path"] / item["path"]
            molecule = Molecule.from_xyz(path)
            molecule.charge = item["charge"]
            molecule.multiplicity = item["multiplicity"]
            molecule.energy = item["energy_hartree"]
            candidates.append((item["energy_hartree"], molecule, trial["trial"], item["rank"]))
    candidates.sort(key=lambda entry: (entry[0], entry[2], entry[3]))
    comparator = GraphRMSDComparator(threshold=threshold_angstrom, atom_mode="heavy")
    basins = []
    for entry in candidates:
        if not any(comparator.compare(entry[1], representative[1]).equivalent is True
                   for representative in basins):
            basins.append(entry)
    if not candidates:
        return {"candidate_count": 0, "distinct_basins": 0, "near_minimum_basins": 0,
                "minimum_energy_hartree": None, "mean_nearest_basin_rmsd_angstrom": None}
    best = candidates[0][0]
    near_best = [entry for entry in basins
                 if (entry[0] - best) * HARTREE_TO_KCAL_MOL <= near_minimum_kcal_mol]
    distances = []
    for index, (_, molecule, *_rest) in enumerate(basins):
        nearest = [comparator.compare(molecule, other[1]).distance
                   for other_index, other in enumerate(basins) if other_index != index]
        nearest = [value for value in nearest if value is not None]
        if nearest:
            distances.append(min(nearest))
    return {
        "candidate_count": len(candidates), "distinct_basins": len(basins),
        "near_minimum_basins_within_1_kcal_mol": len(near_best),
        "minimum_energy_hartree": best,
        "basin_rmsd_threshold_angstrom": threshold_angstrom,
        "near_minimum_window_kcal_mol": near_minimum_kcal_mol,
        "mean_nearest_basin_rmsd_angstrom": float(np.mean(distances)) if distances else None,
        "basin_representatives": [{"energy_hartree": e, "trial": t, "rank": r}
                                   for e, _molecule, t, r in basins],
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workdir", required=True, type=Path)
    parser.add_argument("--pyar", default="pyar")
    parser.add_argument("--xtb", default=shutil.which("xtb") or "xtb")
    parser.add_argument("--solutes", nargs="+", choices=tuple(SOLUTES), default=list(SOLUTES))
    parser.add_argument("--replicates", type=int, default=3)
    parser.add_argument("--orientations", type=int, default=8)
    parser.add_argument("--maximum-seeds", type=int, default=4)
    parser.add_argument("--opt-cycles", type=int, default=150)
    args = parser.parse_args(argv)
    if args.replicates < 2 or args.replicates > len(TRIAL_SEEDS):
        parser.error("replicates must be between 2 and 3 for the versioned seed set")
    pyar = shutil.which(args.pyar)
    xtb = shutil.which(args.xtb)
    if not pyar or not xtb:
        parser.error(f"required executables not found: pyar={pyar}, xtb={xtb}")
    root = args.workdir.expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)
    records = []
    version = subprocess.run([xtb, "--version"], text=True, stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT, check=False).stdout
    for solute in args.solutes:
        charge, multiplicity = SOLUTES[solute]
        for trial, seed in enumerate(TRIAL_SEEDS[:args.replicates], start=1):
            run_dir = root / f"{solute}-water-n4" / f"trial_{trial:02d}"
            run_dir.mkdir(parents=True, exist_ok=False)
            inputs = run_dir / "inputs"
            inputs.mkdir()
            solute_source = ROOT / "fixtures" / f"{solute}.xyz"
            solvent_source = ROOT / "fixtures" / "water.xyz"
            solute_input = inputs / "solute.xyz"
            solvent_input = inputs / "water.xyz"
            shutil.copyfile(solute_source, solute_input)
            rotate_xyz(solvent_source, solvent_input, seed)
            bin_dir = run_dir / "bin"
            bin_dir.mkdir()
            write_xtb_counter(bin_dir, xtb)
            call_log = run_dir / "xtb_calls.log"
            blocked_log = run_dir / "xtb_calls_blocked.log"
            call_log.touch(); blocked_log.touch()
            command = [pyar, "microsolvate", str(solute_input), str(solvent_input),
                       "--count", "4", "--orientations", str(args.orientations),
                       "--maximum-number-of-seeds", str(args.maximum_seeds),
                       "--surface-points-per-atom", "32", "--objective", "cluster",
                       "--placement-score", "energy", "--backend", "xtb", "--xtb-model", "gfn2",
                       "--geometry-optimizer", "geometric", "--opt-cycles", str(args.opt_cycles),
                       "--nprocs", "1", "--charge", str(charge), "0",
                       "--multiplicity", str(multiplicity), "1",
                       "--output", str(run_dir / "pyar-run"), "--json"]
            env = os.environ.copy()
            env["PATH"] = str(bin_dir) + os.pathsep + env.get("PATH", "")
            env["PYTHONPATH"] = str(PROJECT) + os.pathsep + env.get("PYTHONPATH", "")
            env["OMP_NUM_THREADS"] = "1"
            env["XTB_CALL_LOG"] = str(call_log)
            env["XTB_BLOCKED_LOG"] = str(blocked_log)
            started = time.monotonic()
            with (run_dir / "pyar.log").open("w") as log:
                completed = subprocess.run(command, cwd=run_dir, env=env, stdout=log,
                                           stderr=subprocess.STDOUT, check=False)
            record = collect_trial(run_dir, solute, trial, call_log, blocked_log)
            record.update({"run_path": run_dir, "status": record["status"] if completed.returncode == 0
                           else f"failed/exit-{completed.returncode}",
                           "walltime_seconds": round(time.monotonic() - started, 3),
                           "command": command, "orientation_seed": seed,
                           "input_sha256": {"solute": sha256(solute_input), "water": sha256(solvent_input)},
                           "python": sys.version.split()[0], "xtb_version": version.strip()})
            (run_dir / "trial.json").write_text(json.dumps(
                {key: (str(value) if isinstance(value, Path) else value)
                 for key, value in record.items()}, indent=2) + "\n")
            records.append(record)
            print(f"{solute} trial={trial} status={record['status']} candidates="
                  f"{len(record['final_candidates'])} xtb_calls={record['xtb_evaluations']} "
                  f"wall_s={record['walltime_seconds']}", flush=True)
    summaries = {}
    for solute in args.solutes:
        case_records = [record for record in records if record["solute"] == solute]
        for record in case_records:
            record["run_path"] = Path(record["run_path"])
        summaries[solute] = cluster_basins(case_records)
    report = {"protocol": "ms4-pyar-cluster-growth-v1", "count": 4,
              "orientation_seeds": list(TRIAL_SEEDS[:args.replicates]),
              "orientations_per_step": args.orientations, "maximum_seeds": args.maximum_seeds,
              "opt_cycles": args.opt_cycles, "energy_unit": "hartree", "hamiltonian": "GFN2-xTB",
              "basin_method": "heavy-atom graph-constrained Kabsch RMSD",
              "trials": [{key: (str(value) if isinstance(value, Path) else value)
                          for key, value in record.items()} for record in records],
              "case_summaries": summaries}
    (root / "ms4_trials_summary.json").write_text(json.dumps(report, indent=2) + "\n")
    csv_path = root / "ms4_trials.csv"
    with csv_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=["solute", "trial", "status", "charge", "multiplicity",
            "completed_steps", "candidate_count", "xtb_evaluations", "blocked_backend_calls", "walltime_seconds",
            "minimum_energy_hartree", "orientation_seed"])
        writer.writeheader()
        for record in records:
            energies = [candidate["energy_hartree"] for candidate in record["final_candidates"]]
            writer.writerow({"solute": record["solute"], "trial": record["trial"], "status": record["status"],
                "charge": record.get("charge"), "multiplicity": record.get("multiplicity"),
                "completed_steps": record.get("completed_steps", 0),
                "candidate_count": len(record["final_candidates"]),
                "xtb_evaluations": record["xtb_evaluations"],
                "blocked_backend_calls": record["blocked_backend_calls"],
                "walltime_seconds": record.get("walltime_seconds"),
                "minimum_energy_hartree": min(energies) if energies else None,
                "orientation_seed": record.get("orientation_seed")})
    print(f"summary={root / 'ms4_trials_summary.json'}")
    return 1 if any(not record["status"].startswith("completed") for record in records) else 0


if __name__ == "__main__":
    raise SystemExit(main())
