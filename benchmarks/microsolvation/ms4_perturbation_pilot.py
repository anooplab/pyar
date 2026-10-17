#!/usr/bin/env python3
"""External, non-production whole-cluster perturbation pilot for MS4 readiness."""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import shutil
import sys

import numpy as np
from ase import Atoms

ROOT = Path(__file__).resolve().parent
PROJECT = ROOT.parents[1]
sys.path.insert(0, str(PROJECT))

from pyar.core.molecule import Molecule  # noqa: E402
from pyar.optimiser import is_usable, optimise  # noqa: E402
from pyar.structure_comparison import GraphRMSDComparator  # noqa: E402

HARTREE_TO_KCAL_MOL = 627.509474


def rigid_move_solvent(molecule, fragment_index: int, seed: int,
                       translation_range=(1.0, 2.0)) -> dict:
    """Rotate and translate one solvent fragment while preserving its geometry."""
    fragment = np.asarray(molecule.solvent_fragments[fragment_index], dtype=int)
    xyz = np.asarray(molecule.coordinates, dtype=float)
    masses = Atoms(symbols=molecule.atoms_list, positions=xyz).get_masses()
    center = np.average(xyz[fragment], axis=0, weights=masses[fragment])
    rng = np.random.default_rng(seed)
    axis = rng.normal(size=3)
    axis /= np.linalg.norm(axis)
    angle = float(rng.uniform(math.pi / 3.0, math.pi))
    cross = np.array([[0.0, -axis[2], axis[1]],
                      [axis[2], 0.0, -axis[0]],
                      [-axis[1], axis[0], 0.0]])
    rotation = np.eye(3) + math.sin(angle) * cross + (1.0 - math.cos(angle)) * (cross @ cross)
    direction = rng.normal(size=3)
    direction /= np.linalg.norm(direction)
    distance = float(rng.uniform(*translation_range))
    translation = direction * distance
    molecule.coordinates[fragment] = (xyz[fragment] - center) @ rotation.T + center + translation
    return {"fragment_index": fragment_index, "rotation_axis": axis.tolist(),
            "rotation_angle_radians": angle, "translation_angstrom": translation.tolist(),
            "translation_distance_angstrom": distance, "proposal_seed": seed}


def collision_free(molecule, minimum_distance_angstrom: float) -> tuple[bool, float]:
    xyz = np.asarray(molecule.coordinates, dtype=float)
    groups = [np.asarray(molecule.solute_atom_indices, dtype=int),
              *[np.asarray(fragment, dtype=int) for fragment in molecule.solvent_fragments]]
    minimum = math.inf
    for left in range(len(groups)):
        for right in range(left + 1, len(groups)):
            delta = xyz[groups[left], None, :] - xyz[groups[right]][None, :, :]
            minimum = min(minimum, float(np.min(np.linalg.norm(delta, axis=2))))
    return minimum >= minimum_distance_angstrom, minimum


def molecule_from_start(path: Path, metadata: dict, name: str) -> Molecule:
    molecule = Molecule.from_xyz(path)
    molecule.name = name
    molecule.charge = int(metadata["charge"])
    molecule.multiplicity = int(metadata["multiplicity"])
    if metadata.get("energy_hartree") is not None:
        molecule.energy = float(metadata["energy_hartree"])
    molecule.scftype = "rhf" if molecule.multiplicity == 1 else "uhf"
    molecule.solute_atom_indices = tuple(metadata["solute_atom_indices"])
    molecule.solvent_fragments = tuple(tuple(part) for part in metadata["solvent_fragments"])
    molecule.fragments = [list(molecule.solute_atom_indices), *[list(part) for part in molecule.solvent_fragments]]
    molecule.fragments_coordinates = molecule.split_coordinates()
    molecule.fragments_atoms_list = molecule.split_atoms_lists()
    return molecule


def count_lines(path: Path) -> int:
    return sum(1 for _ in path.open()) if path.exists() else 0


def archive_metrics(molecules, threshold=0.5, near_minimum_kcal=1.0):
    molecules = sorted(molecules, key=lambda item: (float(item.energy), item.name))
    comparator = GraphRMSDComparator(threshold=threshold, atom_mode="heavy")
    representatives = []
    for molecule in molecules:
        if not any(comparator.compare(molecule, prior).equivalent is True for prior in representatives):
            representatives.append(molecule)
    if not molecules:
        return {"candidate_count": 0, "distinct_basins": 0, "near_minimum_basins": 0,
                "minimum_energy_hartree": None}
    minimum = float(molecules[0].energy)
    near = [molecule for molecule in representatives
            if (float(molecule.energy) - minimum) * HARTREE_TO_KCAL_MOL <= near_minimum_kcal]
    return {"candidate_count": len(molecules), "distinct_basins": len(representatives),
            "near_minimum_basins_within_1_kcal_mol": len(near),
            "minimum_energy_hartree": minimum,
            "basin_rmsd_threshold_angstrom": threshold,
            "near_minimum_window_kcal_mol": near_minimum_kcal,
            "basin_representatives": [{"name": item.name, "energy_hartree": float(item.energy)}
                                      for item in representatives]}


def count_novel_basins(candidates, starts, threshold=0.5):
    comparator = GraphRMSDComparator(threshold=threshold, atom_mode="heavy")
    initial_representatives = []
    for molecule in sorted(starts, key=lambda item: float(item.energy)):
        if not any(comparator.compare(molecule, prior).equivalent is True
                   for prior in initial_representatives):
            initial_representatives.append(molecule)
    novel = []
    for molecule in sorted(candidates, key=lambda item: float(item.energy)):
        if any(comparator.compare(molecule, prior).equivalent is True
               for prior in initial_representatives):
            continue
        if any(comparator.compare(molecule, prior).equivalent is True for prior in novel):
            continue
        novel.append(molecule)
    return len(initial_representatives), len(novel)


def load_start_pool(growth_summary: Path, solute: str):
    summary = json.loads(growth_summary.read_text())
    records = [item for item in summary["trials"]
               if item["solute"] == solute and item["status"] == "completed"]
    starts = []
    for trial in records:
        run_dir = Path(trial["run_path"])
        if not run_dir.is_absolute():
            run_dir = growth_summary.parent / run_dir
        pool = trial["final_candidates"]
        if not pool:
            continue
        best = min(pool, key=lambda item: item["energy_hartree"])
        ref = {"charge": best["charge"], "multiplicity": best["multiplicity"],
               "energy_hartree": float(best["energy_hartree"]),
               "solute_atom_indices": list(range(len(Molecule.from_xyz(
                   ROOT / "fixtures" / f"{solute}.xyz")))),
               "solvent_fragments": [list(range(len(Molecule.from_xyz(
                   ROOT / "fixtures" / f"{solute}.xyz")) + 3 * index,
                   len(Molecule.from_xyz(ROOT / "fixtures" / f"{solute}.xyz")) + 3 * index + 3))
                   for index in range(4)]}
        starts.append((run_dir / best["path"], ref, int(trial["trial"]), float(best["energy_hartree"])))
    if not starts:
        raise ValueError(f"No completed growth structures for solute {solute}")
    return starts


def run_arm(mode: str, starts, workdir: Path, xtb_wrapper: Path, budget: int,
            minimum_distance: float, opt_cycles: int):
    workdir.mkdir(parents=True, exist_ok=False)
    calls = workdir / "xtb_calls.log"
    blocked = workdir / "xtb_calls_blocked.log"
    calls.touch(); blocked.touch()
    os.environ["XTB_CALL_LOG"] = str(calls)
    os.environ["XTB_BLOCKED_LOG"] = str(blocked)
    os.environ["XTB_CALL_BUDGET"] = str(budget)
    os.environ["PATH"] = str(xtb_wrapper.parent) + os.pathsep + os.environ["PATH"]
    candidates = []
    attempts = []
    attempt_index = 0
    maximum_proposals = budget * 12
    start_molecules = [molecule_from_start(path, metadata, f"start_trial_{trial:02d}")
                       for path, metadata, trial, _energy in starts]
    while count_lines(calls) < budget and attempt_index < maximum_proposals:
        source_path, source_metadata, source_trial, source_energy = starts[attempt_index % len(starts)]
        molecule = molecule_from_start(source_path, source_metadata, f"{mode}_{attempt_index:04d}")
        molecule.name = f"{mode}_{attempt_index:04d}"
        proposal = None
        min_distance = None
        if mode == "perturb":
            fragment_index = int(np.random.default_rng(901001 + attempt_index).integers(
                0, len(molecule.solvent_fragments)))
            proposal = rigid_move_solvent(molecule, fragment_index, 712003 + attempt_index)
            free, min_distance = collision_free(molecule, minimum_distance)
            if not free:
                attempts.append({"attempt": attempt_index, "status": "rejected_collision",
                                 "source_trial": source_trial, "source_energy_hartree": source_energy,
                                 "minimum_interfragment_distance_angstrom": min_distance,
                                 "proposal": proposal})
                attempt_index += 1
                continue
        else:
            proposal = {"kind": "unperturbed-reoptimization"}
        qc = {"software": "xtb", "xtb_model": "gfn2", "geometry_optimizer": "geometric",
              "opt_cycles": opt_cycles, "opt_threshold": "normal", "nprocs": 1,
              "charge": molecule.charge, "multiplicity": molecule.multiplicity,
              "scftype": molecule.scftype, "opt_target": "minimum", "gamma": None}
        attempt_dir = workdir / f"attempt_{attempt_index:04d}"
        attempt_dir.mkdir()
        previous = Path.cwd()
        try:
            os.chdir(attempt_dir)
            status = optimise(molecule, qc)
        except Exception as exc:  # Record experimental backend failures and continue only within the call cap.
            status = False
            error = str(exc)
        else:
            error = None
        finally:
            os.chdir(previous)
        usable = is_usable(status) and molecule.coordinates is not None and molecule.energy is not None
        attempt = {"attempt": attempt_index, "status": "completed" if usable else "failed_optimization",
                   "source_trial": source_trial, "source_energy_hartree": source_energy,
                   "energy_hartree": None if not usable else float(molecule.energy),
                   "backend_calls_total": count_lines(calls), "proposal": proposal,
                   "minimum_interfragment_distance_angstrom": min_distance, "error": error}
        attempts.append(attempt)
        if usable:
            molecule.mol_to_xyz(str(attempt_dir / "final.xyz"))
            candidates.append(molecule.copy())
        attempt_index += 1
    initial_basins, novel_basins = count_novel_basins(candidates, start_molecules)
    result = {"mode": mode, "budget": budget, "backend_calls": count_lines(calls),
              "blocked_backend_calls": count_lines(blocked), "attempt_count": len(attempts),
              "collision_rejections": sum(attempt["status"] == "rejected_collision" for attempt in attempts),
              "failed_optimizations": sum(attempt["status"] == "failed_optimization" for attempt in attempts),
              "initial_basin_count": initial_basins, "novel_basin_count": novel_basins,
              "attempts": attempts, "metrics": archive_metrics([*start_molecules, *candidates]),
              "status": "budget_reached" if count_lines(calls) == budget else "budget_underfilled"}
    (workdir / "arm.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--growth-summary", type=Path, nargs="+", required=True,
                        help="one combined summary, or one summary per selected solute")
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--xtb", default="/usr/bin/xtb")
    parser.add_argument("--solutes", nargs="+", choices=("methane", "ammonia", "acetate"),
                        default=["methane", "ammonia", "acetate"])
    parser.add_argument("--budget", type=int, default=40, help="actual GFN2-xTB gradient calls per arm and case")
    parser.add_argument("--minimum-interfragment-distance", type=float, default=1.3)
    parser.add_argument("--opt-cycles", type=int, default=250)
    args = parser.parse_args(argv)
    if args.budget < 1 or args.minimum_interfragment_distance <= 0 or args.opt_cycles < 1:
        parser.error("budget, collision threshold and optimizer cycles must be positive")
    if len(args.growth_summary) not in (1, len(args.solutes)):
        parser.error("provide one combined growth summary or one summary per selected solute")
    root = args.workdir.expanduser().resolve()
    root.mkdir(parents=True, exist_ok=False)
    wrapper_dir = root / "bin"
    wrapper_dir.mkdir()
    script = ROOT / "run_ms4_trials.py"
    # Reuse the benchmark's transparent xTB invocation counter.
    sys.path.insert(0, str(ROOT))
    from run_ms4_trials import write_xtb_counter
    xtb = shutil.which(args.xtb)
    if not xtb:
        parser.error(f"xTB executable not found: {args.xtb}")
    wrapper = write_xtb_counter(wrapper_dir, xtb)
    all_results = {}
    for index, solute in enumerate(args.solutes):
        summary_index = 0 if len(args.growth_summary) == 1 else index
        starts = load_start_pool(args.growth_summary[summary_index].resolve(), solute)
        case = root / solute
        case.mkdir()
        arms = {}
        for mode in ("baseline", "perturb"):
            arms[mode] = run_arm(mode, starts, case / mode, wrapper, args.budget,
                                 args.minimum_interfragment_distance, args.opt_cycles)
        all_results[solute] = arms
        print(f"{solute}: baseline={arms['baseline']['backend_calls']} calls, "
              f"perturb={arms['perturb']['backend_calls']} calls, "
              f"distinct basins={arms['baseline']['metrics']['distinct_basins']} -> "
              f"{arms['perturb']['metrics']['distinct_basins']}", flush=True)
    output = {"protocol": "ms4-external-rigid-fragment-perturbation-v1",
              "budget_per_arm_per_case": args.budget,
              "collision_cutoff_angstrom": args.minimum_interfragment_distance,
              "optimization_cycles": args.opt_cycles,
              "selection": "lowest energy output per arm; no Monte Carlo acceptance",
              "cases": all_results}
    (root / "ms4_perturbation_summary.json").write_text(json.dumps(output, indent=2) + "\n")
    return 1 if any(arm["backend_calls"] != args.budget
                    for case in all_results.values() for arm in case.values()) else 0


if __name__ == "__main__":
    raise SystemExit(main())
