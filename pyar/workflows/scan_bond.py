"""ORCA-only relaxed bond-scan workflow."""

from __future__ import annotations

import csv
import json
import math
import os
import shutil
from pathlib import Path

import numpy as np

from pyar.backends import write_xyz
from pyar.backends.orca_scan import OrcaBondScanRequest, run_orca_bond_scan
from pyar.core.molecule import Molecule
from pyar.optimiser import optimise, is_success
from pyar.sampling.trial_generator import generate_trial_vectors, merge_two_molecules

DEFAULT_SCAN_END_FACTOR = 0.8


def scan_point_count(start, end, step=None, points=None):
    if not all(math.isfinite(float(value)) for value in (start, end)) or start <= 0 or end <= 0:
        raise ValueError("scan distances must be finite and positive")
    if math.isclose(start, end, rel_tol=0.0, abs_tol=1e-10):
        raise ValueError("scan start and end distances must differ")
    if points is not None and step is not None:
        raise ValueError("--scan-step and --scan-points cannot be used together")
    if points is not None:
        if int(points) != points or points < 2:
            raise ValueError("scan points must be an integer >= 2")
        return int(points)
    if step is None:
        step = 0.10
    if not math.isfinite(float(step)) or step <= 0:
        raise ValueError("scan step must be finite and positive")
    return max(2, int(math.ceil(abs(start - end) / step)) + 1)


def absolute_target_indices(fragment_a, local_i, local_j):
    if local_i < 0 or local_i >= fragment_a.number_of_atoms:
        raise ValueError(f"Atom index {local_i} is out of range for fragment A containing {fragment_a.number_of_atoms} atoms (valid 0..{fragment_a.number_of_atoms - 1}).")
    if local_j < 0:
        raise ValueError("Atom index for fragment B must be non-negative")
    return local_i, fragment_a.number_of_atoms + local_j


def _distance(coordinates, i, j):
    return float(np.linalg.norm(np.asarray(coordinates[i]) - np.asarray(coordinates[j])))


def _fragment_signature(molecule):
    return {"atoms": molecule.atoms_list, "coordinates": np.asarray(molecule.coordinates).tolist(),
            "charge": molecule.charge, "multiplicity": molecule.multiplicity, "scftype": molecule.scftype}


def _make_orientations(fragment_a, fragment_b, n, local_i, local_j):
    vectors = generate_trial_vectors(n, direction_method="fibonacci", rotation_method="halton")
    orientations = []
    for index, vector in enumerate(vectors):
        try:
            orientation = merge_two_molecules(vector, fragment_a, fragment_b, distance_scaling=1.5)
        except Exception as exc:
            orientations.append((index, None, str(exc)))
            continue
        orientation.fragments = [list(range(fragment_a.number_of_atoms)),
                                list(range(fragment_a.number_of_atoms, orientation.number_of_atoms))]
        absolute_i = local_i
        absolute_j = fragment_a.number_of_atoms + local_j
        current = orientation.coordinates[absolute_j] - orientation.coordinates[absolute_i]
        norm = np.linalg.norm(current)
        if norm == 0.0:
            current = np.array([1.0, 0.0, 0.0])
            norm = 1.0
        desired = 1.5 * (orientation.covalent_radius[absolute_i] + orientation.covalent_radius[absolute_j])
        translated = orientation.coordinates.copy()
        if not np.isclose(norm, desired):
            translated[fragment_a.number_of_atoms:] += (
                orientation.coordinates[absolute_i] + desired * current / norm
                - orientation.coordinates[absolute_j]
            )
            valid = True
            for atom_i in range(fragment_a.number_of_atoms):
                for atom_j in range(fragment_a.number_of_atoms, orientation.number_of_atoms):
                    radii = orientation.covalent_radius[atom_i] + orientation.covalent_radius[atom_j]
                    if _distance(translated, atom_i, atom_j) < 0.8 * radii:
                        valid = False
                        break
                if not valid:
                    break
            if valid:
                orientation.coordinates = translated
        orientation.name = f"orientation_{index:03d}"
        orientation.title = orientation.name
        orientations.append((index, orientation, None))
    return orientations


def run_scan_bond(input_a, input_b, atoms, orientations, qc_params, output_dir,
                  scan_end=None, scan_step=None, scan_points=None):
    """Run the complete ORCA scan/relaxation workflow."""
    if str(qc_params.get("software", "orca")).lower() != "orca":
        raise ValueError("scan-bond currently supports only the ORCA backend")
    fragment_a = Molecule.from_xyz(input_a)
    fragment_b = Molecule.from_xyz(input_b)
    for fragment, suffix in ((fragment_a, "a"), (fragment_b, "b")):
        fragment.charge = int(qc_params.get(f"charge_{suffix}", fragment.charge))
        fragment.multiplicity = int(qc_params.get(f"multiplicity_{suffix}", fragment.multiplicity))
        fragment.scftype = qc_params.get(f"scftype_{suffix}", fragment.scftype)
    local_i, local_j = map(int, atoms)
    absolute_i, absolute_j = absolute_target_indices(fragment_a, local_i, local_j)
    if local_j >= fragment_b.number_of_atoms:
        raise ValueError(f"Atom index {local_j} is out of range for fragment B containing {fragment_b.number_of_atoms} atoms (valid 0..{fragment_b.number_of_atoms - 1}).")
    merged = fragment_a.merged_with(fragment_b)
    target_radii = merged.covalent_radius[absolute_i] + merged.covalent_radius[absolute_j]
    default_scan_end = DEFAULT_SCAN_END_FACTOR * target_radii
    root = Path(output_dir)
    if root.exists():
        raise FileExistsError(f"scan-bond output directory already exists: {root}")
    root.mkdir(parents=True)
    request = {
        "inputs": {"A": str(input_a), "B": str(input_b)},
        "fragments": [_fragment_signature(fragment_a), _fragment_signature(fragment_b)],
        "atoms_local": [local_i, local_j], "atoms_absolute": [absolute_i, absolute_j],
        "target_symbols": [merged.atoms_list[absolute_i], merged.atoms_list[absolute_j]],
        "qc_params": dict(qc_params), "orientations": orientations,
        "scan_end": scan_end, "scan_step": scan_step, "scan_points": scan_points,
        "default_scan_end_factor": DEFAULT_SCAN_END_FACTOR,
        "default_scan_end_angstrom": default_scan_end,
    }
    (root / "request.json").write_text(json.dumps(request, indent=2, sort_keys=True, default=str))
    results = []
    def save_result(directory, result):
        (directory / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True, default=str))

    for index, orientation, generation_error in _make_orientations(
            fragment_a, fragment_b, orientations, local_i, local_j):
        directory = root / f"orientation_{index:03d}"
        scan_dir = directory / "scan"
        relax_dir = directory / "relax"
        directory.mkdir()
        if orientation is None:
            result = {"orientation": index, "status": "scan_failed",
                      "scan_status": "orientation_generation_failed",
                      "error": generation_error}
            save_result(directory, result)
            results.append(result)
            continue
        write_xyz(orientation.atoms_list, orientation.coordinates, directory / "start.xyz", job_name=orientation.name, precision=12)
        start_distance = _distance(orientation.coordinates, absolute_i, absolute_j)
        end_distance = float(scan_end if scan_end is not None else default_scan_end)
        try:
            n_points = scan_point_count(start_distance, end_distance, scan_step, scan_points)
        except ValueError as exc:
            result = {"status": "scan_failed", "error": str(exc), "orientation": index}
            save_result(directory, result)
            results.append(result)
            continue
        request_model = OrcaBondScanRequest(absolute_i, absolute_j, start_distance, end_distance, n_points)
        try:
            scan = run_orca_bond_scan(orientation, request_model, scan_dir, qc_params)
        except Exception as exc:
            result = {"orientation": index, "status": "scan_failed", "scan_status": "exception",
                      "error": str(exc), "target_distance_start_angstrom": start_distance,
                      "target_distance_scan_end_angstrom": end_distance, "scan_points": n_points}
            save_result(directory, result)
            results.append(result)
            continue
        result = {"orientation": index, "status": "scan_failed", "scan_status": scan.status,
                  "target_distance_start_angstrom": start_distance, "target_distance_scan_end_angstrom": end_distance,
                  "scan_points": n_points, "scan_input": str(scan.input_path), "scan_output": str(scan.output_path)}
        if not scan.success:
            save_result(directory, result)
            results.append(result)
            continue
        result.update({"scan_status": "success", "trajectory_path": str(scan.trajectory_path),
                       "final_scan_path": str(scan.final_geometry_path),
                       "target_distance_final_scan_angstrom": _distance(scan.final_coordinates, absolute_i, absolute_j)})
        relaxed = Molecule(orientation.atoms_list, scan.final_coordinates, name="relaxed",
                           fragments=orientation.fragments, charge=orientation.charge,
                           multiplicity=orientation.multiplicity, scftype=orientation.scftype)
        relax_dir.mkdir(parents=True, exist_ok=True)
        old_cwd = Path.cwd()
        try:
            os.chdir(relax_dir)
            relax_status = optimise(relaxed, dict(qc_params, geometry_optimizer="native", gamma=0.0))
        except Exception as exc:
            result["status"] = "scan_success_relax_failed"
            result["error"] = str(exc)
            save_result(directory, result)
            results.append(result)
            continue
        finally:
            os.chdir(old_cwd)
        if is_success(relax_status):
            result["status"] = "scan_success_relax_success"
            result["relaxed_path"] = str(directory / "result_relaxed.xyz")
            shutil.copy(relax_dir / "job_relaxed" / "result_relaxed.xyz", directory / "result_relaxed.xyz")
            result["target_distance_relaxed_angstrom"] = _distance(relaxed.coordinates, absolute_i, absolute_j)
            result["target_bond_present_after_relaxation"] = result["target_distance_relaxed_angstrom"] < 1.3 * target_radii
        else:
            result["status"] = "scan_success_relax_failed"
        save_result(directory, result)
        results.append(result)
    (root / "summary.json").write_text(json.dumps(results, indent=2, sort_keys=True, default=str))
    with (root / "summary.csv").open("w", newline="") as handle:
        fields = sorted({key for result in results for key in result})
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader(); writer.writerows(results)
    return {"workflow": "scan-bond", "output_dir": str(root), "results": results}
