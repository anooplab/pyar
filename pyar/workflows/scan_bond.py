"""ORCA-only relaxed bond-scan workflow."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import shutil
from pathlib import Path

import numpy as np

from pyar.backends import write_xyz
from pyar.backends.orca_scan import (
    OrcaBondScanRequest,
    load_orca_bond_scan_result,
    run_orca_bond_scan,
)
from pyar.core.molecule import Molecule
from pyar.optimiser import optimise, is_success
from pyar.reaction_identity import (
    molecule_identity_from_xyz,
    reaction_product_changed,
    separated_reactant_identity,
)
from pyar.sampling.trial_generator import generate_trial_vectors, merge_two_molecules

DEFAULT_SCAN_END_FACTOR = 0.8
TARGET_CONTACT_FACTOR = 1.3


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


def _write_scan_profile(directory, frames, profile):
    """Write a normalized energy profile and scan-maximum geometry candidates."""
    if len(frames) != len(profile) or not profile:
        raise ValueError("scan geometry and energy point counts do not match")
    first_energy = profile[0]["energy_hartree"]
    profile_rows = []
    for frame_index, point in enumerate(profile):
        profile_rows.append({
            "scan_index": frame_index + 1,
            "frame_index": frame_index,
            "target_distance_angstrom": point["target_distance_angstrom"],
            "energy_hartree": point["energy_hartree"],
            "relative_energy_kcal_mol": (point["energy_hartree"] - first_energy) * 627.509474,
            "geometry_path": "scan_trajectory.xyz",
        })
    csv_path = directory / "scan_profile.csv"
    with csv_path.open("w", newline="") as handle:
        fields = list(profile_rows[0])
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(profile_rows)
    json_path = directory / "scan_profile.json"
    json_path.write_text(json.dumps({
        "energy_source": "scan.relaxscanact.dat",
        "relative_energy_reference_scan_index": 1,
        "points": profile_rows,
    }, indent=2, sort_keys=True))

    maximum_frame = max(range(len(profile)), key=lambda index: profile[index]["energy_hartree"])
    candidate_dir = directory.parent / "ts_candidates"
    candidate_dir.mkdir(parents=True, exist_ok=True)
    candidate_files = {}
    for label, frame_index in (
        ("pre_maximum", maximum_frame - 1),
        ("highest_scan_energy", maximum_frame),
        ("post_maximum", maximum_frame + 1),
    ):
        if 0 <= frame_index < len(frames):
            symbols, coordinates, _ = frames[frame_index]
            candidate_path = candidate_dir / f"{label}.xyz"
            write_xyz(symbols, coordinates, candidate_path,
                      job_name=f"{label}_scan_{frame_index + 1:03d}",
                      energy=profile[frame_index]["energy_hartree"], precision=12)
            candidate_files[label] = str(candidate_path)
    maximum_energy = profile[maximum_frame]["energy_hartree"]
    maximum_metadata = {
        "energy_source": "scan.relaxscanact.dat",
        "maximum_scan_index": maximum_frame + 1,
        "maximum_frame_index": maximum_frame,
        "maximum_target_distance_angstrom": profile[maximum_frame]["target_distance_angstrom"],
        "maximum_energy_hartree": maximum_energy,
        "barrier_from_first_scan_point_kcal_mol": (maximum_energy - first_energy) * 627.509474,
        "maximum_is_internal": 0 < maximum_frame < len(profile) - 1,
        "candidate_files": candidate_files,
        "interpretation": "scan maximum candidate; not a confirmed transition state",
    }
    metadata_path = candidate_dir / "metadata.json"
    metadata_path.write_text(json.dumps(maximum_metadata, indent=2, sort_keys=True))
    return {
        "scan_profile_csv": str(csv_path),
        "scan_profile_json": str(json_path),
        "scan_energy_source": "scan.relaxscanact.dat",
        **maximum_metadata,
        "ts_candidate_metadata": str(metadata_path),
    }


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
    if isinstance(orientations, bool) or int(orientations) != orientations or orientations < 1:
        raise ValueError("orientation count must be a positive integer")
    if scan_end is not None and (not math.isfinite(float(scan_end)) or scan_end <= 0):
        raise ValueError("scan end distance must be finite and positive")
    if scan_step is not None and (not math.isfinite(float(scan_step)) or scan_step <= 0):
        raise ValueError("scan step must be finite and positive")
    if scan_points is not None and (int(scan_points) != scan_points or scan_points < 2):
        raise ValueError("scan points must be an integer >= 2")
    if scan_step is not None and scan_points is not None:
        raise ValueError("--scan-step and --scan-points cannot be used together")
    fragment_a = Molecule.from_xyz(input_a)
    fragment_b = Molecule.from_xyz(input_b)
    for fragment, suffix in ((fragment_a, "a"), (fragment_b, "b")):
        fragment.charge = int(qc_params.get(f"charge_{suffix}", fragment.charge))
        fragment.multiplicity = int(qc_params.get(f"multiplicity_{suffix}", fragment.multiplicity))
        fragment.scftype = qc_params.get(f"scftype_{suffix}", fragment.scftype)
    if len(atoms) != 2 or any(isinstance(value, bool) or int(value) != value for value in atoms):
        raise ValueError("exactly two integer fragment-local atom indices are required")
    local_i, local_j = map(int, atoms)
    absolute_i, absolute_j = absolute_target_indices(fragment_a, local_i, local_j)
    if local_j >= fragment_b.number_of_atoms:
        raise ValueError(f"Atom index {local_j} is out of range for fragment B containing {fragment_b.number_of_atoms} atoms (valid 0..{fragment_b.number_of_atoms - 1}).")
    merged = fragment_a.merged_with(fragment_b)
    target_radii = merged.covalent_radius[absolute_i] + merged.covalent_radius[absolute_j]
    default_scan_end = DEFAULT_SCAN_END_FACTOR * target_radii
    try:
        reactant_identity = separated_reactant_identity(fragment_a, fragment_b)
        reactant_identity_error = None
    except Exception as exc:
        reactant_identity = None
        reactant_identity_error = str(exc)
    orientation_items = _make_orientations(
        fragment_a, fragment_b, int(orientations), local_i, local_j
    )
    orientation_definitions = [
        {"orientation": index,
         "molecule": None if orientation is None else _fragment_signature(orientation),
         "generation_error": generation_error}
        for index, orientation, generation_error in orientation_items
    ]
    root = Path(output_dir)
    request = {
        "schema_version": 2,
        "inputs": {"A": str(input_a), "B": str(input_b)},
        "fragments": [_fragment_signature(fragment_a), _fragment_signature(fragment_b)],
        "atoms_local": [local_i, local_j], "atoms_absolute": [absolute_i, absolute_j],
        "target_symbols": [merged.atoms_list[absolute_i], merged.atoms_list[absolute_j]],
        "qc_params": dict(qc_params), "orientations": int(orientations),
        "orientation_definitions": orientation_definitions,
        "scan_end": scan_end, "scan_step": scan_step, "scan_points": scan_points,
        "default_scan_end_factor": DEFAULT_SCAN_END_FACTOR,
        "default_scan_end_angstrom": default_scan_end,
        "reactant_identity": reactant_identity,
        "reactant_identity_error": reactant_identity_error,
    }
    signature_payload = json.dumps(request, sort_keys=True, separators=(",", ":"), default=str)
    request["request_signature"] = hashlib.sha256(signature_payload.encode("utf-8")).hexdigest()
    if root.exists():
        request_path = root / "request.json"
        if not root.is_dir() or not request_path.is_file():
            raise FileExistsError(f"scan-bond output directory has no resumable request state: {root}")
        try:
            previous_request = json.loads(request_path.read_text())
        except (OSError, json.JSONDecodeError) as exc:
            raise FileExistsError(f"scan-bond request state is unreadable: {request_path}") from exc
        if previous_request.get("request_signature") != request["request_signature"]:
            raise FileExistsError(f"scan-bond output directory belongs to a different request: {root}")
    else:
        root.mkdir(parents=True)
    request_path = root / "request.json"
    request_temp = root / "request.json.tmp"
    request_temp.write_text(json.dumps(request, indent=2, sort_keys=True, default=str))
    request_temp.replace(request_path)
    results = []
    def save_result(directory, result):
        (directory / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True, default=str))

    for index, orientation, generation_error in orientation_items:
        directory = root / f"orientation_{index:03d}"
        scan_dir = directory / "scan"
        relax_dir = directory / "relax"
        directory.mkdir(exist_ok=True)
        result_path = directory / "result.json"
        previous_result = None
        if result_path.is_file():
            try:
                previous_result = json.loads(result_path.read_text())
            except (OSError, json.JSONDecodeError):
                previous_result = None
        completed_artifacts = (
            directory / "result_relaxed.xyz",
            scan_dir / "final_scan.xyz",
            scan_dir / "scan_trajectory.xyz",
            scan_dir / "scan_profile.csv",
            scan_dir / "scan_profile.json",
            directory / "ts_candidates" / "metadata.json",
            directory / "ts_candidates" / "highest_scan_energy.xyz",
        )
        if (previous_result
                and previous_result.get("status") == "scan_success_relax_success"
                and all(path.is_file() for path in completed_artifacts)):
            results.append(previous_result)
            continue
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
        scan = None
        if previous_result and previous_result.get("scan_status") == "success":
            scan = load_orca_bond_scan_result(orientation, request_model, scan_dir)
        if scan is None or not scan.success:
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
        try:
            if len(scan.frames or []) != n_points or len(scan.profile or []) != n_points:
                raise ValueError(
                    f"scan returned {len(scan.frames or [])} geometries and "
                    f"{len(scan.profile or [])} energies; expected {n_points} points"
                )
            result.update(_write_scan_profile(scan_dir, scan.frames, scan.profile))
        except Exception as exc:
            result["status"] = "scan_success_analysis_failed"
            result["error"] = str(exc)
            save_result(directory, result)
            results.append(result)
            continue
        relaxed = Molecule(orientation.atoms_list, scan.final_coordinates, name="relaxed",
                           fragments=orientation.fragments, charge=orientation.charge,
                           multiplicity=orientation.multiplicity, scftype=orientation.scftype)
        relax_dir.mkdir(parents=True, exist_ok=True)
        relax_run_dir = relax_dir
        if (relax_dir / "job_relaxed").exists():
            retry_index = 1
            while (relax_dir / f"retry_{retry_index:03d}").exists():
                retry_index += 1
            relax_run_dir = relax_dir / f"retry_{retry_index:03d}"
            relax_run_dir.mkdir()
        old_cwd = Path.cwd()
        try:
            os.chdir(relax_run_dir)
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
            optimized_output = relax_run_dir / "job_relaxed" / "result_relaxed.xyz"
            if not optimized_output.is_file():
                result["status"] = "scan_success_relax_failed"
                result["error"] = f"optimizer reported success but did not write {optimized_output}"
                save_result(directory, result)
                results.append(result)
                continue
            try:
                shutil.copy(optimized_output, directory / "result_relaxed.xyz")
            except OSError as exc:
                result["status"] = "scan_success_relax_failed"
                result["error"] = f"could not preserve relaxed geometry: {exc}"
                save_result(directory, result)
                results.append(result)
                continue
            result["status"] = "scan_success_relax_success"
            result["relaxed_path"] = str(directory / "result_relaxed.xyz")
            result["target_distance_relaxed_angstrom"] = _distance(relaxed.coordinates, absolute_i, absolute_j)
            result["target_contact_threshold_angstrom"] = TARGET_CONTACT_FACTOR * target_radii
            result["target_bond_present_after_relaxation"] = (
                result["target_distance_relaxed_angstrom"]
                < result["target_contact_threshold_angstrom"]
            )
            result["reactant_identity"] = reactant_identity
            result["product_identity_changed"] = None
            if reactant_identity is None:
                result["identity_status"] = "reactant_identity_failed"
                result["identity_error"] = reactant_identity_error
            else:
                try:
                    relaxed_identity = molecule_identity_from_xyz(directory / "result_relaxed.xyz")
                    result["relaxed_identity"] = relaxed_identity
                    result["product_identity_changed"] = reaction_product_changed(
                        reactant_identity, relaxed_identity
                    )
                    result["identity_status"] = "success"
                except Exception as exc:
                    result["identity_status"] = "relaxed_identity_failed"
                    result["identity_error"] = str(exc)
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
