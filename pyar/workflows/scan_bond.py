"""Relaxed bond scans and optional validated reaction-path continuation."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import warnings
import shutil
from pathlib import Path

import numpy as np

from pyar.backends import write_xyz
from pyar.backend_errors import BackendExecutionError
from pyar.backends.orca_scan import (
    OrcaBondScanRequest,
    load_orca_bond_scan_result,
    run_orca_bond_scan,
)
from pyar.core.molecule import Molecule
from pyar.optimiser import optimise, is_success
from pyar.optimization_service import optimize_workflow_once
from pyar.reaction_identity import (
    molecule_identity_from_xyz,
    reaction_product_changed,
    separated_reactant_identity,
)
from pyar.sampling.trial_generator import generate_trial_vectors, merge_two_molecules

DEFAULT_SCAN_END_FACTOR = 0.8
TARGET_CONTACT_FACTOR = 1.3
INTRAMOLECULAR_FORM_FACTOR = 1.05
INTRAMOLECULAR_BREAK_FACTOR = 2.5
INTRAMOLECULAR_SAFE_FACTOR = 2.0


def intramolecular_scan_plan(molecule, atoms, scan_end=None, scan_step=None, scan_points=None):
    """Plan a finite scan using the existing distance-based connectivity rule."""
    if not np.isfinite(molecule.coordinates).all():
        raise ValueError('Intramolecular scan coordinates must all be finite')
    if (len(atoms) != 2 or any(isinstance(i, bool) or not isinstance(i, (int, np.integer))
                              or not 0 <= i < molecule.number_of_atoms for i in atoms)
            or atoms[0] == atoms[1]):
        raise ValueError('Intramolecular scan requires two distinct 0-based atom indices within the molecule')
    i, j = map(int, atoms)
    start = _distance(molecule.coordinates, i, j)
    if not math.isfinite(start) or start <= 0:
        raise ValueError('Target atoms must have a finite nonzero separation')
    radii = float(molecule.covalent_radius[i] + molecule.covalent_radius[j])
    bonded = start < TARGET_CONTACT_FACTOR * radii
    proposed = max(INTRAMOLECULAR_BREAK_FACTOR * radii, start + 1.0) if bonded else INTRAMOLECULAR_FORM_FACTOR * radii
    end = float(proposed if scan_end is None else scan_end)
    if (not math.isfinite(end) or end <= 0 or
            (bonded and end <= start) or (not bonded and end >= start)):
        raise ValueError('Intramolecular --scan-end must elongate a bonded pair or shorten an unbonded pair')
    count = scan_point_count(start, end, scan_step, scan_points)
    if count > 10001:
        raise ValueError('Intramolecular scan exceeds 10001 points; choose a larger --scan-step or fewer --scan-points')
    return dict(scan_mode='intramolecular', scan_direction='break' if bonded else 'form',
                atoms=list(map(int, atoms)), initially_bonded=bool(bonded),
                start_distance_angstrom=start, proposed_end_distance_angstrom=proposed,
                end_distance_angstrom=end, n_points=count, intervals=count - 1,
                step_angstrom=abs(end - start) / (count - 1),
                bond_threshold_angstrom=TARGET_CONTACT_FACTOR * radii,
                safe_separation_threshold_angstrom=INTRAMOLECULAR_SAFE_FACTOR * radii,
                endpoint_source='automatic' if scan_end is None else 'user_explicit',
                connectivity_rule='distance < 1.3 * sum of covalent radii')


def _optimize_with_diagnostics(molecule, qc_params, workflow_role):
    """Run one legacy scan optimization and persist observational ORE data."""
    status, _diagnostic, _report = optimize_workflow_once(
        molecule, qc_params, workflow_role=workflow_role, optimizer=optimise
    )
    return status


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


def _save_orientation_result(directory, result, *, qc_params, intramolecular,
                             completed_artifacts, layout):
    """Persist the derived orientation report and update its layout record."""
    result.setdefault("total_multiplicity", qc_params["total_multiplicity"])
    result.setdefault("total_multiplicity_source", qc_params["total_multiplicity_source"])
    if intramolecular and result.get('status') == 'scan_success_relax_success':
        result['discovery_status'] = ('new_product' if result.get('product_identity_changed') is True else
                                      'no_new_product' if result.get('product_identity_changed') is False else
                                      'identity_unresolved')
        result.setdefault('reactant_product_connection_confirmed', False)
    if result.get("status") == "scan_success_relax_success":
        result["completed_artifact_hashes"] = {
            str(path.relative_to(directory)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in completed_artifacts if path.is_file()
        }
        if layout is not None:
            layout.relaxed(result['orientation'], result, completed_artifacts)
    (directory / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True, default=str))


def _load_or_execute_scan(orientation, request, scan_dir, qc_params, *, software,
                          previous_result, committed):
    """Reuse a valid scan or execute once, preserving committed-scan guards."""
    from pyar.backends.bond_scan import load_bond_scan_result, run_bond_scan
    from pyar.recovery.stages import scan_diagnostic

    scan = None
    if previous_result and previous_result.get("scan_status") == "success":
        if software == "orca":
            loader = load_orca_bond_scan_result
            loader_options = {"publish": False} if committed is not None else {}
        else:
            loader = load_bond_scan_result
            loader_options = {}
        scan = loader(orientation, request, scan_dir, **loader_options)
    if scan is not None and scan.success:
        if scan_diagnostic(orientation, request, scan)["stage_gate"] != "passed":
            scan.success = False
            scan.status = "scan_validation_failed"
    if committed is not None and (scan is None or not scan.success):
        raise ValueError("Committed scan no longer passes native scan validation")
    if scan is None or not scan.success:
        runner = run_orca_bond_scan if software == "orca" else run_bond_scan
        try:
            from pyar.recovery.stages import run_scan_once
            scan = run_scan_once(runner, orientation, request, scan_dir, qc_params)
        except BackendExecutionError:
            raise
        except Exception as exc:
            return None, exc
    return scan, None


def _fragment_signature(molecule):
    return {"atoms": molecule.atoms_list, "coordinates": np.asarray(molecule.coordinates).tolist(),
            "charge": molecule.charge, "multiplicity": molecule.multiplicity, "scftype": molecule.scftype}



def scan_request_fields(input_a, input_b, fragment_a, fragment_b, atoms, orientations,
                        qc_params, scan_end=None, scan_step=None, scan_points=None):
    """Scientific request fields that can be validated without sampling or jobs."""
    if fragment_b is None:
        if orientations != 1:
            raise ValueError('Intramolecular scans use the input geometry and require one orientation')
        plan = intramolecular_scan_plan(fragment_a, atoms, scan_end, scan_step, scan_points)
        total = qc_params.get('total_multiplicity')
        if total is not None and total != fragment_a.multiplicity:
            raise ValueError('Intramolecular total multiplicity must equal the input molecule multiplicity')
        return dict(schema_version=3, scan_mode='intramolecular',
                    inputs={'molecule': str(input_a)}, fragments=[_fragment_signature(fragment_a)],
                    atoms_local=list(map(int, atoms)), atoms_absolute=list(map(int, atoms)),
                    target_symbols=[fragment_a.atoms_list[int(i)] for i in atoms],
                    qc_params=dict(qc_params, total_multiplicity=fragment_a.multiplicity,
                                   total_multiplicity_source='input_molecule'), orientations=1,
                    scan_end=scan_end, scan_step=scan_step, scan_points=scan_points,
                    intramolecular_plan=plan)
    local_i, local_j = map(int, atoms)
    total_multiplicity, spin_source = _resolve_total_multiplicity(
        fragment_a, fragment_b, qc_params.get("total_multiplicity"), warn=False
    )
    spin_source = qc_params.get("total_multiplicity_source", spin_source)
    request_qc_params = dict(qc_params, total_multiplicity=total_multiplicity,
                             total_multiplicity_source=spin_source)
    absolute_i, absolute_j = absolute_target_indices(fragment_a, local_i, local_j)
    target_radii = fragment_a.covalent_radius[local_i] + fragment_b.covalent_radius[local_j]
    return {
        "schema_version": 2,
        "inputs": {"A": str(input_a), "B": str(input_b)},
        "fragments": [_fragment_signature(fragment_a), _fragment_signature(fragment_b)],
        "atoms_local": [local_i, local_j], "atoms_absolute": [absolute_i, absolute_j],
        "target_symbols": [fragment_a.atoms_list[local_i], fragment_b.atoms_list[local_j]],
        "qc_params": request_qc_params, "orientations": int(orientations),
        "scan_end": scan_end, "scan_step": scan_step, "scan_points": scan_points,
        "default_scan_end_factor": DEFAULT_SCAN_END_FACTOR,
        "default_scan_end_angstrom": DEFAULT_SCAN_END_FACTOR * target_radii,
    }


def validate_scan_restart(output_dir, request):
    """Read-only output and request validation, also usable before sampling.

    Preflight supplies static fields; execution additionally supplies generated
    orientations and identity, retaining the full historical signature check.
    """
    root = Path(output_dir).resolve()
    if not root.exists():
        return
    request_path = root / "request.json"
    if not root.is_dir() or not request_path.is_file():
        raise FileExistsError(f"scan-bond output directory has no resumable request state: {root}")
    try:
        previous = json.loads(request_path.read_text())
        if not isinstance(previous, dict):
            raise ValueError("request state must be an object")
        signature = previous.get("request_signature")
        payload = {key: value for key, value in previous.items() if key != "request_signature"}
        digest = hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":"), default=str).encode("utf-8")).hexdigest()
        if not signature or signature != digest:
            raise ValueError("request signature is missing or does not match its contents")
    except (OSError, ValueError) as exc:
        raise FileExistsError(f"scan-bond request state is unreadable or modified: {request_path}") from exc
    comparison_request = dict(request)
    old_qc = previous.get("qc_params", {})
    new_qc = comparison_request.get("qc_params", {})
    if (isinstance(old_qc, dict) and isinstance(new_qc, dict)
            and "total_multiplicity" not in old_qc
            and new_qc.get("total_multiplicity_source") == "historical_heuristic"):
        # v1 requests predate explicit total-spin provenance. Preserve their
        # old heuristic continuation contract, but do not ignore an explicit
        # total-spin choice made by the caller.
        compatible_qc = dict(new_qc)
        compatible_qc.pop("total_multiplicity", None)
        compatible_qc.pop("total_multiplicity_source", None)
        comparison_request["qc_params"] = compatible_qc
        old_orientations = previous.get("orientation_definitions", [])
        new_orientations = comparison_request.get("orientation_definitions", [])
        if len(old_orientations) == len(new_orientations):
            compatible_orientations = []
            for old, new in zip(old_orientations, new_orientations):
                item = dict(new)
                item.pop("total_multiplicity", None)
                item.pop("total_multiplicity_source", None)
                compatible_orientations.append(item)
            comparison_request["orientation_definitions"] = compatible_orientations
        if "request_signature" in comparison_request:
            # Execution supplies a full signed request, unlike preflight.
            # Compare the signature of the same legacy-compatible payload.
            payload = {key: value for key, value in comparison_request.items()
                       if key != "request_signature"}
            comparison_request["request_signature"] = hashlib.sha256(
                json.dumps(payload, sort_keys=True, separators=(",", ":"), default=str).encode("utf-8")
            ).hexdigest()
    if any(previous.get(key) != value for key, value in comparison_request.items()):
        raise FileExistsError(f"scan-bond output directory belongs to a different request: {root}")


def _write_scan_profile(directory, frames, profile, energy_source="scan.relaxscanact.dat"):
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
        "energy_source": energy_source,
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
        "energy_source": energy_source,
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
        "scan_energy_source": energy_source,
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


def _resolve_total_multiplicity(fragment_a, fragment_b, requested, *, warn=True):
    """Validate a coupled total spin, preserving the legacy heuristic by default."""
    heuristic = fragment_a.merged_with(fragment_b).multiplicity
    if requested is None:
        if warn and (fragment_a.multiplicity != 1 or fragment_b.multiplicity != 1):
            warnings.warn(
                "Open-shell fragment coupling is ambiguous; using the historical "
                f"combined multiplicity heuristic ({heuristic}). Pass --total-multiplicity "
                "to select the total spin explicitly.", RuntimeWarning, stacklevel=2,
            )
        return heuristic, "historical_heuristic"
    if isinstance(requested, bool) or int(requested) != requested or requested < 1:
        raise ValueError("total_multiplicity must be a positive integer")
    low = abs(int(fragment_a.multiplicity) - int(fragment_b.multiplicity)) + 1
    high = int(fragment_a.multiplicity) + int(fragment_b.multiplicity) - 1
    allowed = list(range(low, high + 1, 2))
    if int(requested) not in allowed:
        raise ValueError(
            f"total multiplicity {requested} is incompatible with fragment "
            f"multiplicities {fragment_a.multiplicity} and {fragment_b.multiplicity}; "
            f"allowed total multiplicities are {allowed}"
        )
    return int(requested), "user_explicit"


def _run_scan_bond(input_a, input_b, atoms, orientations, qc_params, output_dir,
                  scan_end=None, scan_step=None, scan_points=None, through=None, reaction_options=None,
                  *, characterize_unique_products=False, _layout=None):
    """Run a scan, optionally continuing through a named reaction-path stage.

    ``through=None`` preserves the original Python scan-plus-relaxation API.
    Explicit ``through="scan"`` performs only the constrained scan for two
    fragments. Single-molecule scans always release the constraint and identify
    the relaxed product; later stages require an identified changed product.
    """
    from pyar.backend_capabilities import normalize_backend_name, supported_geometry_backends
    from pyar.workflows.scan_path import validate_continuation, continue_scan_paths

    software = normalize_backend_name(qc_params.get("software", "orca"))
    if software not in supported_geometry_backends():
        raise ValueError(f"scan-bond requires an energy/gradient backend; unsupported: {software}")
    qc_params = dict(qc_params, software=software)
    if software == "xtb":
        from pyar.backends.xtb_utils import canonical_xtb_model
        qc_params["xtb_model"] = canonical_xtb_model(qc_params.get("xtb_model"))
    if characterize_unique_products and through in (None, "scan"):
        raise ValueError("unique-product characterization requires through='neb' or a later stage")
    if characterize_unique_products and (reaction_options or {}).get("interpolation") == "geodesic":
        raise ValueError("unique-product path-seed characterization does not support geodesic interpolation")
    if through is not None:
        validate_continuation(through, reaction_options or {})
        if software == "orca" and through != "scan":
            from pyar.backends.orca_methods import orca_method
            if orca_method(qc_params.get("method", "BP86"))[0] == "g-xTB":
                raise ValueError("ORCA external g-xTB supports native scanning only; "
                                 "use software=xtb, xtb_model=gxtb for continuation")
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
    intramolecular = input_b is None
    if intramolecular:
        if int(orientations) != 1:
            raise ValueError('Intramolecular scans require one orientation; no fragment placement is performed')
        fragment_a = Molecule.from_xyz(input_a)
        fragment_b = None
        fragment_a.charge = int(qc_params.get('charge_a', qc_params.get('charge', fragment_a.charge)))
        fragment_a.multiplicity = int(qc_params.get('multiplicity_a', qc_params.get('multiplicity', fragment_a.multiplicity)))
        fragment_a.scftype = qc_params.get('scftype_a', qc_params.get('scftype', fragment_a.scftype))
        fragment_a.fragments = [list(range(fragment_a.number_of_atoms))]
        intra_plan = intramolecular_scan_plan(fragment_a, atoms, scan_end, scan_step, scan_points)
        absolute_i, absolute_j = map(int, atoms)
        local_i, local_j = absolute_i, absolute_j
        total_multiplicity, spin_source = fragment_a.multiplicity, 'input_molecule'
        if qc_params.get('total_multiplicity') not in (None, total_multiplicity):
            raise ValueError('Intramolecular total multiplicity must equal the input molecule multiplicity')
        qc_params.update(total_multiplicity=total_multiplicity, total_multiplicity_source=spin_source)
        target_radii = fragment_a.covalent_radius[absolute_i] + fragment_a.covalent_radius[absolute_j]
        default_scan_end = intra_plan['end_distance_angstrom']
        # A product cannot be identified against an unknown input. Fail before
        # allocating the run or starting an expensive constrained scan.
        reactant_identity = molecule_identity_from_xyz(input_a)
        if not reactant_identity or not reactant_identity.get('inchi'):
            raise ValueError('Intramolecular input identity is unavailable; no scan was started')
        reactant_identity_error = None
        orientation = fragment_a.copy()
        orientation.name = orientation.title = 'orientation_000'
        orientation_items = [(0, orientation, None)]
    else:
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
        requested_total = qc_params.get("total_multiplicity")
        recorded_spin_source = qc_params.get("total_multiplicity_source")
        heuristic_total = merged.multiplicity
        if recorded_spin_source == "historical_heuristic":
            if requested_total is not None and int(requested_total) != heuristic_total:
                raise ValueError("recorded historical total multiplicity does not match the fragment heuristic")
            total_multiplicity, spin_source = heuristic_total, "historical_heuristic"
        else:
            total_multiplicity, spin_source = _resolve_total_multiplicity(
                fragment_a, fragment_b, requested_total
            )
        qc_params["total_multiplicity"] = total_multiplicity
        qc_params["total_multiplicity_source"] = spin_source
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
        for _, orientation, _ in orientation_items:
            if orientation is not None:
                orientation.multiplicity = total_multiplicity
    orientation_definitions = [
        {"orientation": index,
         "molecule": None if orientation is None else _fragment_signature(orientation),
         "total_multiplicity": total_multiplicity,
         "total_multiplicity_source": spin_source,
         "generation_error": generation_error}
        for index, orientation, generation_error in orientation_items
    ]
    root = Path(output_dir).resolve()
    request = scan_request_fields(input_a, input_b, fragment_a, fragment_b, atoms,
                                  orientations, qc_params, scan_end, scan_step, scan_points)
    request.update(orientation_definitions=orientation_definitions,
                   reactant_identity=reactant_identity, reactant_identity_error=reactant_identity_error)
    signature_payload = json.dumps(request, sort_keys=True, separators=(",", ":"), default=str)
    request["request_signature"] = hashlib.sha256(signature_payload.encode("utf-8")).hexdigest()
    if _layout is not None:
        request.pop('request_signature')
        request.update(request_schema_version=1, workflow='scan-bond', run_id=_layout.run_id)
        signature_payload = json.dumps(request, sort_keys=True, separators=(",", ":"), default=str)
        request['request_signature'] = hashlib.sha256(signature_payload.encode('utf-8')).hexdigest()
    if _layout is None or (root.exists() and any(root.iterdir())):
        validate_scan_restart(root, request)
    root.mkdir(parents=True, exist_ok=True)
    request_path = root / "request.json"
    if request_path.exists():
        # Validation above established compatibility. Keep the original
        # request/signature as the authority for existing native/v2 commits.
        request = json.loads(request_path.read_text())
    elif _layout is not None:
        from pyar.persistence import write_json
        write_json(request_path, request, overwrite=False)
    else:
        request_temp = root / "request.json.tmp"
        request_temp.write_text(json.dumps(request, indent=2, sort_keys=True, default=str))
        request_temp.replace(request_path)
    if _layout is not None:
        _layout.initialize((input_a,) if intramolecular else (input_a, input_b))
    results = []
    for index, orientation, generation_error in orientation_items:
        directory = (root / f"orientation_{index:03d}" if _layout is None else
                     _layout.layout.path(f"orientation_{index:03d}"))
        execution = (directory if _layout is None else
                     _layout.layout.path(f"orientation_{index:03d}/attempts/attempt_000/execution"))
        scan_dir = execution / "scan" if _layout is None else _layout.layout.path(
            f"orientation_{index:03d}/attempts/attempt_000/execution/scan")
        relax_dir = execution / "relax" if _layout is None else _layout.layout.path(
            f"orientation_{index:03d}/attempts/attempt_000/execution/relax")
        directory.mkdir(exist_ok=True)
        execution.mkdir(parents=True, exist_ok=True)
        result_path = directory / "result.json"
        previous_result = None
        if result_path.is_file():
            try:
                previous_result = json.loads(result_path.read_text())
            except (OSError, json.JSONDecodeError):
                previous_result = None
        committed = _layout.restore(index, request) if _layout is not None else None
        if committed is not None:
            # Immutable records, not the derived per-orientation report, are
            # the recovery authority for a committed scan.
            previous_result = _layout.relaxed(index) or committed
        completed_artifacts = (
            directory / "start.xyz",
            execution / "result_relaxed.xyz",
            scan_dir / "final_scan.xyz",
            scan_dir / "scan_trajectory.xyz",
            scan_dir / "scan_profile.csv",
            scan_dir / "scan_profile.json",
            execution / "ts_candidates" / "metadata.json",
            execution / "ts_candidates" / "highest_scan_energy.xyz",
        )
        if (previous_result
                and previous_result.get("status") == "scan_success_relax_success"
                and all(path.is_file() for path in completed_artifacts)
                and previous_result.get("completed_artifact_hashes") == {
                    str(path.relative_to(directory)): hashlib.sha256(path.read_bytes()).hexdigest()
                    for path in completed_artifacts
                }):
            _save_orientation_result(
                directory, previous_result, qc_params=qc_params, intramolecular=intramolecular,
                completed_artifacts=completed_artifacts, layout=_layout,
            )
            results.append(previous_result)
            continue
        if orientation is None:
            result = {"orientation": index, "status": "scan_failed",
                      "scan_status": "orientation_generation_failed",
                      "error": generation_error}
            _save_orientation_result(
                directory, result, qc_params=qc_params, intramolecular=intramolecular,
                completed_artifacts=completed_artifacts, layout=_layout,
            )
            results.append(result)
            continue
        if committed is None:
            write_xyz(orientation.atoms_list, orientation.coordinates, directory / "start.xyz", job_name=orientation.name, precision=12)
        start_distance = _distance(orientation.coordinates, absolute_i, absolute_j)
        end_distance = float(scan_end if scan_end is not None else default_scan_end)
        try:
            n_points = scan_point_count(start_distance, end_distance, scan_step, scan_points)
        except ValueError as exc:
            result = {"status": "scan_failed", "error": str(exc), "orientation": index}
            _save_orientation_result(
                directory, result, qc_params=qc_params, intramolecular=intramolecular,
                completed_artifacts=completed_artifacts, layout=_layout,
            )
            results.append(result)
            continue
        request_model = OrcaBondScanRequest(absolute_i, absolute_j, start_distance, end_distance, n_points)
        scan, scan_error = _load_or_execute_scan(
            orientation, request_model, scan_dir, qc_params, software=software,
            previous_result=previous_result, committed=committed,
        )
        if scan_error is not None:
            result = {"orientation": index, "status": "scan_failed", "scan_status": "exception",
                      "error": str(scan_error), "target_distance_start_angstrom": start_distance,
                      "target_distance_scan_end_angstrom": end_distance, "scan_points": n_points}
            _save_orientation_result(
                directory, result, qc_params=qc_params, intramolecular=intramolecular,
                completed_artifacts=completed_artifacts, layout=_layout,
            )
            results.append(result)
            continue
        result = {"orientation": index, "status": "scan_failed", "scan_status": scan.status,
                  "target_distance_start_angstrom": start_distance, "target_distance_scan_end_angstrom": end_distance,
                  "scan_points": n_points, "scan_input": str(scan.input_path), "scan_output": str(scan.output_path)}
        from pyar.recovery.stages import scan_diagnostic
        result['optimization_diagnostics'] = scan_diagnostic(orientation, request_model, scan)
        if intramolecular:
            result.update(intra_plan)
            result['orientation'] = index
        if not scan.success:
            _save_orientation_result(
                directory, result, qc_params=qc_params, intramolecular=intramolecular,
                completed_artifacts=completed_artifacts, layout=_layout,
            )
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
            if committed is not None:
                result.update(committed)
            else:
                result.update(_write_scan_profile(
                    scan_dir, scan.frames, scan.profile,
                    "scan.relaxscanact.dat" if software == "orca" else "scan_summary.json"))
            if _layout is not None:
                _layout.commit(index, request, result, directory, execution)
        except Exception as exc:
            result["status"] = "scan_success_analysis_failed"
            result["error"] = str(exc)
            _save_orientation_result(
                directory, result, qc_params=qc_params, intramolecular=intramolecular,
                completed_artifacts=completed_artifacts, layout=_layout,
            )
            results.append(result)
            continue
        if through is not None and not characterize_unique_products and not intramolecular:
            result["status"] = "scan_success"
            result["requested_through"] = through
            _save_orientation_result(
                directory, result, qc_params=qc_params, intramolecular=intramolecular,
                completed_artifacts=completed_artifacts, layout=_layout,
            )
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
            relax_status = _optimize_with_diagnostics(
                relaxed,
                dict(qc_params, geometry_optimizer="native" if software == "orca" else "geometric", gamma=0.0),
                "bond_scan_unbiased_relaxation",
            )
        except BackendExecutionError:
            raise
        except Exception as exc:
            result["status"] = "scan_success_relax_failed"
            result["error"] = str(exc)
            _save_orientation_result(
                directory, result, qc_params=qc_params, intramolecular=intramolecular,
                completed_artifacts=completed_artifacts, layout=_layout,
            )
            results.append(result)
            continue
        finally:
            os.chdir(old_cwd)
        if is_success(relax_status):
            optimized_output = relax_run_dir / "job_relaxed" / "result_relaxed.xyz"
            if not optimized_output.is_file():
                result["status"] = "scan_success_relax_failed"
                result["error"] = f"optimizer reported success but did not write {optimized_output}"
                _save_orientation_result(
                    directory, result, qc_params=qc_params, intramolecular=intramolecular,
                    completed_artifacts=completed_artifacts, layout=_layout,
                )
                results.append(result)
                continue
            try:
                shutil.copy(optimized_output, execution / "result_relaxed.xyz")
            except OSError as exc:
                result["status"] = "scan_success_relax_failed"
                result["error"] = f"could not preserve relaxed geometry: {exc}"
                _save_orientation_result(
                    directory, result, qc_params=qc_params, intramolecular=intramolecular,
                    completed_artifacts=completed_artifacts, layout=_layout,
                )
                results.append(result)
                continue
            result["status"] = "scan_success_relax_success"
            result["relaxed_path"] = str(execution / "result_relaxed.xyz")
            result["target_distance_relaxed_angstrom"] = _distance(relaxed.coordinates, absolute_i, absolute_j)
            result["target_contact_threshold_angstrom"] = TARGET_CONTACT_FACTOR * target_radii
            result["target_bond_present_after_relaxation"] = (
                result["target_distance_relaxed_angstrom"]
                < result["target_contact_threshold_angstrom"]
            )
            if intramolecular:
                target_present = result['target_bond_present_after_relaxation']
                result['target_change_achieved'] = bool(
                    not target_present and result['target_distance_relaxed_angstrom'] >=
                    intra_plan['safe_separation_threshold_angstrom']
                    if intra_plan['scan_direction'] == 'break' else target_present)
            result["reactant_identity"] = reactant_identity
            result["product_identity_changed"] = None
            if reactant_identity is None:
                result["identity_status"] = "reactant_identity_failed"
                result["identity_error"] = reactant_identity_error
            else:
                try:
                    relaxed_identity = molecule_identity_from_xyz(execution / "result_relaxed.xyz")
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
        _save_orientation_result(
            directory, result, qc_params=qc_params, intramolecular=intramolecular,
            completed_artifacts=completed_artifacts, layout=_layout,
        )
        results.append(result)
    if through not in (None, "scan") and not characterize_unique_products:
        continue_scan_paths(root.resolve(), results, qc_params, through, reaction_options or {})
    if _layout is not None:
        registry = _layout.product_registry(results, request)
        if characterize_unique_products:
            if through in (None, "scan"):
                raise ValueError("--characterize-unique-products requires --through neb or a later stage")
            from pyar.workflows.scan_path import characterize_scan_products
            characterization = characterize_scan_products(
                root.resolve(), request, qc_params, through,
                reaction_options or {}, _layout,
            )
            for result in results:
                if result.get("product_role") == "representative":
                    route = characterization.get(result['product_id'], {})
                    result["continuation_status"] = route.get("status", "not_started")
                    if intramolecular:
                        result['reactant_product_connection_confirmed'] = bool(
                            route.get('path_result', {}).get('reactant_product_connection_confirmed', False))
                    (root / f"orientation_{int(result['orientation']):03d}" / "result.json").write_text(
                        json.dumps(result, indent=2, sort_keys=True, default=str))
                elif intramolecular:
                    result['continuation_status'] = ('not_started_no_new_product'
                        if result.get('product_identity_changed') is False else 'blocked_by_product_identity')
        _layout.manifest(results)
    (root / "summary.json").write_text(json.dumps(results, indent=2, sort_keys=True, default=str))
    with (root / "summary.csv").open("w", newline="") as handle:
        fields = sorted({key for result in results for key in result})
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader(); writer.writerows(results)
    if characterize_unique_products:
        complete = all(row.get("status") in {"scan_success_relax_success", "scan_success"}
                       for row in results)
        index = json.loads((root / "products/index.json").read_text())
        if intramolecular:
            complete = complete and all(row.get('identity_status') == 'success' for row in results)
        complete = complete and all(item.get("characterization_status") in {
                                        "reaction_connected", "path_complete_connection_unconfirmed"}
                                    for item in index.get("products", []))
    elif through is None:
        complete = all(row.get("status") == "scan_success_relax_success" for row in results)
    elif through == "scan":
        complete = all(row.get("status") in {"scan_success", "scan_success_relax_success"} for row in results)
    else:
        complete = all(row.get("continuation_status") == "complete" for row in results)
    if intramolecular:
        complete = complete and all(row.get('identity_status') == 'success' for row in results)
    return {"workflow": "scan-bond", "output_dir": str(root), "results": results,
            "status": "complete" if complete else "failed", "requested_through": through,
            "scan_status": ("scan_finished" if all(row.get("scan_status") == "success" for row in results)
                            else "scan_incomplete")}


def run_scan_bond(input_a, input_b, atoms, orientations, qc_params, output_dir,
                  scan_end=None, scan_step=None, scan_points=None, through=None,
                  reaction_options=None, *, characterize_unique_products=False, recover_stale=False):
    """Route existing v1 scans unchanged and guard fresh v2 scan publication."""
    from pyar.persistence import detect_version, RunLock
    from pyar.workflows.scan_layout import ScanLayout
    if input_b is None and through not in (None, 'scan'):
        characterize_unique_products = True
    root = Path(output_dir).resolve()
    def legacy(path):
        validate_scan_restart(path, {})
        if 'request_schema_version' in json.loads((path / 'request.json').read_text()):
            raise ValueError('Incomplete v2 scan publication; recover explicitly')
        return True
    arguments = (input_a, input_b, atoms, orientations, qc_params, output_dir,
                 scan_end, scan_step, scan_points, through, reaction_options)
    version = detect_version(root, 'scan-bond', validate_v1=legacy)
    if characterize_unique_products and version == 1:
        raise ValueError('--characterize-unique-products requires a fresh or v2 scan-bond output directory')
    if version == 1:
        return _run_scan_bond(*arguments)
    # The sibling lock must be creatable before the run root; support nested
    # output paths while still taking the lock before inspecting/creating root.
    root.parent.mkdir(parents=True, exist_ok=True)
    with RunLock(root, recover_stale=recover_stale):
        detect_version(root, 'scan-bond', validate_v1=legacy)
        return _run_scan_bond(*arguments, characterize_unique_products=characterize_unique_products,
                              _layout=ScanLayout(root))
