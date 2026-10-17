"""Auditable per-product reports and validated reaction-path exports."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re

import numpy as np


def _sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_trajectory(path, expected_symbols):
    from pyar.xyz_trajectory import parse_xyz_trajectory

    frames = parse_xyz_trajectory(path, len(expected_symbols), expected_symbols)
    result = []
    for symbols, coordinates, comment in frames:
        match = re.search(r"(?:^|\s)energy_hartree=([-+0-9.eE]+)", comment)
        result.append({"coordinates": np.asarray(coordinates, dtype=float),
                       "energy_hartree": float(match.group(1)) if match else None})
    return result


def _aligned_rmsd(first, second):
    first, second = np.asarray(first, dtype=float), np.asarray(second, dtype=float)
    if first.shape != second.shape:
        return float("inf")
    first, second = first - first.mean(axis=0), second - second.mean(axis=0)
    left, _, right = np.linalg.svd(first.T @ second)
    correction = np.eye(3)
    correction[-1, -1] = np.linalg.det(left @ right)
    return float(np.sqrt(np.mean(np.sum((first @ left @ correction @ right - second) ** 2, axis=1))))


def assemble_validated_path(execution_directory, symbols, endpoint_stage, assignment,
                            destination, *, join_tolerance=1.0e-5):
    """Write reactant→TS→product frames only for an explicit branch assignment."""
    if (not isinstance(assignment, (list, tuple)) or len(assignment) != 2
            or sorted(assignment) != [0, 1]):
        raise ValueError("Validated reaction path requires a bijective reactant/product IRC assignment")
    execution = Path(execution_directory)
    reactant_index, product_index = map(int, assignment)
    branches = ("forward", "backward")
    reactant_direction, product_direction = branches[reactant_index], branches[product_index]
    reactant_frequency = endpoint_stage[f"{reactant_direction}_frequency"]
    product_frequency = endpoint_stage[f"{product_direction}_frequency"]
    reactant_endpoint = np.asarray(reactant_frequency["evaluated_coordinates_angstrom"], dtype=float)
    product_endpoint = np.asarray(product_frequency["evaluated_coordinates_angstrom"], dtype=float)
    from pyar.neb import read_xyz
    ts_symbols, ts_coordinates, _ = read_xyz(execution / "frequency_geometry.xyz")
    if list(ts_symbols) != list(symbols):
        raise ValueError("Validated TS atom mapping differs from endpoint mapping")
    def branch_to_minimum(direction):
        irc_branch = _read_trajectory(execution / f"irc_{direction}_path.xyz", symbols)
        relaxation = _read_trajectory(execution / f"irc_{direction}_relaxed_path.xyz", symbols)
        if (_aligned_rmsd(irc_branch[0]["coordinates"], ts_coordinates) > join_tolerance
                or _aligned_rmsd(irc_branch[-1]["coordinates"], relaxation[0]["coordinates"]) > join_tolerance):
            raise ValueError(f"IRC {direction} path does not join its TS or endpoint-relaxation trajectory")
        verified_endpoint = np.asarray(
            endpoint_stage[f"{direction}_frequency"]["evaluated_coordinates_angstrom"], dtype=float)
        if _aligned_rmsd(relaxation[-1]["coordinates"], verified_endpoint) > 1.0e-4:
            raise ValueError(f"IRC {direction} endpoint relaxation does not join its frequency geometry")
        return irc_branch, relaxation

    reactant_irc, reactant_relaxation = branch_to_minimum(reactant_direction)
    product_irc, product_relaxation = branch_to_minimum(product_direction)
    if _aligned_rmsd(reactant_endpoint, reactant_relaxation[-1]["coordinates"]) > 1.0e-4:
        raise ValueError("Reactant frequency geometry differs from its branch endpoint")
    if _aligned_rmsd(product_endpoint, product_relaxation[-1]["coordinates"]) > 1.0e-4:
        raise ValueError("Product frequency geometry differs from its branch endpoint")
    frames = []
    pieces = (
        (reversed(reactant_relaxation), f"endpoint_relax_{reactant_direction}"),
        (reversed(reactant_irc), f"irc_{reactant_direction}"),
        (product_irc, f"irc_{product_direction}"),
        (product_relaxation, f"endpoint_relax_{product_direction}"),
    )
    for piece, source in pieces:
        for frame in piece:
            if frames and _aligned_rmsd(frames[-1]["coordinates"], frame["coordinates"]) <= join_tolerance:
                continue
            frames.append({**frame, "source": source})
    if len(frames) < 3:
        raise ValueError("Validated path has too few frames after removing duplicate joins")
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="utf-8") as stream:
        for index, frame in enumerate(frames):
            coordinates = frame["coordinates"]
            if coordinates.shape != (len(symbols), 3) or not np.isfinite(coordinates).all():
                raise ValueError("Validated path contains invalid coordinates")
            energy = frame.get("energy_hartree")
            energy_text = "" if energy is None else f" energy_hartree={float(energy):.14f}"
            stream.write(f"{len(symbols)}\nframe={index} source={frame['source']} units=angstrom{energy_text}\n")
            for symbol, xyz in zip(symbols, coordinates):
                stream.write(f"{symbol:<2s} {xyz[0]: .14f} {xyz[1]: .14f} {xyz[2]: .14f}\n")
    return {"path": destination.name, "sha256": _sha256(destination),
            "frame_count": len(frames), "units": "angstrom",
            "energy_unit": "hartree", "energy_missing_frames": sum(
                frame.get("energy_hartree") is None for frame in frames),
            "reactant_branch": reactant_direction, "product_branch": product_direction,
            "endpoint_assignment": list(assignment)}


def _stage_details(path_result, route_root):
    details, artifacts, energies, optimizer_steps, hessian_evaluations = {}, {}, {}, {}, {}
    artifacts_verified = True
    for name in ("relax", "neb", "ts", "frequency", "irc", "endpoint-relax", "endpoint-frequency"):
        stage = path_result.get(name)
        if not isinstance(stage, dict):
            continue
        detail = {key: value for key, value in stage.items() if key not in {"artifacts", "inputs"}}
        details[name] = detail
        detail["artifact_refs"] = []
        detail["input_refs"] = []
        energies[name] = stage.get("backend_energy_gradient_evaluations")
        optimizer_steps[name] = stage.get("optimizer_steps")
        hessian_evaluations[name] = stage.get("hessian_evaluations")
        for path, digest in (stage.get("artifacts") or {}).items():
            artifact = Path(path)
            if artifact.is_file():
                resolved = artifact.resolve()
                try:
                    relative = resolved.relative_to(Path(route_root).resolve()).as_posix()
                except ValueError:
                    relative = str(resolved)
                actual = _sha256(resolved)
                artifacts[relative] = {"sha256": actual, "stage_recorded_sha256": digest,
                                       "verified": actual == digest}
                detail["artifact_refs"].append(relative)
                artifacts_verified = artifacts_verified and actual == digest
            else:
                artifacts_verified = False
                artifacts[str(artifact)] = {"sha256": None, "stage_recorded_sha256": digest,
                                             "verified": False}
                detail["artifact_refs"].append(str(artifact))
        for path, digest in (stage.get("inputs") or {}).items():
            source = Path(path)
            if source.is_file():
                actual = _sha256(source)
                try:
                    relative = source.resolve().relative_to(Path(route_root).resolve()).as_posix()
                except ValueError:
                    relative = str(source.resolve())
                detail["input_refs"].append({"path": relative, "sha256": actual,
                                             "verified": actual == digest})
                artifacts_verified = artifacts_verified and actual == digest
            else:
                detail["input_refs"].append({"path": str(source), "sha256": None, "verified": False})
                artifacts_verified = False
    required = {"relax", "neb", "ts", "frequency", "irc", "endpoint-relax", "endpoint-frequency"}
    stage_records_complete = all(name in details and details[name].get("artifact_refs") for name in required)
    artifacts_verified = artifacts_verified and bool(artifacts) and stage_records_complete
    return details, artifacts, energies, optimizer_steps, hessian_evaluations, artifacts_verified, stage_records_complete


def _artifact_named(artifacts, filename):
    return next((value for path, value in artifacts.items() if Path(path).name == filename), None)


def _verify_scan_evidence(evidence):
    if not evidence:
        return True
    verified = True
    for name in ("trajectory", "profile", "scan_commit", "relax_commit"):
        path, expected = evidence.get(name), evidence.get(name + "_sha256")
        if path is None and expected is None:
            continue
        if not path or not expected or not Path(path).is_file() or _sha256(path) != expected:
            verified = False
    return verified


def build_reaction_characterization_report(*, route_root, product, provenance,
                                           path_result, route_result, alternatives,
                                           reference_mode="minima", reactant_reference=None,
                                           discovery=None, scan_evidence=None):
    """Build a JSON-safe report; validated energetics require compatible gates."""
    route_root = Path(route_root).resolve()
    execution = Path(path_result.get("execution_directory", route_root)).resolve()
    (stage_results, artifacts, evaluation_counts, optimizer_steps,
     hessian_evaluations, artifacts_verified, stage_records_complete) = _stage_details(path_result, route_root)
    endpoint_stage = stage_results.get("endpoint-frequency", {})
    frequency = stage_results.get("frequency", {})
    association = path_result.get("association_endpoint_validation", {})
    scan_evidence_verified = _verify_scan_evidence(scan_evidence)
    connected = bool(route_result.get("reactant_product_connection_confirmed"))
    assignment = association.get("reactant_product_assignment") or endpoint_stage.get("endpoint_assignment")
    physical = [stage.get("qc_params") for stage in stage_results.values() if stage.get("qc_params")]
    physical_compatible = bool(physical) and all(item == physical[0] for item in physical[1:])
    if reference_mode == "separated-reactants" and reactant_reference:
        reference_method = reactant_reference.get("request", {}).get("physical_method", {})
        if physical:
            stage_method = physical[0]
            keys = ("software", "method", "basis", "xtb_model")
            ref_settings = {key: reference_method.get(key) for key in keys}
            stage_settings = {key: stage_method.get(key) for key in keys}
            if str(ref_settings.get("software", "")).lower() == "xtb":
                from pyar.backends.xtb_utils import canonical_xtb_model
                ref_settings["xtb_model"] = canonical_xtb_model(ref_settings.get("xtb_model"))
                stage_settings["xtb_model"] = canonical_xtb_model(stage_settings.get("xtb_model"))
                ref_settings["method"] = stage_settings["method"] = None
                ref_settings["basis"] = stage_settings["basis"] = None
            physical_compatible = physical_compatible and ref_settings == stage_settings
    saddle = frequency.get("first_order_saddle_confirmed") is True
    irc_ok = stage_results.get("irc", {}).get("irc_converged") is True
    endpoints_ok = endpoint_stage.get("endpoints_are_minima") is True
    gates_passed = (connected and saddle and irc_ok and endpoints_ok and physical_compatible
                    and artifacts_verified and stage_records_complete and scan_evidence_verified)
    energies = {"unit": "hartree", "reference_mode": reference_mode,
                "E_R": None, "E_TS": frequency.get("energy_hartree"), "E_P": None,
                "delta_E_rxn": None, "delta_E_forward_activation": None,
                "delta_E_reverse_activation": None, "validated": False}
    if gates_passed and assignment and sorted(assignment) == [0, 1]:
        reactant_index, product_index = map(int, assignment)
        directions = ("forward", "backward")
        reactant_frequency = endpoint_stage.get(f"{directions[reactant_index]}_frequency", {})
        product_frequency = endpoint_stage.get(f"{directions[product_index]}_frequency", {})
        reactant_energy = (reactant_reference.get("reference_energy_hartree")
                           if reference_mode == "separated-reactants" and reactant_reference else
                           reactant_frequency.get("energy_hartree"))
        product_energy = product_frequency.get("energy_hartree")
        ts_energy = frequency.get("energy_hartree")
        if all(value is not None and np.isfinite(float(value))
               for value in (reactant_energy, product_energy, ts_energy)):
            energies.update(E_R=float(reactant_energy), E_TS=float(ts_energy), E_P=float(product_energy),
                            delta_E_rxn=float(product_energy - reactant_energy),
                            delta_E_forward_activation=float(ts_energy - reactant_energy),
                            delta_E_reverse_activation=float(ts_energy - product_energy), validated=True,
                            physical_settings_compatible=True)
    report = {
        "schema_version": 1, "product_id": provenance.get("product_id"),
        "pathway_id": provenance.get("pathway_id"),
        "status": route_result.get("association_status") or route_result.get("status"),
        "through": route_result.get("through"), "connection_confirmed": connected,
        "discovery": discovery or {}, "representative": provenance,
        "alternative_routes": alternatives or [],
        "trajectory_evidence": {"discovery": scan_evidence or {
            "trace_directory": provenance.get("trace_directory"),
            "trace_file": (str(Path(provenance["trace_directory"]) / "trace.jsonl")
                           if provenance.get("trace_directory") else None),
            "trace_sha256": provenance.get("inputs", {}).get("reaction_trace_sha256")},
            "neb_approximation": _artifact_named(artifacts, "neb_path.xyz"),
            "irc": _artifact_named(artifacts, "irc_path.xyz")},
        "stage_parameters_and_results": stage_results,
        "ts_imaginary_frequencies_cm-1": frequency.get("imaginary_frequencies_cm-1"),
        "irc_branches_and_endpoint_assignment": association or {
            "endpoint_assignment": endpoint_stage.get("endpoint_assignment"),
            "forward_frequency": endpoint_stage.get("forward_frequency"),
            "backward_frequency": endpoint_stage.get("backward_frequency")},
        "evaluation_costs": {"backend_energy_gradient_evaluations_by_stage": evaluation_counts,
                             "total_backend_energy_gradient_evaluations": sum(
                                 int(value) for value in evaluation_counts.values()
                                 if isinstance(value, (int, float))),
                             "wall_seconds_by_stage": {key: value.get("wall_seconds")
                                                        for key, value in stage_results.items()},
                             "optimizer_steps_by_stage": optimizer_steps,
                             "hessian_evaluations_by_stage": hessian_evaluations},
        "physical_settings_compatible": physical_compatible,
        "stage_records_complete": stage_records_complete,
        "scan_evidence_verified": scan_evidence_verified,
        "artifacts_verified": artifacts_verified,
        "artifacts": artifacts, "electronic_energetics": energies,
        "validated_reaction_path": None,
    }
    if gates_passed and assignment:
        try:
            report["validated_reaction_path"] = assemble_validated_path(
                execution, product.get("symbols") or report.get("symbols") or
                stage_results["relax"].get("reactant_symbols") or
                _symbols_from_file(execution / "frequency_geometry.xyz"), endpoint_stage, assignment,
                route_root / "validated_reaction_path.xyz",
            )
        except (OSError, ValueError, KeyError, IndexError) as exc:
            report["validated_path_error"] = f"{type(exc).__name__}: {exc}"
            report["validated_reaction_path"] = None
    report["report_artifacts_verified"] = artifacts_verified and scan_evidence_verified
    return report


def _symbols_from_file(path):
    from pyar.neb import read_xyz
    return read_xyz(path)[0]


def write_reaction_characterization_report(path, report):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False, default=str) + "\n",
                         encoding="utf-8")
    temporary.replace(path)
    return path
