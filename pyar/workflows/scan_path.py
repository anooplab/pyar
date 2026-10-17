"""Stage orchestration from a relaxed scan to independent reaction validation."""

from __future__ import annotations

import inspect
import json
from pathlib import Path

import numpy as np

from pyar.backends import write_xyz
from pyar.backends.orca_scan import parse_xyz_trajectory
from pyar.backends.orca_methods import orca_method


THROUGH_STAGES = ("scan", "neb", "ts", "frequency", "irc", "endpoints", "endpoint-frequency", "all")
_CONTINUATION = ("relax", "neb", "ts", "frequency", "irc", "endpoint-relax", "endpoint-frequency")
_LAST_STAGE = {"neb": "neb", "ts": "ts", "frequency": "frequency", "irc": "irc",
               "endpoints": "endpoint-relax", "endpoint-frequency": "endpoint-frequency",
               "all": "endpoint-frequency"}
REACTION_THROUGH_STAGES = ("react", "neb", "ts", "frequency", "irc", "endpoints",
                           "endpoint-frequency", "all")
REACTION_OPTION_NAMES = {
    "images", "max_cycles", "max_gradient", "average_gradient", "spring", "climb", "align",
    "interpolation", "idpp_fmax", "idpp_steps", "geodesic_tol", "geodesic_max_iter",
    "product_relaxation_fmax", "product_relaxation_max_steps", "ts_max_cycles", "ts_optimizer",
    "ts_fmax", "sella_fmax", "sella_internal_coordinates", "irc_max_cycles", "endpoint_max_cycles",
    "imaginary_frequency_threshold", "irc_endpoint_rmsd_tolerance",
}


def validate_continuation(through, options):
    """Reject ambiguous stage/option requests before creating any scan jobs."""
    from pyar.neb import canonical_neb_parameters, canonical_ts_parameters
    if through not in THROUGH_STAGES:
        raise ValueError(f"Unknown --through stage: {through}; choose {', '.join(THROUGH_STAGES)}")
    unknown = set(options) - REACTION_OPTION_NAMES
    if unknown:
        raise ValueError(f"Unsupported reaction options: {', '.join(sorted(unknown))}")
    canonical_ts_parameters({name: options[name] for name in (
        "ts_optimizer", "ts_max_cycles", "ts_fmax", "sella_fmax", "sella_internal_coordinates") if name in options})
    # Full numerical option validation is also performed by run_neb before
    # each stage; these checks catch user errors before expensive scanning.
    from pyar.neb import _run_neb_in_directory, _STAGE_OPTIONS
    defaults = {name: parameter.default for name, parameter in
                inspect.signature(_run_neb_in_directory).parameters.items()}
    effective = dict(defaults, **options)
    canonical_neb_parameters({name: effective[name] for name in _STAGE_OPTIONS["neb"]})
    images = effective["images"]
    if isinstance(images, bool) or not isinstance(images, int) or images < 3 or images % 2 != 1:
        raise ValueError("images must be an odd integer of at least 3")
    if effective["interpolation"] not in {"linear", "idpp", "geodesic"}:
        raise ValueError("interpolation must be linear, idpp, or geodesic")
    if effective["interpolation"] in {"idpp", "geodesic"}:
        from pyar.neb import _validate_idpp_options
        _validate_idpp_options(effective["idpp_fmax"], effective["idpp_steps"])
    if effective["interpolation"] == "geodesic":
        from pyar.neb import _validate_geodesic_options
        _validate_geodesic_options(effective["geodesic_tol"], effective["geodesic_max_iter"])
    for name in ("max_cycles", "ts_max_cycles", "irc_max_cycles", "endpoint_max_cycles", "product_relaxation_max_steps"):
        value = effective[name]
        if isinstance(value, bool) or not isinstance(value, int) or value < 1:
            raise ValueError(f"{name} must be a positive integer")
    for name in ("max_gradient", "average_gradient", "spring", "climb", "imaginary_frequency_threshold", "irc_endpoint_rmsd_tolerance", "product_relaxation_fmax"):
        value = effective[name]
        if isinstance(value, bool) or not np.isfinite(value) or value <= 0:
            raise ValueError(f"{name} must be positive and finite")


def run_path_stages(*, start, end, ts_guess=None, path_seed=None, output, qc_params,
                    through, options, layout_dependency=None):
    """Run the canonical cumulative path stages for one mapped endpoint pair.

    Both scan-bond and reaction characterization use this adapter; all
    chemistry, artifacts, restart hashes, and scientific gates remain in
    :func:`pyar.neb.run_neb`.
    """
    from pyar.neb import run_neb, _stage_gate_passed

    if through not in _LAST_STAGE:
        raise ValueError(f"Unknown path continuation stage: {through}")
    if (ts_guess is None) == (path_seed is None):
        raise ValueError("Path continuation requires exactly one of ts_guess or path_seed")
    layout = None
    if layout_dependency is not None:
        from pyar.workflows.path_layout import PathLayout
        layout = PathLayout(output, start, end, ts_guess, qc_params, options,
                            layout_dependency, path_seed=path_seed)
        start, end, route_seed = layout.inputs
        if path_seed is None:
            ts_guess = route_seed
        else:
            path_seed = route_seed
        output = layout.execution
    base = {key: qc_params.get(key) for key in
            ("software", "method", "basis", "charge", "multiplicity", "nprocs", "xtb_model")}
    base["nprocs"] = base["nprocs"] or 1
    base["xtb_model"] = base["xtb_model"] or "gfn2"
    base["backend_options"] = {key: qc_params[key] for key in
                                ("scf_cycles", "scftype", "total_multiplicity_source") if key in qc_params}
    base["output"] = output
    base["reuse"] = True
    stage_options = dict(options)
    association_endpoint_mode = bool(stage_options.pop("association_endpoint_mode", False))
    seed_options = {}
    if path_seed is not None:
        for key in ("path_refinement", "path_seed_maximum_frame_index"):
            if key in stage_options:
                seed_options[key] = stage_options.pop(key)
    base.update(stage_options)
    last = _CONTINUATION.index(_LAST_STAGE[through])
    results = {"status": "running", "completed_stages": [],
               "execution_directory": str(Path(output).resolve())}
    if layout is not None:
        results["path_layout_root"] = str(layout.root)
    for stage in _CONTINUATION[:last + 1]:
        try:
            if layout is not None:
                layout.verify(stage)
            result = run_neb(start=start, end=end,
                             ts_guess=ts_guess if stage == "neb" and path_seed is None else None,
                             path_seed=path_seed if stage == "neb" else None,
                             stage=stage, **(seed_options if stage == "neb" else {}), **base)
            results[stage] = result
            gate_passed = _stage_gate_passed(result, stage)
            if association_endpoint_mode and stage == "endpoint-frequency":
                # Keep the IRC, endpoint optimization/frequency, and distinct
                # endpoint gates. Fragment-aware chemical assignment is done
                # by reaction_characterization after this stage is committed.
                gate_passed = all(result.get(key) is True for key in (
                    "irc_converged", "endpoints_are_minima", "reference_endpoints_distinct",
                    "observed_endpoints_distinct",
                ))
                result["association_endpoint_prerequisites_passed"] = gate_passed
            if not gate_passed:
                results.update(status="scientific_gate_failed", failed_stage=stage)
                return results
            if layout is not None:
                layout.commit(stage, result)
            results["completed_stages"].append(stage)
        except Exception as exc:
            results.update(status="failed", failed_stage=stage, error=str(exc))
            return results
    results["status"] = "complete"
    results["reactant_product_connection_confirmed"] = bool(
        results.get("endpoint-frequency", {}).get("reactant_product_connection_confirmed", False)
    )
    return results


def continue_scan_paths(root, results, qc_params, through, options):
    """Run/reuse stages in order and stop at the first failed scientific gate."""
    if (qc_params["software"] == "orca" and orca_method(qc_params.get("method", "BP86"))[0] == "g-xTB"):
        raise ValueError("ORCA external g-xTB is available for native scans only; "
                         "use --software xtb --xtb-model gxtb for reaction-path continuation")
    root = Path(root).resolve()
    request = json.loads((root / "request.json").read_text())
    for result in results:
        if result.get("scan_status") != "success" or not result.get("scan_profile_json"):
            result["continuation_status"] = "blocked_by_scan"
            continue
        directory = root / f"orientation_{result['orientation']:03d}"
        try:
            path_dir = directory / "reaction_path"
            path_dir.mkdir(exist_ok=True)
            metadata = next(item["molecule"] for item in request["orientation_definitions"]
                            if item["orientation"] == result["orientation"])
            path_qc = dict(qc_params, charge=metadata["charge"],
                           multiplicity=metadata["multiplicity"],
                           scftype=metadata.get("scftype", "rhf"),
                           total_multiplicity_source=next(
                               item.get("total_multiplicity_source", "historical_heuristic")
                               for item in request["orientation_definitions"]
                               if item["orientation"] == result["orientation"]
                           ))
            frames = parse_xyz_trajectory(result["trajectory_path"], len(metadata["atoms"]), metadata["atoms"])
            profile = json.loads(Path(result["scan_profile_json"]).read_text())["points"]
            if len(frames) > 2:
                waypoint_index = max(range(1, len(frames) - 1), key=lambda index: profile[index]["energy_hartree"])
                waypoint = frames[waypoint_index][1]
                waypoint_source = "highest_energy_interior_scan_frame"
            else:
                waypoint_index = None
                waypoint = (frames[0][1] + frames[-1][1]) / 2
                waypoint_source = "interpolated_scan_midpoint"
            waypoint_path = directory / "scan_waypoint.xyz"
            write_xyz(metadata["atoms"], waypoint, waypoint_path, job_name="NEB initializer", precision=12)
            result.update(requested_through=through, scan_waypoint_source=waypoint_source,
                          scan_waypoint_frame_index=waypoint_index, reaction_path_dir=str(path_dir),
                          reactant_product_connection_confirmed=False)
            for stale in ("failed_stage", "continuation_error"):
                result.pop(stale, None)
            extra = {}
            if (root / 'layout.json').exists():
                from pyar.persistence import sha256_file
                extra['layout_dependency'] = sha256_file(directory / 'scan_commit.json')
            path_result = run_path_stages(start=directory / "start.xyz", end=result["final_scan_path"],
                                          ts_guess=waypoint_path, output=path_dir,
                                          qc_params=path_qc, through=through, options=options, **extra)
            result["path_stage_results"] = {key: value for key, value in path_result.items()
                                            if key not in {"status", "completed_stages"}}
            result["continuation_status"] = path_result["status"]
            result["completed_reaction_stages"] = path_result["completed_stages"]
            if path_result["status"] != "complete":
                result["failed_stage"] = path_result.get("failed_stage")
                result["continuation_error"] = path_result.get("error")
            result["reactant_product_connection_confirmed"] = path_result.get(
                "reactant_product_connection_confirmed", False)
        except Exception as exc:
            result.update(continuation_status="failed", failed_stage="initialization",
                          continuation_error=str(exc), reactant_product_connection_confirmed=False)
        (directory / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True))


def characterize_scan_products(root, request, qc_params, through, options, layout):
    """Characterize one scan-derived full trajectory per unique product.

    Product discovery and its immutable scan/relaxation commits are complete
    before this function is called. Each product is isolated so one failed
    route cannot erase or block its siblings.
    """
    import hashlib
    from pyar.persistence import ArtifactRef, read_json, write_json

    root = Path(root).resolve()
    if through == "scan":
        raise ValueError("unique-product characterization requires --through neb or later")
    if qc_params.get("software") == "orca" and orca_method(qc_params.get("method", "BP86"))[0] == "g-xTB":
        raise ValueError("ORCA external g-xTB supports native scanning only; use standalone xTB for path characterization")
    definition_rows = {int(item["orientation"]): item
                       for item in request.get("orientation_definitions", [])}
    summary = {}
    index_path = root / "products/index.json"
    index = read_json(index_path)
    for product in index.get("products", []):
        product_id = product["product_id"]
        representative = int(product["representative_orientation"])
        orientation = root / f"orientation_{representative:03d}"
        execution = orientation / "attempts/attempt_000/execution"
        start = orientation / "start.xyz"
        end = execution / "result_relaxed.xyz"
        seed = execution / "scan/scan_trajectory.xyz"
        profile_path = execution / "scan/scan_profile.json"
        state_path = root / "products" / product_id / "characterization.json"
        state = {"product_id": product_id, "representative_orientation": representative,
                 "status": "product_identified", "product_record_sha256": product["product_record_sha256"],
                 "status_history": ["scan_finished", "product_identified"], "through": through}
        try:
            product_record_path = root / product["product_record_path"]
            read_json(product_record_path)
            if hashlib.sha256(product_record_path.read_bytes()).hexdigest() != product["product_record_sha256"]:
                raise ValueError(f"Product record checksum mismatch for {product_id}")
            if not all(path.is_file() for path in (start, end, seed, profile_path)):
                raise ValueError("representative scan trajectory or relaxed product artifact is missing")
            scan_commit = read_json(layout.commit_path(representative))
            relax_commit = read_json(layout.layout.path(f"orientation_{representative:03d}/relax_commit.json"))
            for record in (scan_commit, relax_commit):
                for ref in record["artifacts"]:
                    ArtifactRef.from_dict(ref).verify(root)
            definition_row = definition_rows.get(representative)
            definition = None if definition_row is None else definition_row.get("molecule")
            if definition is None:
                raise ValueError("representative orientation is absent from the scan request")
            path_qc = dict(qc_params, charge=definition.get("charge", 0),
                           multiplicity=definition.get("multiplicity", 1),
                           scftype=definition.get("scftype", "rhf"),
                           total_multiplicity=definition_row.get("total_multiplicity", qc_params.get("total_multiplicity")),
                           total_multiplicity_source=definition_row.get(
                               "total_multiplicity_source", "historical_heuristic"))
            path_options = dict(options)
            requested_interpolation = path_options.get("interpolation", "linear")
            if requested_interpolation == "geodesic":
                raise ValueError("geodesic interpolation is not supported with scan path-seed characterization")
            path_options["interpolation"] = "linear"
            path_options["path_refinement"] = ("idpp" if requested_interpolation == "idpp" else "none")
            profile = json.loads(profile_path.read_text()).get("points", [])
            maximum_index = max(range(len(profile)), key=lambda i: profile[i]["energy_hartree"])
            path_options["path_seed_maximum_frame_index"] = maximum_index
            state["status"] = "path_attempted"
            state["status_history"].append("path_attempted")
            state["scan_commit_sha256"] = hashlib.sha256(layout.commit_path(representative).read_bytes()).hexdigest()
            state["relax_commit_sha256"] = hashlib.sha256(
                layout.layout.path(f"orientation_{representative:03d}/relax_commit.json").read_bytes()).hexdigest()
            state["physical_method"] = path_qc
            state["path_options"] = path_options
            state["scan_path_seed_sha256"] = hashlib.sha256(seed.read_bytes()).hexdigest()
            path_root = root / "products" / product_id / "characterization"
            metadata_path = path_root / "inputs/metadata.json"
            if metadata_path.is_file():
                expected = {"inputs/reactant.xyz": hashlib.sha256(start.read_bytes()).hexdigest(),
                            "inputs/product.xyz": hashlib.sha256(end.read_bytes()).hexdigest(),
                            "inputs/path_seed.xyz": state["scan_path_seed_sha256"]}
                saved_refs = read_json(metadata_path)["inputs"]
                for ref in saved_refs:
                    ArtifactRef.from_dict(ref).verify(path_root)
                saved = {ref["path"]: ref["sha256"] for ref in saved_refs}
                if saved != expected:
                    # A newly recovered orientation may become representative.
                    # Preserve the earlier representative's immutable inputs.
                    path_root = path_root / "representatives" / f"orientation_{representative:03d}"
            state["path_directory"] = path_root.relative_to(root).as_posix()
            write_json(state_path, state)
            result = run_path_stages(
                start=start, end=end, path_seed=seed,
                output=path_root,
                qc_params=path_qc, through=through, options=path_options,
                layout_dependency=product["product_record_sha256"],
            )
            state["path_result"] = result
            if result.get("status") == "scientific_gate_failed":
                state["status"] = "path_scientific_gate_failed"
            elif result.get("status") == "complete" and result.get("reactant_product_connection_confirmed"):
                state["status"] = "reaction_connected"
            elif result.get("status") == "complete":
                state["status"] = "path_complete_connection_unconfirmed"
            else:
                state["status"] = "path_failed"
            from pyar.neb import read_xyz
            from pyar.workflows.reaction_report import (
                build_reaction_characterization_report, write_reaction_characterization_report,
            )
            report = build_reaction_characterization_report(
                route_root=path_root,
                product={"product_id": product_id, "symbols": read_xyz(end)[0]},
                provenance={"product_id": product_id,
                            "representative_orientation": representative,
                            "discovery_product_id": product_id,
                            "product_record_sha256": product["product_record_sha256"],
                            "physical_method": path_qc, "path_options": path_options,
                            "scan_commit_sha256": state["scan_commit_sha256"],
                            "relax_commit_sha256": state["relax_commit_sha256"]},
                path_result=result,
                route_result={"status": state["status"], "through": through,
                              "reactant_product_connection_confirmed": bool(
                                  result.get("reactant_product_connection_confirmed"))},
                alternatives=[{"orientation": int(value)} for value in
                              product.get("alternative_orientations", []) if int(value) != representative],
                reference_mode="minima",
                discovery={"product_record": product,
                           "representative_orientation": representative,
                           "member_orientations": product.get("orientation_ids", [])},
                scan_evidence={"trajectory": str(seed), "trajectory_sha256": state["scan_path_seed_sha256"],
                               "profile": str(profile_path),
                               "profile_sha256": hashlib.sha256(profile_path.read_bytes()).hexdigest(),
                               "scan_commit": str(layout.commit_path(representative)),
                               "scan_commit_sha256": state["scan_commit_sha256"],
                               "relax_commit": str(layout.layout.path(
                                   f"orientation_{representative:03d}/relax_commit.json")),
                               "relax_commit_sha256": state["relax_commit_sha256"]},
            )
            report_path = root / "products" / product_id / "reaction_characterization.json"
            report["report_pathway_directory"] = str(path_root)
            write_reaction_characterization_report(report_path, report)
            state["reaction_characterization_report"] = report_path.relative_to(root).as_posix()
            state["reaction_characterization_report_sha256"] = hashlib.sha256(
                report_path.read_bytes()).hexdigest()
        except Exception as exc:
            state.update(status="path_failed", error=str(exc))
        if not state["status_history"] or state["status_history"][-1] != state["status"]:
            state["status_history"].append(state["status"])
        write_json(state_path, state)
        product["characterization_status"] = state["status"]
        product["characterization_path"] = state_path.relative_to(root).as_posix()
        product["characterization_sha256"] = hashlib.sha256(state_path.read_bytes()).hexdigest()
        summary[product_id] = state
        write_json(index_path, index)
    return summary
