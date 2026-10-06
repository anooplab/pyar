"""Bridge accepted reaction products to the shared unbiased path engine."""

from __future__ import annotations

import hashlib
import json
import os
import tempfile
from pathlib import Path

import numpy as np

from pyar.backends import write_xyz
from pyar.reaction_trace import load_trace_records
from pyar.reaction_analysis import _persistent_transition_index
from pyar.state.reaction import ReactionRunState
from pyar.workflows.scan_path import (
    REACTION_THROUGH_STAGES,
    run_path_stages,
    validate_continuation,
)


WAYPOINT_FILES = {
    "highest-backend-energy": "highest_backend_energy.xyz",
    "first-topology-change": "first_topology_change.xyz",
    "pre-product": "pre_product_geometry.xyz",
    "highest-total-energy": "highest_total_energy.xyz",
}


def validate_path_request(through, options):
    """Validate the continuation selector and all shared NEB stage options."""
    if through not in REACTION_THROUGH_STAGES:
        raise ValueError(f"Unknown --through value {through!r}; choose {', '.join(REACTION_THROUGH_STAGES)}")
    if through != "react":
        validate_continuation(through, options)
    source = options.get("pathway_ts_source", "highest-backend-energy")
    if source not in WAYPOINT_FILES:
        raise ValueError(f"Unknown pathway TS source {source!r}; choose {', '.join(WAYPOINT_FILES)}")


def preflight_path_request(request, through, options):
    """Preflight the physical, unbiased stages before reaction discovery starts."""
    validate_path_request(through, options)
    if through == "react":
        return []
    from pyar.optimization_request import preflight

    settings = {key: request.qc_params.get(key) for key in (
        "software", "method", "basis", "nprocs", "scf_cycles", "opt_cycles", "opt_threshold")}
    settings.update(geometry_optimizer="geometric", gamma=0.0, opt_target="minimum",
                    xtb_model=request.qc_params.get("xtb_model", "gfn2"))
    molecule = request.reactants[0].merged_with(request.reactants[1])
    requirements = preflight(settings, [molecule], check_example=f"pyar react A.xyz B.xyz --backend {settings['software']} --through {through} --check")
    from pyar.workflows.scan_path import _CONTINUATION, _LAST_STAGE
    includes_ts = _CONTINUATION.index(_LAST_STAGE[through]) >= _CONTINUATION.index("ts")
    if includes_ts and options.get("ts_optimizer", "geometric") == "sella":
        try:
            import sella  # noqa: F401
        except ImportError as exc:
            raise ImportError("Sella TS optimization requires the optional Sella dependency. "
                              "Install with `python -m pip install 'pyar-chem[sella]'`.") from exc
        requirements.append("sella")
    if through != "react" and options.get("interpolation", "linear") == "idpp":
        from pyar.neb import _ase_idpp_api
        _ase_idpp_api()
        requirements.append("ASE IDPP")
    elif through != "react" and options.get("interpolation", "linear") == "geodesic":
        from pyar.neb import _ase_idpp_api, _geodesic_api
        _ase_idpp_api()
        _geodesic_api()
        requirements.extend(("ASE IDPP", "geodesic-interpolate"))
    return sorted(set(requirements))


def _sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _write_json_atomic(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=path.parent,
                                     prefix="." + path.name, suffix=".tmp", delete=False) as stream:
        json.dump(data, stream, indent=2, sort_keys=True, allow_nan=False, default=str)
        stream.write("\n")
        temporary = Path(stream.name)
    os.replace(temporary, path)


def _xyz_data(path):
    from pyar.neb import read_xyz
    elements, coordinates, _comment = read_xyz(path)
    coordinates = np.asarray(coordinates, dtype=float)
    if not elements or coordinates.shape != (len(elements), 3) or not np.isfinite(coordinates).all():
        raise ValueError(f"Invalid or non-finite pathway geometry: {path}")
    return tuple(elements), coordinates


def _trace_reactant_frame(trace_dir, elements):
    """Return the latest traced frame retaining initial coordinate topology."""
    trace_dir = Path(trace_dir)
    records = sorted(load_trace_records(trace_dir), key=lambda item: int(item["step_index"]))
    if not records:
        return None, None
    first_bonds = {tuple(pair) for pair in records[0].get("current_bonds", [])}
    transition = _persistent_transition_index(records)
    candidates = [record for position, record in enumerate(records)
                  if position < transition
                  and {tuple(pair) for pair in record.get("current_bonds", [])} == first_bonds]
    if not candidates:
        return None, None
    record = candidates[-1]
    frame_elements = tuple(record.get("symbols", ()))
    coordinates = np.asarray(record.get("coordinates_angstrom"), dtype=float)
    if frame_elements != elements or coordinates.shape != (len(elements), 3) or not np.isfinite(coordinates).all():
        return None, None
    index = int(record["step_index"])
    path = trace_dir / "steps" / f"step_{index:06d}.xyz"
    return (path if path.is_file() else None), index


def _select_waypoint(candidate_dir, requested):
    candidate_dir = Path(candidate_dir)
    preferred = [WAYPOINT_FILES[requested]]
    if requested == "highest-backend-energy":
        preferred.extend((WAYPOINT_FILES["first-topology-change"], WAYPOINT_FILES["pre-product"]))
    for name in preferred:
        candidate = candidate_dir / name
        if candidate.is_file():
            return candidate, name.removesuffix(".xyz")
    raise ValueError("No usable reaction-trace TS guess; expected one of " + ", ".join(preferred))


def _alternative_routes(root, state, product):
    """Retain raw duplicate-product pathway provenance without characterizing it."""
    identity = (product.get("inchi"), product.get("smiles"))
    if not any(identity):
        return []
    alternatives = []
    for job in state.data.get("completed_jobs", []):
        other = job.get("product_identity") or {}
        same = ((identity[0] and other.get("inchi") == identity[0]) or
                (identity[1] and other.get("smiles") == identity[1]))
        if not same or job.get("status") != "duplicate_product":
            continue
        matches = list((Path(root).resolve() / "reaction").rglob(f"job_{job['job_name']}"))
        alternatives.append({"reaction_job": job["job_name"], "gamma": job.get("gamma"),
                             "trace_directory": str(matches[0] / "reaction_trace") if matches else None,
                             "characterized": False})
    return alternatives


def _prepare_route(root, product, product_number, options, qc_params, previous=None):
    root = Path(root).resolve()
    reaction_dir = root / "reaction"
    product_path = reaction_dir / product["path"]
    summary = product.get("trace_summary") or {}
    candidate_dir = Path(summary.get("candidate_ts_directory", ""))
    if not candidate_dir.is_absolute():
        candidate_dir = (reaction_dir / candidate_dir).resolve()
    job_dir = candidate_dir.parent
    if not product_path.is_file():
        raise ValueError(f"Accepted product geometry is missing: {product_path}")
    waypoint, waypoint_source = _select_waypoint(candidate_dir, options.get("pathway_ts_source", "highest-backend-energy"))
    product_elements, product_coordinates = _xyz_data(product_path)
    waypoint_elements, waypoint_coordinates = _xyz_data(waypoint)
    if waypoint_elements != product_elements:
        raise ValueError("Trace TS guess does not have the accepted product's atom count/order")

    trace_dir = job_dir / "reaction_trace"
    source_path, trace_step = _trace_reactant_frame(trace_dir, product_elements)
    reactant_source = "latest_valid_pre_transition_trace_frame"
    if source_path is None:
        job_name = product["job_name"]
        initial = job_dir.parent / f"trial_{job_name}.xyz"
        if not initial.is_file():
            raise ValueError(f"No reactant-side trace frame or oriented starting geometry for {product['job_name']}")
        source_path = initial
        reactant_source = "initial_oriented_reaction_job_geometry"
        trace_step = None
    start_elements, start_coordinates = _xyz_data(source_path)
    if start_elements != product_elements:
        raise ValueError("Reactant-side geometry does not have the accepted product's atom count/order")

    pathway_id = "route_001"
    route_dir = reaction_dir / "pathways" / f"product_{product_number:03d}" / pathway_id
    start_out, end_out, waypoint_out = route_dir / "reactant.xyz", route_dir / "product.xyz", route_dir / "ts_guess.xyz"
    source_hashes = {"reactant_sha256": _sha256(source_path), "product_sha256": _sha256(product_path),
                     "ts_guess_sha256": _sha256(waypoint),
                     "reaction_trace_sha256": _sha256(trace_dir / "trace.jsonl"),
                     "candidate_metadata_sha256": _sha256(candidate_dir / "metadata.json")}
    if previous is not None and previous.get("inputs") != source_hashes:
        raise ValueError(f"Saved pathway inputs for product_{product_number:03d} changed; refusing unsafe stage reuse")
    expected_options = {key: value for key, value in options.items() if key not in {"through"}}
    if previous is not None and previous.get("path_options") != expected_options:
        raise ValueError("Saved pathway stage settings differ from this invocation")
    route_dir.mkdir(parents=True, exist_ok=True)
    if not start_out.exists():
        write_xyz(start_elements, start_coordinates, start_out, job_name="reactant-side pathway endpoint", precision=12)
    if not end_out.exists():
        write_xyz(product_elements, product_coordinates, end_out, job_name="accepted unbiased product endpoint", precision=12)
    if not waypoint_out.exists():
        write_xyz(waypoint_elements, waypoint_coordinates, waypoint_out,
                  job_name=f"candidate path waypoint from {waypoint_source}", precision=12)
    provenance = {
        "schema_version": 1, "product_id": f"product_{product_number:03d}", "pathway_id": pathway_id,
        "reaction_job": product["job_name"], "gamma": product.get("gamma"),
        "product_identity": {"inchi": product.get("inchi"), "smiles": product.get("smiles")},
        "product_path": str(product_path), "trace_directory": str(trace_dir),
        "stage_output_directory": str(route_dir),
        "candidate_ts_source": waypoint_source, "candidate_ts_is_validated_transition_state": False,
        "reactant_geometry_source": reactant_source, "reactant_trace_step": trace_step,
        "backend_settings": {key: value for key, value in qc_params.items()
                             if key in {"software", "xtb_model", "method", "basis", "charge", "multiplicity", "scftype", "nprocs"}},
        "through": options["through"], "path_options": expected_options, "inputs": source_hashes,
        "prepared_geometry_sha256": {"reactant": _sha256(start_out), "product": _sha256(end_out),
                                     "ts_guess": _sha256(waypoint_out)},
        "status": "prepared",
    }
    _write_json_atomic(route_dir / "pathway.json", provenance)
    return {"product": f"product_{product_number:03d}", "pathway": pathway_id,
            "directory": str(route_dir), "provenance": provenance,
            "start": start_out, "end": end_out, "ts_guess": waypoint_out}


def characterize_reaction(root, request, through, options):
    """Characterize one selected trace pathway for each accepted product."""
    if through == "react":
        return {"status": "not_requested", "pathways": []}
    state = ReactionRunState.read_completed(root, request.restart_request)
    if state is None:
        raise ValueError("Path continuation requires a completed reaction/state.json")
    state_path = Path(root).resolve() / "reaction" / "pathways" / "pathways.json"
    if state_path.is_file():
        existing = json.loads(state_path.read_text(encoding="utf-8"))
        if existing.get("reaction_request") != state.data.get("request"):
            raise ValueError("Pathway state does not match the completed reaction request")
        if existing.get("reaction_state_sha256") not in (None, _sha256(Path(root).resolve() / "reaction" / "state.json")):
            raise ValueError("Completed reaction state changed since pathway characterization")
    else:
        existing = {}
    path_options = {key: value for key, value in options.items() if key in {
        "pathway_ts_source", "images", "max_cycles", "max_gradient", "average_gradient", "spring", "climb", "align",
        "interpolation", "idpp_fmax", "idpp_steps", "geodesic_tol", "geodesic_max_iter",
        "product_relaxation_fmax", "product_relaxation_max_steps", "ts_max_cycles", "ts_optimizer",
        "ts_fmax", "sella_fmax", "sella_internal_coordinates", "irc_max_cycles", "endpoint_max_cycles",
        "imaginary_frequency_threshold", "irc_endpoint_rmsd_tolerance"}}
    if existing is not None and existing.get("path_options") not in (None, path_options):
        raise ValueError("Pathway continuation settings differ from this invocation")
    if existing.get("path_options") not in (None, path_options):
        raise ValueError("Pathway continuation settings differ from this invocation")
    previous_routes = existing.get("routes", {})
    records = []
    failures = []
    summary = {"schema_version": 1, "reaction_request": state.data.get("request"),
               "reaction_state_sha256": _sha256(Path(root).resolve() / "reaction" / "state.json"),
               "through": through, "path_options": path_options, "selected_routes": records,
               "preparation_failures": failures,
               "alternatives": {f"product_{i:03d}": _alternative_routes(root, state, product)
                                for i, product in enumerate(state.data.get("products", []), start=1)},
               "routes": previous_routes,
               "status": "running"}
    _write_json_atomic(state_path, summary)
    for number, product in enumerate(state.data.get("products", []), start=1):
        try:
            product_key = f"product_{number:03d}"
            route = _prepare_route(root, product, number, dict(options, through=through), request.qc_params,
                                   previous=previous_routes.get(product_key))
        except (OSError, ValueError, KeyError) as exc:
            failures.append({"product": product.get("job_name"), "status": "preparation_failed", "error": str(exc)})
            _write_json_atomic(state_path, summary)
            continue
        summary["routes"][route["product"]] = route["provenance"]
        _write_json_atomic(state_path, summary)
        run_options = {key: value for key, value in path_options.items() if key != "pathway_ts_source"}
        result = run_path_stages(start=route["start"], end=route["end"], ts_guess=route["ts_guess"],
                                 output=route["directory"], qc_params=request.qc_params,
                                 through=through, options=run_options)
        route_result = {"product": route["product"], "pathway": route["pathway"],
                        "directory": route["directory"], "through": through,
                        "status": result["status"], "failed_stage": result.get("failed_stage"),
                        "error": result.get("error"),
                        "reactant_product_connection_confirmed": result.get("reactant_product_connection_confirmed", False)}
        route["provenance"].update(status=result["status"], stage_status=result)
        _write_json_atomic(Path(route["directory"]) / "pathway.json", route["provenance"])
        summary["routes"][route["product"]] = route["provenance"]
        records.append(route_result)
        summary["selected_routes"] = records
        summary["preparation_failures"] = failures
        _write_json_atomic(state_path, summary)
    status = "complete" if not failures and all(item["status"] == "complete" for item in records) else "failed"
    if not state.data.get("products"):
        status = "complete_no_products"
    summary["status"] = status
    _write_json_atomic(state_path, summary)
    return summary


def validate_characterization_restart(root, request, through, options):
    """Validate pathway provenance without creating or modifying any files."""
    state_path = Path(root).resolve() / "reaction" / "pathways" / "pathways.json"
    pathway_root = state_path.parent
    reaction_state = ReactionRunState.read_completed(root, request.restart_request)
    if reaction_state is None:
        raise ValueError("Path continuation requires completed reaction discovery")
    if pathway_root.exists() and not state_path.is_file():
        raise ValueError(f"Existing pathway directory has no restart state: {pathway_root}")
    existing = json.loads(state_path.read_text(encoding="utf-8")) if state_path.is_file() else None
    if existing is not None and existing.get("reaction_request") != reaction_state.data.get("request"):
        raise ValueError("Pathway restart state does not match completed reaction discovery")
    if (existing is not None and existing.get("reaction_state_sha256") !=
            _sha256(Path(root).resolve() / "reaction" / "state.json")):
        raise ValueError("Completed reaction state changed since pathway characterization; refusing stage reuse")
    path_options = {key: value for key, value in options.items() if key in {
        "pathway_ts_source", "images", "max_cycles", "max_gradient", "average_gradient", "spring", "climb", "align",
        "interpolation", "idpp_fmax", "idpp_steps", "geodesic_tol", "geodesic_max_iter",
        "product_relaxation_fmax", "product_relaxation_max_steps", "ts_max_cycles", "ts_optimizer",
        "ts_fmax", "sella_fmax", "sella_internal_coordinates", "irc_max_cycles", "endpoint_max_cycles",
        "imaginary_frequency_threshold", "irc_endpoint_rmsd_tolerance"}}
    for route in (existing or {}).get("routes", {}).values():
        if route.get("path_options") != path_options:
            raise ValueError("Pathway restart settings differ from this invocation")
        directory = pathway_root / route.get("product_id", "") / route.get("pathway_id", "")
        route_file = directory / "pathway.json"
        if not route_file.is_file():
            raise ValueError(f"Pathway restart metadata is missing: {route_file}")
        saved = json.loads(route_file.read_text(encoding="utf-8"))
        if saved.get("path_options") != path_options:
            raise ValueError("Pathway stage settings differ from this invocation")
        for name, filename in (("reactant", "reactant.xyz"), ("product", "product.xyz"),
                               ("ts_guess", "ts_guess.xyz")):
            artifact = directory / filename
            if not artifact.is_file() or _sha256(artifact) != saved.get("prepared_geometry_sha256", {}).get(name):
                raise ValueError(f"Pathway input geometry was changed or removed: {artifact}")
    for index, product in enumerate(reaction_state.data.get("products", []), start=1):
        product_key = f"product_{index:03d}"
        saved = (existing or {}).get("routes", {}).get(product_key, {})
        reaction_dir = Path(root).resolve() / "reaction"
        product_path = reaction_dir / product["path"]
        candidate_dir = Path((product.get("trace_summary") or {}).get("candidate_ts_directory", ""))
        if not candidate_dir.is_absolute():
            candidate_dir = (reaction_dir / candidate_dir).resolve()
        waypoint_source = saved.get("candidate_ts_source") or options.get(
            "pathway_ts_source", "highest-backend-energy")
        waypoint, _ = _select_waypoint(candidate_dir, waypoint_source.replace("_", "-"))
        job_dir = candidate_dir.parent
        trace_dir = job_dir / "reaction_trace"
        trace_start, _step = _trace_reactant_frame(trace_dir, _xyz_data(product_path)[0])
        if trace_start is None:
            trace_start = job_dir.parent / f"trial_{product['job_name']}.xyz"
        expected_hashes = {"reactant_sha256": _sha256(trace_start), "product_sha256": _sha256(product_path),
                           "ts_guess_sha256": _sha256(waypoint),
                           "reaction_trace_sha256": _sha256(trace_dir / "trace.jsonl"),
                           "candidate_metadata_sha256": _sha256(candidate_dir / "metadata.json")}
        if saved and expected_hashes != saved.get("inputs"):
            raise ValueError(f"Reaction evidence for {product_key} changed since its pathway was prepared")
        symbols, _coords = _xyz_data(product_path)
        if _xyz_data(trace_start)[0] != symbols or _xyz_data(waypoint)[0] != symbols:
            raise ValueError(f"Pathway atom ordering changed for {product_key}")
        if saved:
            continue
        # Validate request-specific source names/settings before any pathway
        # directory is created; the actual path-stage engine performs its own
        # artifact/dependency-hash checks when execution begins.
        if not trace_dir.joinpath("trace.jsonl").is_file():
            raise ValueError(f"Reaction trace is missing for {product_key}: {trace_dir}")
        if not (candidate_dir / "metadata.json").is_file():
            raise ValueError(f"Candidate waypoint metadata is missing: {candidate_dir / 'metadata.json'}")
    return existing
