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


def continue_scan_paths(root, results, qc_params, through, options):
    """Run/reuse stages in order and stop at the first failed scientific gate."""
    from pyar.neb import run_neb, _stage_gate_passed

    if (qc_params["software"] == "orca" and orca_method(qc_params.get("method", "BP86"))[0] == "g-xTB"):
        raise ValueError("ORCA external g-xTB is available for native scans only; "
                         "use --software xtb --xtb-model gxtb for reaction-path continuation")
    root = Path(root).resolve()
    request = json.loads((root / "request.json").read_text())
    last = _CONTINUATION.index(_LAST_STAGE[through])
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
            base = dict(software=qc_params["software"], output=path_dir,
                        method=qc_params.get("method"), basis=qc_params.get("basis"),
                        charge=metadata["charge"], multiplicity=metadata["multiplicity"],
                        nprocs=qc_params.get("nprocs", 1), xtb_model=qc_params.get("xtb_model", "gxtb"),
                        backend_options={"scftype": metadata.get("scftype", "rhf"),
                                         **{name: qc_params[name] for name in ("scf_cycles",) if name in qc_params}},
                        reuse=True, **options)
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
            completed = []
            for stage in _CONTINUATION[:last + 1]:
                try:
                    stage_result = run_neb(start=directory / "start.xyz", end=result["final_scan_path"],
                                           ts_guess=waypoint_path if stage == "neb" else None,
                                           stage=stage, **base)
                    if not _stage_gate_passed(stage_result, stage):
                        result.update(continuation_status="scientific_gate_failed", failed_stage=stage)
                        break
                    completed.append(stage)
                    if stage == "endpoint-frequency":
                        result["reactant_product_connection_confirmed"] = stage_result["reactant_product_connection_confirmed"]
                except Exception as exc:
                    result.update(continuation_status="failed", failed_stage=stage, continuation_error=str(exc))
                    break
            else:
                result["continuation_status"] = "complete"
            result["completed_reaction_stages"] = completed
        except Exception as exc:
            result.update(continuation_status="failed", failed_stage="initialization",
                          continuation_error=str(exc), reactant_product_connection_confirmed=False)
        (directory / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True))
