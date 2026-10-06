"""Task-oriented NEB, TS, and IRC interfaces over the existing path engine."""

from __future__ import annotations

from pyar.cli_runtime import prepared_run, started_run, finished_run

import argparse
import json
from pathlib import Path
from types import SimpleNamespace

from pyar.backend_capabilities import get_backend_capabilities, normalize_backend_name
from pyar.neb import run_neb
from pyar.optimization_request import preflight, resolve_settings


def build_parser(command, prog=None):
    descriptions = {
        "neb": "Run a validated geomeTRIC reaction-path search between endpoint structures.",
        "ts": "Optimize and validate a transition-state geometry.",
        "irc": "Continue a validated transition-state run along both IRC directions.",
    }
    parser = argparse.ArgumentParser(prog=prog or f"pyar {command}", description=descriptions[command])
    if command == "neb":
        parser.add_argument("start")
        parser.add_argument("end")
        parser.add_argument("--ts-guess", required=True, help="Initial waypoint geometry for the path")
    elif command == "ts":
        parser.add_argument("geometry", help="Initial TS guess geometry")
    else:
        parser.add_argument("run_directory", help="Existing run directory containing a validated TS frequency stage")
    parser.add_argument("--backend", required=True)
    parser.add_argument("--method")
    parser.add_argument("--basis")
    parser.add_argument("--xtb-model", choices=("gfn2", "gxtb"))
    parser.add_argument("--charge", type=int, default=0)
    parser.add_argument("--multiplicity", type=int, help="Infer singlet/doublet from electron parity when omitted")
    parser.add_argument("--nprocs", type=int, default=1)
    parser.add_argument("--output", default={"neb": "neb_run", "ts": "ts_run", "irc": None}[command])
    parser.add_argument("--check", action="store_true")
    advanced = parser.add_argument_group("Advanced stage controls")
    if command == "neb":
        advanced.add_argument("--images", type=int, default=11)
        advanced.add_argument("--interpolation", choices=("linear", "idpp", "geodesic"), default="linear")
        advanced.add_argument("--max-cycles", type=int, default=100)
        advanced.add_argument("--max-gradient", type=float, default=0.05)
        advanced.add_argument("--average-gradient", type=float, default=0.025)
        advanced.add_argument("--spring", type=float, default=1.0)
        advanced.add_argument("--climb", type=float, default=0.5)
        advanced.add_argument("--align", action="store_true")
    if command == "ts":
        advanced.add_argument("--ts-optimizer", choices=("geometric", "sella"), default="geometric")
        advanced.add_argument("--ts-max-cycles", type=int, default=200)
        advanced.add_argument("--ts-fmax", type=float)
        advanced.add_argument("--sella-fmax", type=float, default=0.05)
        advanced.add_argument("--sella-internal-coordinates", action="store_true")
    if command == "irc":
        advanced.add_argument("--irc-max-cycles", type=int, default=200)
    return parser


def _load_for_state(path, charge, multiplicity):
    from pyar.modern_workflow import load_fragments, resolve_states
    molecules = load_fragments([path])
    return resolve_states(molecules, [charge], None if multiplicity is None else [multiplicity])[0]


def _resolve_settings(args):
    backend = normalize_backend_name(args.backend.lower())
    if args.xtb_model is not None and backend != "xtb":
        raise ValueError("--xtb-model requires --backend xtb")
    cap = get_backend_capabilities(backend)
    if not cap.energy_gradient:
        raise ValueError(f"Backend {backend!r} does not provide the Cartesian gradients required for {args.command}")
    if cap.family == "dft_qc":
        if not args.method:
            raise ValueError(f"Backend {backend!r} requires an explicit --method for this path task")
        if backend == "orca":
            from pyar.backends.orca_methods import orca_method
            canonical, is_xtb = orca_method(args.method)
            builtin = canonical.lower() in {"r2scan-3c", "b97-3c", "pbeh-3c", "hf-3c"}
            if not is_xtb and not builtin and not args.basis:
                raise ValueError("ORCA DFT methods require --basis; ORCA xTB/composite methods define their own basis")
            if is_xtb or builtin:
                if args.basis:
                    raise ValueError("The selected ORCA method defines its own basis; omit --basis")
        elif not args.basis:
            raise ValueError(f"Backend {backend!r} requires explicit --basis for this path task")
    namespace = argparse.Namespace(
        backend=backend, geometry_optimizer="geometric", opt_target="minimum", method=args.method,
        basis=args.basis, nprocs=args.nprocs, opt_cycles=100, opt_threshold="normal",
        scf_cycles=None, scf_threshold=None, custom_keywords=None,
        model=None, xtb_model=args.xtb_model if backend == "xtb" else None, gxtb_wrapper=None,
    )
    return resolve_settings(namespace)


def _validate_inputs(command, args):
    if command == "neb":
        paths = [args.start, args.end, args.ts_guess]
        from pyar.neb import read_xyz
        frames = [read_xyz(path) for path in paths]
        if any(frame[0] != frames[0][0] for frame in frames[1:]):
            raise ValueError("NEB start, end, and TS guess must have identical ordered elements")
        return [_load_for_state(args.start, args.charge, args.multiplicity)]
    if command == "ts":
        return [_load_for_state(args.geometry, args.charge, args.multiplicity)]
    run_dir = Path(args.run_directory).expanduser().resolve()
    frequency_summary = run_dir / "frequency_summary.json"
    geometry = run_dir / "frequency_geometry.xyz"
    if not frequency_summary.is_file() or not geometry.is_file():
        raise ValueError("IRC requires an existing frequency-validated TS run directory")
    summary = json.loads(frequency_summary.read_text(encoding="utf-8"))
    if summary.get("first_order_saddle_confirmed") is not True:
        raise ValueError("IRC requires a TS confirmed as a first-order saddle by frequency analysis")
    return [_load_for_state(geometry, args.charge, args.multiplicity)]


def main(command, argv=None, *, prog=None):
    parser = build_parser(command, prog)
    args = parser.parse_args(argv)
    args.command = command
    try:
        from pyar.neb import validate_stage_input_paths, _load_stage, _stage_gate_passed
        from pyar.workflows.scan_path import validate_continuation, REACTION_OPTION_NAMES
        validate_continuation(command, {key: value for key, value in vars(args).items()
                                        if key in REACTION_OPTION_NAMES})
        molecules = _validate_inputs(command, args)
        settings = _resolve_settings(args)
        requirements = preflight(settings, molecules, check_example=f"pyar {command} INPUT --backend {args.backend} --check")
        if command == "neb" and args.interpolation == "geodesic":
            from pyar.neb import _geodesic_api, _ase_idpp_api
            _geodesic_api()
            _ase_idpp_api()
            requirements.append("geodesic-interpolate")
        if command == "neb" and args.interpolation == "idpp":
            from pyar.neb import _ase_idpp_api
            _ase_idpp_api()
            requirements.append("ASE IDPP")
        if command == "ts" and args.ts_optimizer == "sella":
            try:
                __import__("sella")
            except ImportError as exc:
                raise ImportError("Sella TS optimization requires the optional Sella dependency. "
                                  "Install it with `python -m pip install 'pyar-chem[sella]'`.") from exc
            requirements.append("sella")
        if command == "irc":
            output = Path(args.run_directory).expanduser().resolve()
        else:
            output = Path(args.output).expanduser().resolve()
        if command == "irc":
            if args.output and Path(args.output).expanduser().resolve() != output:
                raise ValueError("IRC continues in RUN_DIRECTORY; omit --output")
            # Same physical signature as the canonical engine, without constructing
            # its calculator or touching state. Artifact/dependency hashes are checked.
            qc = dict(software=settings["software"], method=settings.get("method"),
                      basis=settings.get("basis"), charge=args.charge,
                      multiplicity=molecules[0].multiplicity, nprocs=args.nprocs, gamma=0.0)
            if settings["software"] == "xtb":
                qc["xtb_model"] = settings["xtb_model"]
            frequency = _load_stage(output, "frequency", SimpleNamespace(qc_params=qc))
            if not _stage_gate_passed(frequency, "frequency"):
                raise ValueError("IRC requires a frequency-validated first-order saddle")
        else:
            for stage in (("relax", "neb") if command == "neb" else ("ts", "frequency")):
                validate_stage_input_paths(output, stage, start=getattr(args, "start", None),
                                           end=getattr(args, "end", None),
                                           ts_guess=getattr(args, "ts_guess", None),
                                           ts_geometry=getattr(args, "geometry", None))
    except (ValueError, OSError, ImportError, RuntimeError, KeyError) as exc:
        parser.error(f"Preflight failed for {command}: {exc}\nNo calculations were started.")
    prepared_run(args, request=dict(vars(args)), backend=settings,
                 molecules=molecules, requirements=requirements, outputs=[output], state_lists=False)
    if args.check:
        print(f"Preflight: {command}\nBackend: {settings['software']}\nOutput: {output}")
        print("Requirements: " + ", ".join(requirements))
        print("Ready to run.\n--check specified; no calculations were performed.")
        return
    try:
        started_run()
        if command == "neb":
            common = dict(start=args.start, end=args.end, ts_guess=args.ts_guess, output=output,
                             software=settings["software"], method=settings.get("method"), basis=settings.get("basis"),
                             charge=args.charge, multiplicity=molecules[0].multiplicity, nprocs=args.nprocs,
                             xtb_model=settings.get("xtb_model", "gfn2"), images=args.images,
                             interpolation=args.interpolation, max_cycles=args.max_cycles,
                             max_gradient=args.max_gradient, average_gradient=args.average_gradient,
                             spring=args.spring, climb=args.climb, align=args.align)
            # The engine defaults to stage=all. Keep this public task at NEB.
            relaxed = run_neb(stage="relax", **common)
            if not _stage_gate_passed(relaxed, "relax"):
                result = {"status": "scientific_gate_failed", "failed_stage": "relax", "relax": relaxed}
            else:
                result = run_neb(stage="neb", **common)
        elif command == "ts":
            common = dict(software=settings["software"], method=settings.get("method"), basis=settings.get("basis"),
                          charge=args.charge, multiplicity=molecules[0].multiplicity, nprocs=args.nprocs,
                          xtb_model=settings.get("xtb_model", "gfn2"), ts_optimizer=args.ts_optimizer,
                          ts_max_cycles=args.ts_max_cycles, ts_fmax=args.ts_fmax,
                          sella_fmax=args.sella_fmax,
                          sella_internal_coordinates=args.sella_internal_coordinates)
            optimized = run_neb(ts_geometry=args.geometry, output=output, stage="ts", **common)
            if optimized.get("ts_optimization_converged") is not True:
                result = {"status": "scientific_gate_failed", "failed_stage": "ts", "ts": optimized}
            else:
                frequency = run_neb(output=output, stage="frequency", **common)
                result = {"status": "complete" if frequency.get("first_order_saddle_confirmed") else "scientific_gate_failed",
                          "ts": optimized, "frequency": frequency}
        else:
            result = run_neb(output=output, stage="irc", software=settings["software"],
                             method=settings.get("method"), basis=settings.get("basis"),
                             charge=args.charge, multiplicity=molecules[0].multiplicity, nprocs=args.nprocs,
                             xtb_model=settings.get("xtb_model", "gfn2"), irc_max_cycles=args.irc_max_cycles)
    except (ValueError, OSError, ImportError, RuntimeError) as exc:
        parser.error(f"{command} workflow failed after start: {exc}")
    finished_run(result)
    print(json.dumps(result, indent=2, sort_keys=True, default=str))
    failed = result.get("status") in {"failed", "scientific_gate_failed"}
    if command == "neb":
        failed = failed or not _stage_gate_passed(result, "neb")
    elif command == "irc":
        failed = failed or result.get("irc_converged") is not True
    if failed:
        raise SystemExit(1)


def modern_neb(argv=None, *, prog=None):
    return main("neb", argv, prog=prog)


def modern_ts(argv=None, *, prog=None):
    return main("ts", argv, prog=prog)


def modern_irc(argv=None, *, prog=None):
    return main("irc", argv, prog=prog)
