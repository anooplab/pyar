"""Modern CLI for solute-centred first-shell microsolvation."""

import argparse

from pyar.cli_runtime import prepared_run, started_run, finished_run
from pyar.backend_errors import BackendExecutionError
from pyar.microsolvation import MicrosolvationRequest
from pyar.microsolvation.surface import build_solute_surface
from pyar.modern_workflow import add_backend_arguments, add_fragment_arguments, load_fragments, optional_backend, resolve_states
from pyar.optimization_request import preflight
from pyar.state.microsolvation import MicrosolvationRunState, MicrosolvationStateError
from pyar.workflows.microsolvation import microsolvate


def build_parser(prog=None):
    parser = argparse.ArgumentParser(
        prog=prog or "pyar microsolvate",
        description="Build a bounded explicit solvent ensemble around the original solute surface.",
    )
    parser.add_argument("solute", help="Original solute XYZ; remains the placement target")
    parser.add_argument("solvent", help="One solvent XYZ molecule to add repeatedly")
    parser.add_argument("--count", type=int, required=True, help="Number of solvent molecules to add")
    add_backend_arguments(parser)
    add_fragment_arguments(parser)
    parser.add_argument("--orientations", "-N", type=int, default=8)
    parser.add_argument("--maximum-number-of-seeds", type=int, default=12)
    parser.add_argument("--site", nargs="+", type=int, metavar="ATOM",
                        help="0-based original-solute atom(s) whose accessible surface is targeted")
    parser.add_argument("--surface-points-per-atom", type=int, default=96,
                        help="Deterministic Fibonacci samples per solute atom")
    parser.add_argument("--probe-radius", type=float, default=1.4,
                        help="Probe expansion in Angstrom for the discretized vdW placement surface")
    parser.add_argument("--shell-tolerance", type=float, default=3.5,
                        help="Maximum solvent-centre distance from the target surface after relaxation (Angstrom)")
    parser.add_argument("--output", default="microsolvation")
    parser.add_argument("--check", action="store_true", help="Validate request and preflight without creating files")
    return parser


def resolve_request(args):
    molecules = resolve_states(
        load_fragments([args.solute, args.solvent]), args.charge, args.multiplicity, args.scftype,
    )
    request = MicrosolvationRequest(
        molecules[0], molecules[1], args.count, args.orientations, optional_backend(args),
        args.maximum_number_of_seeds, None if args.site is None else tuple(args.site),
        args.surface_points_per_atom, args.probe_radius, args.shell_tolerance,
    )
    return request


def main(argv=None, *, prog=None):
    parser = build_parser(prog)
    args = parser.parse_args(argv)
    try:
        request = resolve_request(args)
        surface = build_solute_surface(
            request.solute, points_per_atom=request.surface_points_per_atom,
            probe_radius=request.probe_radius,
        )
        target_surface = surface.restrict(request.site)
        if not len(target_surface.points):
            raise ValueError("The selected solute atoms have no accessible surface points")
        MicrosolvationRunState.load(args.output, request.to_state_dict())
        qc = dict(request.backend_parameters)
        requirements = preflight(
            qc, [request.solute, request.solvent],
            check_example="pyar microsolvate solute.xyz solvent.xyz --count N --backend BACKEND --check",
        ) if qc else []
        prepared_run(args, request=request.to_state_dict(), backend=qc,
                     molecules=[request.solute, request.solvent], requirements=requirements, outputs=[args.output])
        if args.check:
            site = "whole accessible solute surface" if request.site is None else ", ".join(map(str, request.site))
            print("Preflight: microsolvate")
            print(f"Solute: {request.solute.name} ({len(request.solute)} atoms)")
            print(f"Target surface: {site}")
            print(f"Solvent: {request.solvent.name} ({len(request.solvent)} atoms)")
            print(f"Count: {request.count}; orientations: {request.number_of_orientations}")
            print(f"Surface: Fibonacci, {request.surface_points_per_atom} points/atom, probe {request.probe_radius:g} Å")
            print(f"Maximum survivors: {request.maximum_number_of_seeds}")
            print(f"Shell filter: {request.shell_tolerance:g} Å post-optimization; no restraint potential applied")
            print(f"Backend: {qc.get('software', 'none (geometry-only)')}")
            print("Requirements: " + (", ".join(requirements) if requirements else "none"))
            print("Ready to run.\n--check specified; no calculations were performed.")
            return None
        started_run()
        result = microsolvate(request, output=args.output)
    except (ValueError, OSError, ImportError, MicrosolvationStateError, BackendExecutionError) as exc:
        parser.error(str(exc))
    finished_run(result)
    print(
        f"Microsolvation {result.status}.\nCompleted solvents: "
        f"{result.metadata['completed_solvents']}/{request.count}\n"
        f"Selected final structures: {len(result.selected_paths)}\n"
        f"Run directory: {result.run_directory}"
    )
    if result.status not in {"completed", "stopped"}:
        raise SystemExit(1)


def modern_main(argv=None, *, prog=None):
    return main(argv, prog=prog)
