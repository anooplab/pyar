"""Modern fixed-seed growth interface; the workflow owns sequential growth."""
import argparse

from pyar.cli_runtime import prepared_run, started_run, finished_run
from pyar.backend_errors import BackendExecutionError
from pyar.growth.request import GrowRequest
from pyar.modern_workflow import (add_backend_arguments, add_fragment_arguments,
                                 add_selection_arguments, load_fragments, optional_backend,
                                 resolve_states, print_check)
from pyar.optimization_request import preflight
from pyar.state.grow import GrowRunState, GrowStateError
from pyar.workflows.grow import grow


def build_parser(prog=None):
    parser = argparse.ArgumentParser(prog=prog or 'pyar grow', description=__doc__)
    parser.add_argument('seed')
    parser.add_argument('monomer')
    parser.add_argument('--count', type=int, required=True)
    add_backend_arguments(parser)
    add_fragment_arguments(parser)
    add_selection_arguments(parser, maximum_seeds=12)
    parser.add_argument('--output', default='grow')
    return parser


def resolve_request(args):
    molecules = resolve_states(load_fragments([args.seed, args.monomer]),
                               args.charge, args.multiplicity, args.scftype)
    qc = optional_backend(args)
    return GrowRequest(*molecules, args.count, args.orientations, qc, args.maximum_number_of_seeds,
                       None if args.site is None else tuple(args.site), args.connectivity_policy,
                       args.selection_feature, args.selection_algorithm, args.selection_distance,
                       args.selection_system_type)


def main(argv=None, *, prog=None):
    parser = build_parser(prog)
    args = parser.parse_args(argv)
    try:
        request = resolve_request(args)
        GrowRunState.load(args.output, request.to_state_dict())
        qc = request.backend_parameters
        requirements = preflight(qc, [request.seed, request.monomer],
                                 check_example='pyar grow SEED.xyz MONOMER.xyz --count N --backend BACKEND --check') if qc else []
        prepared_run(args, request=request.to_state_dict(), backend=qc,
                     molecules=[request.seed, request.monomer], requirements=requirements, outputs=[args.output])
        if args.check:
            print_check('grow', qc, requirements, [('Seed', args.seed), ('Addend', args.monomer),
                        ('Additions', request.count), ('Orientations', request.number_of_orientations),
                        ('Max survivors', request.maximum_number_of_seeds),
                        ('Electronic states', [(m.charge, m.multiplicity) for m in (request.seed, request.monomer)])])
            return None
        started_run()
        result = grow(request, output=args.output)
    except (ValueError, OSError, ImportError, GrowStateError, BackendExecutionError) as exc:
        parser.error(str(exc))
    finished_run(result)
    print(f'Growth {result.status}.\nCompleted additions: {result.metadata["completed_additions"]}/{request.count}\n'
          f'Selected final structures: {len(result.selected_paths)}\nRun directory: {result.run_directory}')
    if result.status not in {'completed', 'stopped'}:
        raise SystemExit(1)
