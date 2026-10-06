"""Modern composition search; no input is a permanently privileged seed."""
import argparse
from pathlib import Path

from pyar.aggregation.request import AggregateRequest
from pyar.cli_runtime import prepared_run, started_run, finished_run
from pyar.backend_errors import BackendExecutionError
from pyar.data.defualt_parameters import values
from pyar.modern_workflow import (add_backend_arguments, add_fragment_arguments,
                                 add_selection_arguments, load_fragments, optional_backend,
                                 resolve_states, print_check)
from pyar.optimization_request import preflight
from pyar.state.aggregate import AggregateRunState, AggregateStateError
from pyar.workflows.aggregate import aggregate
from pyar.growth.service import expand_formula_to_aggregate_inputs, read_old_path


def build_parser(prog=None):
    parser = argparse.ArgumentParser(prog=prog or 'pyar aggregate', description=__doc__)
    parser.add_argument('inputs', nargs='*', help='XYZ files or existing element/formula input specifications')
    parser.add_argument('--size', '--aggregate-size', nargs='+', type=int)
    parser.add_argument('--formula', help='Expand a composition formula into atomic building blocks, not bonds')
    add_backend_arguments(parser)
    add_fragment_arguments(parser)
    add_selection_arguments(parser, maximum_seeds=values['maximum_number_of_seeds'])
    parser.add_argument('--number-of-pathways', type=int, default=values['number_of_pathways'])
    parser.add_argument('--first-pathway', type=int, default=values['first_pathway'])
    return parser


def resolve_request(args):
    specs, sizes = args.inputs, args.size
    if args.formula:
        if specs or sizes is not None:
            raise ValueError('--formula cannot be combined with positional inputs or --size')
        specs, sizes = expand_formula_to_aggregate_inputs(args.formula)
    if not specs:
        raise ValueError('Provide building-block inputs or --formula')
    if sizes is None:
        if len(specs) == 1:
            raise ValueError('A single building block requires --size (for example --size 6)')
        sizes = [1] * len(specs)
    if len(sizes) != len(specs) or any(size < 1 for size in sizes) or sum(sizes) < 2:
        raise ValueError('--size requires one positive integer per input and at least two units in total')
    if args.orientations < 1:
        raise ValueError('--orientations must be positive')
    molecules = resolve_states(load_fragments(specs, allow_formula=True), args.charge, args.multiplicity, args.scftype)
    site = args.site
    if site is not None:
        if len(molecules) != 2 or any(not 0 <= index < len(molecule) for index, molecule in zip(site, molecules)):
            raise ValueError('--site requires valid 0-based indices in two input fragments')
        # Preserve the legacy aggregate site's absolute-index contract.
        site = [site[0], len(molecules[0]) + site[1]]
    return AggregateRequest.from_options(molecules, sizes, args.orientations, optional_backend(args),
                                         args.maximum_number_of_seeds, args.first_pathway,
                                         args.number_of_pathways, site, args.connectivity_policy,
                                         args.selection_feature, args.selection_algorithm,
                                         args.selection_distance, args.selection_system_type)


def validate_restart(request):
    state = AggregateRunState.load(Path.cwd(), request.to_state_dict())
    directory = Path('aggregates')
    if directory.exists() and not directory.is_dir():
        raise ValueError('aggregates must be a directory')
    if state is None and directory.is_dir() and any(directory.iterdir()) and not read_old_path():
        raise AggregateStateError('Existing aggregates directory has no resumable state; start in a new directory.')


def main(argv=None, *, prog=None):
    parser = build_parser(prog)
    args = parser.parse_args(argv)
    try:
        request = resolve_request(args)
        validate_restart(request)
        qc = dict(request.backend_parameters)
        requirements = preflight(qc, request.molecules,
                                 check_example='pyar aggregate INPUT.xyz --size N --backend BACKEND --check') if qc else []
        prepared_run(args, request=request.to_state_dict(), backend=qc,
                     molecules=request.molecules, requirements=requirements, outputs=['aggregates'])
        if args.check:
            print_check('aggregate', qc, requirements, [('Fragments (count, charge, multiplicity)',
                        [(m.name, n, m.charge, m.multiplicity) for m, n in zip(request.molecules, request.aggregate_sizes)]),
                        ('Orientations', request.orientations), ('Pathways', request.number_of_pathways),
                        ('Max survivors', request.maximum_number_of_seeds), ('Connectivity', request.connectivity_policy)])
            return
        started_run()
        result = aggregate(list(request.molecules), list(request.aggregate_sizes), request.orientations,
                           qc, request.maximum_number_of_seeds, request.first_pathway,
                           request.number_of_pathways, request.site, request.connectivity_policy,
                           request.selection_feature, request.selection_algorithm,
                           request.selection_distance, request.selection_system_type)
    except (ValueError, OSError, ImportError, AggregateStateError, BackendExecutionError) as exc:
        parser.error(str(exc))
    finished_run(result)
    print(f'Aggregation {result.status}.\nSelected structures: {len(result.selected_paths)}\nRun directory: {result.run_directory}')
    if result.status.startswith('failed'):
        raise SystemExit(1)
