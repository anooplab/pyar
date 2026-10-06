"""Modern reaction interface; the canonical workflow owns all reaction science."""

import argparse
import json
from dataclasses import replace
from pathlib import Path

from pyar.cli_runtime import prepared_run, started_run, finished_run
from pyar.backend_errors import BackendExecutionError
from pyar.reaction_request import resolve_reaction_request, preflight_reaction, validate_restart
from pyar.state.reaction import ReactionStateError
from pyar.workflows import reaction as reaction_workflow
from pyar.workflows.reaction_characterization import (
    characterize_reaction, preflight_path_request, validate_characterization_restart,
)
from pyar.workflows.scan_path import REACTION_OPTION_NAMES


def build_parser(prog=None):
    parser = argparse.ArgumentParser(prog=prog or 'pyar react', description=__doc__)
    parser.add_argument('input_files', nargs=2, metavar='XYZ')
    parser.add_argument('--backend', help='Explicit physical energy/gradient backend')
    parser.add_argument('--bias-max', type=float, help='Required scientific bias ceiling (kJ/mol)')
    parser.add_argument('--orientations', '-N', type=int, help='Trial orientations (default: 8)')
    parser.add_argument('--charge', nargs='+', type=int)
    parser.add_argument('--multiplicity', nargs='+', type=int)
    parser.add_argument('--through', choices=('react', 'neb', 'ts', 'frequency', 'irc', 'endpoints',
                                               'endpoint-frequency', 'all'), default='react',
                        help='Continue each unique product through unbiased path validation (default: react only)')
    parser.add_argument('--bias-potential', choices=('afir', 'softmin'), help='Default: afir')
    parser.add_argument('--bias-controller', choices=('adaptive', 'fixed', 'scheduled'), help='Default: adaptive')
    parser.add_argument('--check', action='store_true', help='Validate both phases and restart; run no calculations')
    advanced = parser.add_argument_group('Advanced reaction/backend overrides')
    advanced.add_argument('--bias-min', type=float, help='Required for fixed/scheduled; ignored for adaptive')
    advanced.add_argument('--softmin-beta', type=float)
    advanced.add_argument('--geometry-optimizer', choices=('geometric', 'native'))
    advanced.add_argument('--site', type=int, nargs=2, metavar=('I', 'J'),
                          help='0-based fragment-local atom index in A and B')
    advanced.add_argument('--proximity-factor', type=float)
    advanced.add_argument('--pathway-ts-source', choices=('highest-backend-energy', 'first-topology-change',
                              'pre-product', 'highest-total-energy'), default='highest-backend-energy',
                          help='Trace geometry used only as a NEB waypoint/TS guess')
    advanced.add_argument('--scftype', nargs='+', choices=('rhf', 'uhf'))
    for name in ('method', 'basis', 'opt-threshold'):
        advanced.add_argument('--'+name)
    for name in ('nprocs', 'opt-cycles', 'scf-cycles', 'release-retry-limit'):
        advanced.add_argument('--'+name, type=int)
    for name in ('bias-alpha-min', 'bias-alpha-margin', 'bias-alpha-smoothing', 'bias-alpha-epsilon',
                 'bias-scheduled-alpha', 'release-margin-factor', 'release-distance-fraction'):
        advanced.add_argument('--'+name, type=float)
    path = parser.add_argument_group('Advanced unbiased path continuation options')
    for name, kind in (('images', int), ('max-cycles', int), ('max-gradient', float),
                       ('average-gradient', float), ('spring', float), ('climb', float),
                       ('interpolation', str), ('idpp-fmax', float), ('idpp-steps', int),
                       ('geodesic-tol', float), ('geodesic-max-iter', int),
                       ('product-relaxation-fmax', float), ('product-relaxation-max-steps', int),
                       ('ts-max-cycles', int), ('ts-fmax', float), ('sella-fmax', float),
                       ('irc-max-cycles', int), ('endpoint-max-cycles', int),
                       ('imaginary-frequency-threshold', float), ('irc-endpoint-rmsd-tolerance', float)):
        kwargs = {'type': kind}
        if name == 'interpolation':
            kwargs['choices'] = ('linear', 'idpp', 'geodesic')
        if name == 'max-cycles':
            path.add_argument('--neb-max-cycles', '--max-cycles', dest='max_cycles', default=None, **kwargs)
        else:
            path.add_argument('--' + name, default=None, **kwargs)
    path.add_argument('--ts-optimizer', choices=('geometric', 'sella'))
    path.add_argument('--align', action=argparse.BooleanOptionalAction, default=None)
    path.add_argument('--sella-internal-coordinates', action=argparse.BooleanOptionalAction, default=None)
    return parser


def main(argv=None, *, prog=None):
    parser = build_parser(prog)
    args = parser.parse_args(argv)
    try:
        request = resolve_reaction_request(args.input_files, vars(args))
        path_options = {key: value for key, value in vars(args).items()
                        if key in REACTION_OPTION_NAMES and value is not None}
        requirements = sorted(set(preflight_reaction(request) +
                                  preflight_path_request(request, args.through, path_options)))
        reaction_state_path = Path.cwd() / 'reaction' / 'state.json'
        completed_state = None
        if args.through != 'react' and reaction_state_path.is_file():
            try:
                raw_status = json.loads(reaction_state_path.read_text(encoding='utf-8')).get('status')
            except (OSError, ValueError):
                raw_status = None
            if raw_status in {'completed', 'completed_products_found', 'completed_no_products',
                              'completed_no_candidates'}:
                from pyar.state.reaction import ReactionRunState
                completed_state = ReactionRunState.read_completed(Path.cwd(), request.restart_request)
                validate_characterization_restart(Path.cwd(), request, args.through,
                                                   dict(path_options, through=args.through,
                                                        pathway_ts_source=args.pathway_ts_source))
            else:
                validate_restart(request)
        else:
            validate_restart(request)
    except (ValueError, OSError, ImportError, ReactionStateError) as exc:
        parser.error(f'Preflight failed for react: {exc}\nNo calculations were started.')
    prepared_run(args, request=request.restart_request, backend=request.qc_params,
                 molecules=request.reactants, requirements=requirements, outputs=['reaction'],
                 options={'orientations': request.orientations, 'bias_min': request.bias_min,
                          'bias_max': request.bias_max, 'proximity_factor': request.proximity_factor,
                          'through': args.through, 'pathway_ts_source': args.pathway_ts_source})
    if args.check:
        print('Preflight: react')
        for path, molecule in zip(args.input_files, request.reactants):
            print(f'  {path}: charge {molecule.charge}, multiplicity {molecule.multiplicity}')
        print(f"Backend: {request.qc_params['software']}\n"
              f"Protocol: {request.qc_params['bias_potential']}, {request.qc_params['bias_controller']}, "
              f"bias ceiling {request.bias_max}, {request.orientations} orientations, geomeTRIC/TRIC")
        if args.through != 'react':
            print(f"Continuation: through {args.through}\n"
                  f"Pathway policy: one selected route per unique product\n"
                  f"TS waypoint: {args.pathway_ts_source} (initializer only)\n"
                  "Planned stages: discovery, endpoint relaxation, NEB, TS, frequency, IRC, endpoint validation")
        print('Requirements: ' + ', '.join(requirements))
        print('Ready to run.\n--check specified; no calculations were performed.')
        return
    try:
        started_run()
        if completed_state is None:
            result = reaction_workflow.react(*request.reactants, request.bias_min, request.bias_max,
                                             request.orientations, request.qc_params, request.site,
                                             request.proximity_factor)
        else:
            products = completed_state.data.get('products', [])
            result = reaction_workflow.ReactionResult(
                workflow='reaction', status=completed_state.data['status'],
                run_directory=str(Path.cwd() / 'reaction'),
                state_path=str(reaction_state_path),
                selected_paths=tuple(str(Path.cwd() / 'reaction' / item['path']) for item in products),
                metadata={'products': tuple(products), 'resumed_discovery': True})
        if args.through != 'react' and not result.status.startswith('failed'):
            characterization = characterize_reaction(
                Path.cwd(), request, args.through, dict(path_options, pathway_ts_source=args.pathway_ts_source))
            result = replace(result, status=(result.status if characterization['status'] in
                              {'complete', 'complete_no_products'} else 'completed_path_failed'),
                             metadata={**dict(result.metadata), 'characterization': characterization})
    except (ValueError, OSError, ReactionStateError, BackendExecutionError) as exc:
        parser.error(str(exc))
    finished_run(result)
    print(f'Reaction search: {result.status}\nProducts found: {len(result.selected_paths)}\n'
          f'Run directory: {result.run_directory}\nState: {result.state_path}')
    for path in result.selected_paths:
        print('  ' + path)
    if result.status.startswith('failed'):
        raise SystemExit(1)
    characterization = result.metadata.get('characterization')
    if characterization is not None:
        print(f"Path characterization: {characterization['status']} (through {args.through})")
        for pathway in characterization.get('selected_routes', []):
            print(f"  {pathway['product']} {pathway['pathway']}: {pathway['status']} — {pathway['directory']}")
        for failure in characterization.get('preparation_failures', []):
            print(f"  {failure['product']}: {failure['status']} — {failure['error']}")
        if characterization['status'] == 'failed':
            raise SystemExit(1)
