"""Modern reaction interface; the canonical workflow owns all reaction science."""

import argparse

from pyar.backend_errors import BackendExecutionError
from pyar.reaction_request import resolve_reaction_request, preflight_reaction, validate_restart
from pyar.state.reaction import ReactionStateError
from pyar.workflows import reaction as reaction_workflow


def build_parser(prog=None):
    parser = argparse.ArgumentParser(prog=prog or 'pyar react', description=__doc__)
    parser.add_argument('input_files', nargs=2, metavar='XYZ')
    parser.add_argument('--backend', help='Explicit physical energy/gradient backend')
    parser.add_argument('--bias-max', type=float, help='Required scientific bias ceiling (kJ/mol)')
    parser.add_argument('--orientations', '-N', type=int, help='Trial orientations (default: 8)')
    parser.add_argument('--charge', nargs='+', type=int)
    parser.add_argument('--multiplicity', nargs='+', type=int)
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
    advanced.add_argument('--scftype', nargs='+', choices=('rhf', 'uhf'))
    for name in ('method', 'basis', 'opt-threshold'):
        advanced.add_argument('--'+name)
    for name in ('nprocs', 'opt-cycles', 'scf-cycles', 'release-retry-limit'):
        advanced.add_argument('--'+name, type=int)
    for name in ('bias-alpha-min', 'bias-alpha-margin', 'bias-alpha-smoothing', 'bias-alpha-epsilon',
                 'bias-scheduled-alpha', 'release-margin-factor', 'release-distance-fraction'):
        advanced.add_argument('--'+name, type=float)
    return parser


def main(argv=None, *, prog=None):
    parser = build_parser(prog)
    args = parser.parse_args(argv)
    try:
        request = resolve_reaction_request(args.input_files, vars(args))
        requirements = preflight_reaction(request)
        validate_restart(request)
    except (ValueError, OSError, ImportError, ReactionStateError) as exc:
        parser.error(f'Preflight failed for react: {exc}\nNo calculations were started.')
    if args.check:
        print('Preflight: react')
        for path, molecule in zip(args.input_files, request.reactants):
            print(f'  {path}: charge {molecule.charge}, multiplicity {molecule.multiplicity}')
        print(f"Backend: {request.qc_params['software']}\n"
              f"Protocol: {request.qc_params['bias_potential']}, {request.qc_params['bias_controller']}, "
              f"bias ceiling {request.bias_max}, {request.orientations} orientations, geomeTRIC/TRIC")
        print('Requirements: ' + ', '.join(requirements))
        print('Ready to run.\n--check specified; no calculations were performed.')
        return
    try:
        result = reaction_workflow.react(*request.reactants, request.bias_min, request.bias_max,
                                         request.orientations, request.qc_params, request.site,
                                         request.proximity_factor)
    except (ValueError, OSError, ReactionStateError, BackendExecutionError) as exc:
        parser.error(str(exc))
    print(f'Reaction search: {result.status}\nProducts found: {len(result.selected_paths)}\n'
          f'Run directory: {result.run_directory}\nState: {result.state_path}')
    for path in result.selected_paths:
        print('  ' + path)
    if result.status.startswith('failed'):
        raise SystemExit(1)
