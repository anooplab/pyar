"""Modern bulk minimum optimization; computational execution stays in optimiser."""

import argparse
import sys

from pyar.cli_runtime import prepared_run, started_run, finished_run
from pyar import optimiser
from pyar.backend_errors import BackendExecutionError
from pyar.optimization_request import load_inputs, preflight, resolve_settings


def build_parser(prog=None):
    parser = argparse.ArgumentParser(prog=prog, description=(
        "Optimize each supplied XYZ independently using the same computational settings."
    ))
    parser.add_argument('input_files', nargs='+', metavar='XYZ')
    parser.add_argument('--backend', help='Explicit computational backend, e.g. xtb or orca')
    parser.add_argument('--check', action='store_true', help='Validate the entire request without calculation')
    parser.add_argument('-c', '--charge', type=int, default=0)
    parser.add_argument('-m', '--multiplicity', type=int, help='Default: infer singlet/doublet from electron parity')
    parser.add_argument('--scftype', choices=['rhf', 'uhf'])
    parser.add_argument('--geometry-optimizer', choices=['native', 'geometric'],
                        help='Default: native, or geomeTRIC for Gaussian')
    parser.add_argument('--opt-target', choices=['minimum', 'ts'], default='minimum',
                        help='Only minimum optimization is supported by this command')
    for option in ('method', 'basis', 'custom-keywords'):
        parser.add_argument('--' + option)
    for option in ('nprocs', 'opt-cycles', 'scf-cycles'):
        parser.add_argument('--' + option, type=int)
    for option in ('opt-threshold', 'scf-threshold'):
        parser.add_argument('--' + option, choices=['loose', 'normal', 'tight'])
    return parser


def main(argv=None, *, prog=None):
    parser = build_parser(prog)
    args = parser.parse_args(argv)
    try:
        settings = resolve_settings(args)
        molecules = load_inputs(args)
        requirements = preflight(settings, molecules)
    except (ValueError, FileNotFoundError, ImportError) as exc:
        parser.error(f'Preflight failed for: optimize\n{exc}\nNo calculations were started.')
    prepared_run(args, request={'inputs': args.input_files}, backend=settings,
                 molecules=molecules, requirements=requirements, state_lists=False,
                 outputs=[f'job_{m.name}' for m in molecules])
    print('Preflight: optimize')
    print(f'Inputs: {len(molecules)} XYZ files; backend: {settings["software"]}; '
          f'optimizer: {settings["geometry_optimizer"]}')
    print(f'Charge: {args.charge}; multiplicities: '
          + ', '.join(str(mol.multiplicity) for mol in molecules))
    print(f'Method: {settings.get("method") or settings.get("xtb_model") or "backend standard"}; '
          f'basis: {settings.get("basis") or "not applicable"}; nprocs: {settings["nprocs"]}')
    print('Requirements satisfied: ' + ', '.join(requirements))
    if args.check:
        print('Ready to run. --check specified; no calculations were performed.')
        return
    started_run()
    results = []
    failures = []
    for molecule in molecules:
        molecule_settings = settings
        if settings.get('_electronic_state_from_molecule'):
            molecule_settings = {
                **settings,
                'charge': molecule.charge,
                'multiplicity': molecule.multiplicity,
                'scftype': molecule.scftype,
                'xtb_unpaired_electrons': molecule.multiplicity - 1,
            }
        try:
            status = optimiser.optimise(molecule, molecule_settings)
        except (BackendExecutionError, FileNotFoundError) as exc:
            failures.append({'input': molecule.name, 'error': str(exc)})
            continue
        if optimiser.is_usable(status):
            results.append(molecule)
        else:
            failures.append({'input': molecule.name, 'error': f'Optimization returned unusable status: {status!r}'})

    if len(results) == len(molecules):
        status = 'completed'
    elif results:
        status = 'partial_failed'
    else:
        status = 'failed'
    report = {
        'workflow': 'optimize', 'status': status,
        'input_count': len(molecules), 'optimized_count': len(results),
        'failed_count': len(failures), 'failures': failures,
    }
    finished_run(report)
    print(f'Optimized {len(results)} of {len(molecules)} structures ({status}).')
    if failures:
        for failure in failures:
            print(f"Failed {failure['input']}: {failure['error']}", file=sys.stderr)
        parser.exit(1, f'Optimization {status}: {len(failures)} input(s) failed.\n')


if __name__ == '__main__':
    main()
