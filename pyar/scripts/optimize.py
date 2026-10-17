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
    parser.add_argument('--recovery', choices=['off', 'safe'], default='off',
                        help='Opt-in bounded recovery pilot (default: off)')
    parser.add_argument('--recovery-max-attempts', type=int, default=3)
    parser.add_argument('--recovery-max-cycles', type=int, default=700)
    parser.add_argument('--recovery-max-evaluations', type=int, default=1000)
    parser.add_argument('--recovery-wall-seconds', type=float, default=3600.)
    parser.add_argument('--recovery-resume', action='store_true')
    parser.add_argument('--recovery-recover-lock', action='store_true',
                        help='Explicitly recover stale lock metadata after obtaining OS lock')
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
        if args.recovery == 'safe':
            from pyar.recovery.policy import RecoveryPolicy
            from pyar.recovery.controller import preflight_recovery
            policy = RecoveryPolicy('safe', args.recovery_max_attempts, args.recovery_max_cycles,
                                    args.recovery_max_evaluations, args.recovery_wall_seconds)
            preflight_recovery(molecules, settings, policy, resume=args.recovery_resume)
        elif args.recovery_resume or args.recovery_recover_lock:
            raise ValueError('Recovery resume/lock recovery requires --recovery safe')
    except (ValueError, OSError, ImportError) as exc:
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
    recoveries = []
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
            if args.recovery == 'safe':
                from pyar.recovery.controller import run_recovery
                from pyar.persistence import LockError
                try:
                    recovery = run_recovery(molecule, molecule_settings, policy,
                        resume=args.recovery_resume, recover_lock=args.recovery_recover_lock)
                except (ValueError, OSError, LockError) as exc:
                    raise BackendExecutionError(str(exc)) from exc
                recoveries.append({'input': molecule.name, **recovery.to_dict()})
                if not recovery.result.converged:
                    raise BackendExecutionError('SAFE recovery stopped: ' + recovery.stop_reason)
                import numpy as np
                frame = recovery.result.accepted_frame
                molecule.coordinates = np.array(frame.coordinates_angstrom)
                molecule.energy = frame.physical_energy_hartree
                status = True
            else:
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
    if args.recovery == 'safe':
        report['recovery'] = recoveries
    finished_run(report)
    print(f'Optimized {len(results)} of {len(molecules)} structures ({status}).')
    if failures:
        for failure in failures:
            print(f"Failed {failure['input']}: {failure['error']}", file=sys.stderr)
        parser.exit(1, f'Optimization {status}: {len(failures)} input(s) failed.\n')


if __name__ == '__main__':
    main()
