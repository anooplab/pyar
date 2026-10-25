"""Batch single-point energy, gradient and minimum qualification commands."""

import argparse
import sys

from pyar.cli_runtime import prepared_run, started_run, finished_run
from pyar.singlepoint.request import resolve, preflight
from pyar.singlepoint.service import run_batch


def build_parser(command="energy", *, prog=None):
    command = {"gradient": "engrad", "engra": "engrad"}.get(command, command)
    parser = argparse.ArgumentParser(
        prog=prog,
        description=f"Calculate {command} for each fixed XYZ geometry independently.",
    )
    parser.add_argument("input_files", nargs="+", metavar="XYZ")
    parser.add_argument(
        "--backend",
        default="xtb",
        help="Default: xtb; also orca, gaussian, turbomole, aimnet_2, ani",
    )
    parser.add_argument(
        "--output", help=f"New output root (default: {command})", default=command
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="Preflight all inputs without calculations or output folders",
    )
    parser.add_argument("-c", "--charge", type=int, default=0)
    parser.add_argument(
        "-m",
        "--multiplicity",
        type=int,
        help="Default: singlet/doublet from electron parity",
    )
    for name in ("method", "basis", "xtb-model", "model", "solvent", "solvent-model"):
        parser.add_argument("--" + name)
    parser.add_argument(
        "--dispersion",
        choices=("none", "D3BJ", "D3ZERO", "D4"),
        default="none",
        help="Additional dispersion correction; default: none (intrinsic method terms remain)",
    )
    parser.add_argument("--nprocs", type=int, default=1)
    parser.add_argument(
        "--memory-mb",
        type=int,
        help="ORCA: per process; Gaussian: total (default: 2000 MB)",
    )
    parser.add_argument(
        "--scf-cycles", type=int, help="Maximum electronic iterations (default: 250)"
    )
    parser.add_argument(
        "--accuracy", type=float, help="xTB accuracy parameter (default: 0.1)"
    )
    if command == "frequency":
        parser.add_argument(
            "--step",
            type=float,
            default=0.005,
            help="Central-difference displacement in bohr (default: 0.005)",
        )
        parser.add_argument(
            "--max-gradient",
            type=float,
            default=4.5e-4,
            help="Stationarity maximum per-atom gradient, Eh/bohr",
        )
        parser.add_argument(
            "--rms-gradient",
            type=float,
            default=3e-4,
            help="Stationarity RMS per-atom gradient, Eh/bohr",
        )
        parser.add_argument(
            "--imaginary-threshold",
            type=float,
            default=20.0,
            help="Significant imaginary frequency magnitude, cm^-1 (default: 20)",
        )
    return parser


def frequency_controls(args):
    names = ("step", "max_gradient", "rms_gradient", "imaginary_threshold")
    return {name: getattr(args, name) for name in names if hasattr(args, name)}


def report_warnings(report):
    for result in report["results"]:
        if result.get("warning"):
            print(f"Warning [{result['name']}]: {result['warning']}", file=sys.stderr)
        if result.get("error"):
            print(f"Failed [{result['name']}]: {result['error']}", file=sys.stderr)


def prepare(parser, args, command):
    try:
        settings, molecules, root = resolve(args, command)
        requirements = preflight(settings)
    except (ValueError, OSError, ImportError) as exc:
        parser.error(f"Preflight failed: {exc}\nNo calculations were started.")
    controls = frequency_controls(args)
    options = {key: value for key, value in settings.items() if key in vars(args)}
    prepared_run(
        args,
        request={"inputs": args.input_files, "frequency_controls": controls},
        backend=settings,
        molecules=molecules,
        requirements=requirements,
        outputs=[root],
        state_lists=False,
        options=options,
    )
    print(
        f"Preflight: {command}; {len(molecules)} geometries; backend: {settings['software']}"
    )
    if command == "frequency":
        print(
            "Central-difference gradient evaluations: "
            + str(sum(6 * m.number_of_atoms + 1 for m in molecules))
        )
    return settings, molecules, root, controls


def main(argv=None, *, prog=None):
    command = (prog or "pyar energy").split()[-1]
    command = {"gradient": "engrad", "engra": "engrad"}.get(command, command)
    parser = build_parser(command, prog=prog)
    args = parser.parse_args(argv)
    settings, molecules, root, controls = prepare(parser, args, command)
    if args.check:
        print("Ready. No calculations or output folders were created.")
        return
    started_run()
    report = run_batch(command, args.input_files, molecules, settings, root, controls)
    finished_run(report)
    report_warnings(report)
    if report["status"] != "completed":
        parser.exit(1, f"{command}: {report['status']}; see {root / 'summary.json'}\n")
    return report


if __name__ == "__main__":
    main()
