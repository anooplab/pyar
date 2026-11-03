"""Backend-independent relaxed bond scans and reaction-path continuation."""
from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from pyar.cli_runtime import prepared_run, started_run, finished_run
import argparse
import os
from pathlib import Path
from pyar.backends.orca_methods import orca_method
from pyar.backend_capabilities import normalize_backend_name, supported_geometry_backends
from pyar.workflows.scan_path import THROUGH_STAGES, validate_continuation
from pyar.workflows.scan_bond import run_scan_bond_request
from pyar.persistence import LockError

if TYPE_CHECKING:
    from pyar.core.molecule import Molecule


class ScanArgumentParser(argparse.ArgumentParser):
    """Validate the one-molecule or two-fragment input shape."""

    def parse_args(self, args=None, namespace=None):
        parsed = super().parse_args(args, namespace)
        if len(parsed.inputs) not in (1, 2):
            self.error('Provide one molecule or two fragment XYZ files')
        return parsed

def build_parser(*, prog=None):
    parser = ScanArgumentParser(prog=prog, description='Relaxed bond scan with optional NEB, TS, frequency, IRC and endpoint validation.')
    parser.add_argument('inputs', nargs='+', metavar='XYZ', help='one molecule for an intramolecular scan, or fragment A and B XYZ files')
    parser.add_argument('--atoms', nargs=2, type=int, required=True, metavar=('I', 'J'), help='0-based atom indices in the molecule, or local to fragments A and B')
    parser.add_argument('-N', '--orientations', type=int, default=None)
    parser.add_argument('--backend', dest='software', type=lambda name: normalize_backend_name(name.lower()), choices=supported_geometry_backends())
    parser.add_argument('--check', action='store_true', help='Validate the complete requested path without calculation')
    parser.add_argument('--method', help='electronic-structure method (ORCA default: BP86)')
    parser.add_argument('--xtb-model', choices=('gxtb', 'gfn2'), default=None, help='standalone xTB model (default: gfn2)')
    parser.add_argument('--basis', help='basis set; required for DFT methods, not used by xTB methods')
    parser.add_argument('--gxtb-wrapper', help='executable ORCA external-method wrapper (oet_gxtb) for --method g-xTB')
    parser.add_argument('--nprocs', type=int, default=1)
    parser.add_argument('--scf-cycles', type=int, default=1000)
    parser.add_argument('--opt-cycles', type=int, default=100)
    parser.add_argument('--opt-threshold', choices=['loose', 'normal', 'tight'], default='normal')
    parser.add_argument('--scan-end', type=float)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--scan-step', type=float)
    group.add_argument('--scan-points', type=int)
    parser.add_argument('-c', '--charge', type=int, nargs='+', default=[0])
    parser.add_argument('-m', '--multiplicity', type=int, nargs='+', default=None)
    parser.add_argument('--total-multiplicity', type=int, help='total multiplicity of the combined A+B system (defaults to the historical coupling heuristic)')
    parser.add_argument('--scftype', nargs='+', default=None)
    parser.add_argument('--output', default='scan_bond')
    parser.add_argument('--recover-stale-lock', action='store_true', help='Explicit stale v2 writer recovery; never breaks an active lock')
    parser.add_argument('--through', choices=THROUGH_STAGES, default='scan', help='last stage to execute; all includes endpoint frequencies (default: scan)')
    parser.add_argument('--characterize-unique-products', action='store_true', help='after scan and free relaxation, run one scan-seeded path per unique product')
    parser.add_argument('--images', type=int, default=11)
    parser.add_argument('--neb-max-cycles', type=int, default=100, dest='max_cycles')
    parser.add_argument('--max-gradient', type=float, default=0.05)
    parser.add_argument('--average-gradient', type=float, default=0.025)
    parser.add_argument('--spring', type=float, default=1.0)
    parser.add_argument('--climb', type=float, default=0.5)
    parser.add_argument('--align', action='store_true')
    parser.add_argument('--interpolation', choices=('linear', 'idpp', 'geodesic'), default='linear')
    parser.add_argument('--idpp-fmax', type=float, default=0.1)
    parser.add_argument('--idpp-steps', type=int, default=100)
    parser.add_argument('--geodesic-tol', type=float, default=0.002)
    parser.add_argument('--geodesic-max-iter', type=int, default=15)
    parser.add_argument('--product-relaxation-fmax', type=float, default=0.05)
    parser.add_argument('--product-relaxation-max-steps', type=int, default=200)
    parser.add_argument('--ts-optimizer', choices=('geometric', 'sella'), default='geometric')
    parser.add_argument('--ts-max-cycles', type=int, default=200)
    parser.add_argument('--ts-fmax', type=float)
    parser.add_argument('--sella-fmax', type=float, default=0.05)
    parser.add_argument('--sella-internal-coordinates', action='store_true')
    parser.add_argument('--irc-max-cycles', type=int, default=200)
    parser.add_argument('--endpoint-max-cycles', type=int, default=300)
    parser.add_argument('--imaginary-frequency-threshold', type=float, default=20.0)
    parser.add_argument('--irc-endpoint-rmsd-tolerance', type=float, default=0.5)
    return parser

def _backend_parameters_from_args(args, intramolecular, parser):
    """Translate parsed CLI settings into the scan workflow's parameter mapping."""

    def pair_value(values, label, index):
        if len(values) == 1:
            return values[0]
        if len(values) == 2:
            return values[index]
        parser.error(f'{label} accepts one value or one value per fragment')
    params = {'software': args.software, 'method': args.method, 'basis': args.basis, 'gxtb_wrapper': args.gxtb_wrapper, 'nprocs': args.nprocs, 'scf_cycles': args.scf_cycles, 'opt_cycles': args.opt_cycles, 'opt_threshold': args.opt_threshold, 'charge_a': pair_value(args.charge, 'charge', 0), 'charge_b': pair_value(args.charge, 'charge', 1), 'multiplicity_a': pair_value(args.multiplicity, 'multiplicity', 0), 'multiplicity_b': pair_value(args.multiplicity, 'multiplicity', 1), 'total_multiplicity': args.total_multiplicity, 'scftype_a': pair_value(args.scftype, 'scftype', 0), 'scftype_b': pair_value(args.scftype, 'scftype', 1)}
    if args.software == 'xtb':
        params['xtb_model'] = args.xtb_model
    if intramolecular:
        for key in ('charge_b', 'multiplicity_b', 'scftype_b'):
            params.pop(key)
    return params

def _scan_invocation(args, parameters, reaction_options):
    """Use one invocation mapping for preflight and execution."""
    return dict(
        input_a=args.inputs[0],
        input_b=None if len(args.inputs) == 1 else args.inputs[1],
        atoms=args.atoms,
        orientations=args.orientations,
        qc_params=parameters,
        output_dir=args.output,
        scan_end=args.scan_end,
        scan_step=args.scan_step,
        scan_points=args.scan_points,
        through=args.through,
        reaction_options=reaction_options,
        characterize_unique_products=args.characterize_unique_products,
    )


@dataclass(frozen=True, slots=True)
class _PreparedScanCLI:
    """Carry validated inputs and preflight evidence to reporting and execution."""

    invocation: dict
    molecules: list[Molecule]
    request_fields: dict
    requirements: list[str]


def _resolve_scan_controls(args, parser):
    """Resolve input-mode defaults and common controls in the established order."""
    if len(args.inputs) not in (1, 2):
        parser.error("Provide one molecule or two fragment XYZ files")
    intramolecular = len(args.inputs) == 1
    if args.orientations is None:
        if intramolecular:
            args.orientations = 1
        else:
            from pyar.data.defualt_parameters import values

            args.orientations = values["how_many_orientations"]
    if intramolecular and args.orientations != 1:
        parser.error(
            "Intramolecular scans preserve the input geometry; --orientations must be 1"
        )
    if args.software is None:
        parser.error(
            "A backend is required.\nExample: pyar scan-bond A.xyz B.xyz --atoms 1 2 --backend xtb"
        )
    if args.xtb_model is not None and args.software != "xtb":
        parser.error("--xtb-model requires --backend xtb")
    if args.software == "xtb":
        args.xtb_model = args.xtb_model or "gfn2"
    if args.orientations < 1:
        parser.error("orientation count must be at least 1")
    if args.nprocs < 1 or args.opt_cycles < 1 or args.scf_cycles < 1:
        parser.error("nprocs, opt-cycles and scf-cycles must be positive")
    return intramolecular


def _resolve_orca_options(args, parser):
    """Canonicalize the existing ORCA method and external-wrapper contract."""
    args.method = args.method or "BP86"
    canonical, is_xtb = orca_method(args.method)
    args.method = canonical
    if is_xtb and args.basis is not None:
        parser.error("--basis is not applicable to ORCA xTB methods")
    if not is_xtb and (not args.basis):
        parser.error("--basis is required for ORCA DFT methods")
    if canonical == "g-xTB":
        if not args.gxtb_wrapper:
            parser.error("--gxtb-wrapper is required for --method g-xTB")
        wrapper_path = Path(args.gxtb_wrapper).expanduser().resolve()
        if not wrapper_path.is_file() or not os.access(wrapper_path, os.X_OK):
            parser.error("--gxtb-wrapper must point to an executable wrapper file")
        args.gxtb_wrapper = str(wrapper_path)
        if args.through != "scan":
            backend_option = "--backend"
            parser.error(
                f"ORCA external g-xTB supports native scanning only; use {backend_option} xtb --xtb-model gxtb for continuation"
            )
    elif args.gxtb_wrapper:
        parser.error("--gxtb-wrapper is only valid with --method g-xTB")


def _validate_backend_options(args, parser):
    if args.software == "orca":
        _resolve_orca_options(args, parser)
    elif args.software == "gaussian":
        if not args.method or not args.basis:
            parser.error("Gaussian requires --method and --basis")
    elif args.software in {"xtb", "aimnet_2"} and (args.method or args.basis):
        parser.error(
            "--method and --basis are only used with ORCA or Gaussian; standalone xTB uses --xtb-model"
        )
    if args.software != "orca" and args.gxtb_wrapper:
        parser.error("--gxtb-wrapper is only used with ORCA")


def _resolve_continuation_options(args, intramolecular, parser):
    from pyar.workflows.scan_path import REACTION_OPTION_NAMES

    reaction_options = {name: getattr(args, name) for name in REACTION_OPTION_NAMES}
    try:
        validate_continuation(args.through, reaction_options)
    except ValueError as exc:
        parser.error(str(exc))
    if args.characterize_unique_products and args.through == "scan":
        parser.error(
            "--characterize-unique-products requires --through neb or a later stage"
        )
    if (
        args.characterize_unique_products or (intramolecular and args.through != "scan")
    ) and args.interpolation == "geodesic":
        parser.error(
            "--characterize-unique-products currently supports linear or idpp path-seed refinement"
        )
    return reaction_options


def _prepare_scan_cli(args, intramolecular, reaction_options, parser):
    """Verify chemistry, restart authority and capabilities before preparation."""
    if not all((Path(path).is_file() for path in args.inputs)):
        parser.error("all input XYZ files must exist")
    from pyar.scan_request import validate_inputs

    try:
        molecules = validate_inputs(args)
    except ValueError as exc:
        parser.error(str(exc))
    params = _backend_parameters_from_args(args, intramolecular, parser)
    invocation = _scan_invocation(args, params, reaction_options)
    from pyar.scan_request import preflight_scan

    try:
        from pyar.workflows.scan_bond import scan_request_fields, validate_scan_restart

        request_fields = scan_request_fields(
            invocation, molecules[0], None if intramolecular else molecules[1]
        )
        validate_scan_restart(args.output, request_fields)
        requirements = preflight_scan(args, params, molecules)
    except (ValueError, ImportError, OSError, RuntimeError) as exc:
        parser.error(
            f"Preflight failed for: scan-bond\n{exc}\nNo calculations were started."
        )
    return _PreparedScanCLI(invocation, molecules, request_fields, requirements)


def _publish_scan_preflight(args, prepared, intramolecular):
    params = prepared.invocation["qc_params"]
    reaction_options = prepared.invocation["reaction_options"]
    molecules, request_fields = prepared.molecules, prepared.request_fields
    requirements = prepared.requirements
    prepared_run(
        args,
        request={
            "scan": vars(args),
            "continuation": reaction_options,
            **(
                {"scan_plan": request_fields["intramolecular_plan"]}
                if intramolecular
                else {}
            ),
        },
        backend=params,
        molecules=molecules,
        requirements=requirements,
        outputs=[args.output],
    )
    print(
        f"Preflight: scan-bond; backend: {args.software}; orientations: {args.orientations}; through: {args.through}"
    )
    if intramolecular:
        plan = request_fields["intramolecular_plan"]
        print(
            f"Intramolecular {plan['scan_direction']}: {args.inputs[0]}; atoms {args.atoms}"
        )
        print(
            f"Distance: {plan['start_distance_angstrom']:.4f} -> {plan['end_distance_angstrom']:.4f} A; {plan['n_points']} points, {plan['step_angstrom']:.4f} A per interval; then free relaxation"
        )
    else:
        print(
            f"Inputs: {args.inputs[0]}, {args.inputs[1]}; bond: A[{args.atoms[0]}] -> B[{args.atoms[1]}]"
        )
    print("Requirements satisfied: " + ", ".join(requirements))


def _execute_scan_cli(args, invocation, parser):
    started_run()
    try:
        result = run_scan_bond_request(
            invocation, **({"recover_stale": True} if args.recover_stale_lock else {})
        )
    except LockError as exc:
        parser.error(str(exc))
    return result


def _publish_scan_result(args, result, parser):
    finished_run(result)
    if isinstance(result, dict) and "status" in result:
        print(
            f"scan-bond: {result['status']} through {args.through}; summaries: {result['output_dir']}"
        )
        if result["status"] != "complete":
            parser.exit(
                1,
                "One or more orientations failed; inspect summary.json for the failed stage and gate.\n",
            )
    return result


def main(argv=None, *, prog=None):
    parser = build_parser(prog=prog)
    args = parser.parse_args(argv)
    intramolecular = _resolve_scan_controls(args, parser)
    _validate_backend_options(args, parser)
    reaction_options = _resolve_continuation_options(args, intramolecular, parser)
    prepared = _prepare_scan_cli(args, intramolecular, reaction_options, parser)
    _publish_scan_preflight(args, prepared, intramolecular)
    if args.check:
        print("Ready to run. --check specified; no calculations were performed.")
        return
    result = _execute_scan_cli(args, prepared.invocation, parser)
    return _publish_scan_result(args, result, parser)


if __name__ == "__main__":
    main()
