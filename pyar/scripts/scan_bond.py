"""Backend-independent relaxed bond scans and reaction-path continuation."""

import argparse
import os
from pathlib import Path

from pyar.backends.orca_methods import orca_method
from pyar.backend_capabilities import supported_geometry_backends
from pyar.workflows.scan_path import THROUGH_STAGES, validate_continuation
from pyar.workflows.scan_bond import run_scan_bond


def main(argv=None):
    parser = argparse.ArgumentParser(description="Relaxed bond scan with optional NEB, TS, frequency, IRC and endpoint validation.")
    parser.add_argument("inputs", nargs=2, metavar="XYZ", help="fragment A and fragment B XYZ files")
    parser.add_argument("--atoms", nargs=2, type=int, required=True, metavar=("I", "J"),
                        help="0-based atom indices local to fragments A and B")
    parser.add_argument("-N", "--number-of-orientations", type=int, required=True, dest="orientations")
    parser.add_argument("--software", type=str.lower, choices=supported_geometry_backends(), default="orca")
    parser.add_argument("--method", help="electronic-structure method (ORCA default: BP86)")
    parser.add_argument("--xtb-model", choices=("gxtb", "gfn2"), default="gfn2",
                        help="standalone xTB model (default: gfn2)")
    parser.add_argument("--basis", help="basis set; required for DFT methods, not used by xTB methods")
    parser.add_argument("--gxtb-wrapper", help="executable ORCA external-method wrapper (oet_gxtb) for --method g-xTB")
    parser.add_argument("--nprocs", type=int, default=1)
    parser.add_argument("--scf-cycles", type=int, default=1000)
    parser.add_argument("--opt-cycles", type=int, default=100)
    parser.add_argument("--opt-threshold", choices=["loose", "normal", "tight"], default="normal")
    parser.add_argument("--scan-end", type=float)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--scan-step", type=float)
    group.add_argument("--scan-points", type=int)
    parser.add_argument("-c", "--charge", type=int, nargs="+", default=[0])
    parser.add_argument("-m", "--multiplicity", type=int, nargs="+", default=[1])
    parser.add_argument("--scftype", nargs="+", default=["rhf"])
    parser.add_argument("--output", default="scan_bond")
    parser.add_argument("--through", choices=THROUGH_STAGES, default="scan",
                        help="last stage to execute; all includes endpoint frequencies (default: scan)")
    parser.add_argument("--images", type=int, default=11)
    parser.add_argument("--neb-max-cycles", type=int, default=100, dest="max_cycles")
    parser.add_argument("--max-gradient", type=float, default=0.05)
    parser.add_argument("--average-gradient", type=float, default=0.025)
    parser.add_argument("--spring", type=float, default=1.0)
    parser.add_argument("--climb", type=float, default=0.5)
    parser.add_argument("--align", action="store_true")
    parser.add_argument("--interpolation", choices=("linear", "idpp", "geodesic"), default="linear")
    parser.add_argument("--idpp-fmax", type=float, default=0.1)
    parser.add_argument("--idpp-steps", type=int, default=100)
    parser.add_argument("--geodesic-tol", type=float, default=0.002)
    parser.add_argument("--geodesic-max-iter", type=int, default=15)
    parser.add_argument("--product-relaxation-fmax", type=float, default=0.05)
    parser.add_argument("--product-relaxation-max-steps", type=int, default=200)
    parser.add_argument("--ts-optimizer", choices=("geometric", "sella"), default="geometric")
    parser.add_argument("--ts-max-cycles", type=int, default=200)
    parser.add_argument("--ts-fmax", type=float)
    parser.add_argument("--sella-fmax", type=float, default=0.05)
    parser.add_argument("--sella-internal-coordinates", action="store_true")
    parser.add_argument("--irc-max-cycles", type=int, default=200)
    parser.add_argument("--endpoint-max-cycles", type=int, default=300)
    parser.add_argument("--imaginary-frequency-threshold", type=float, default=20.0)
    parser.add_argument("--irc-endpoint-rmsd-tolerance", type=float, default=0.5)
    args = parser.parse_args(argv)
    if args.orientations < 1:
        parser.error("orientation count must be at least 1")
    if args.nprocs < 1 or args.opt_cycles < 1 or args.scf_cycles < 1:
        parser.error("nprocs, opt-cycles and scf-cycles must be positive")
    if args.software == "orca":
        args.method = args.method or "BP86"
        canonical, is_xtb = orca_method(args.method)
        args.method = canonical
        if is_xtb and args.basis is not None:
            parser.error("--basis is not applicable to ORCA xTB methods")
        if not is_xtb and not args.basis:
            parser.error("--basis is required for ORCA DFT methods")
        if canonical == "g-xTB":
            if not args.gxtb_wrapper:
                parser.error("--gxtb-wrapper is required for --method g-xTB")
            wrapper_path = Path(args.gxtb_wrapper).expanduser().resolve()
            if not wrapper_path.is_file() or not os.access(wrapper_path, os.X_OK):
                parser.error("--gxtb-wrapper must point to an executable wrapper file")
            args.gxtb_wrapper = str(wrapper_path)
            if args.through != "scan":
                parser.error("ORCA external g-xTB supports native scanning only; use --software xtb --xtb-model gxtb for continuation")
        elif args.gxtb_wrapper:
            parser.error("--gxtb-wrapper is only valid with --method g-xTB")
    elif args.software == "gaussian":
        if not args.method or not args.basis:
            parser.error("Gaussian requires --method and --basis")
    elif args.software in {"xtb", "aimnet_2"} and (args.method or args.basis):
        parser.error("--method and --basis are only used with ORCA or Gaussian; standalone xTB uses --xtb-model")
    if args.software != "orca" and args.gxtb_wrapper:
        parser.error("--gxtb-wrapper is only used with ORCA")
    from pyar.workflows.scan_path import REACTION_OPTION_NAMES
    reaction_options = {name: getattr(args, name) for name in REACTION_OPTION_NAMES}
    try:
        validate_continuation(args.through, reaction_options)
    except ValueError as exc:
        parser.error(str(exc))
    if not all(Path(path).is_file() for path in args.inputs):
        parser.error("both input XYZ files must exist")
    def pair_value(values, label, index):
        if len(values) == 1: return values[0]
        if len(values) == 2: return values[index]
        parser.error(f"{label} accepts one value or one value per fragment")
    params = {"software": args.software, "method": args.method, "basis": args.basis,
              "gxtb_wrapper": args.gxtb_wrapper,
              "nprocs": args.nprocs, "scf_cycles": args.scf_cycles,
              "opt_cycles": args.opt_cycles, "opt_threshold": args.opt_threshold,
              "charge_a": pair_value(args.charge, "charge", 0),
              "charge_b": pair_value(args.charge, "charge", 1),
              "multiplicity_a": pair_value(args.multiplicity, "multiplicity", 0),
              "multiplicity_b": pair_value(args.multiplicity, "multiplicity", 1),
              "scftype_a": pair_value(args.scftype, "scftype", 0),
              "scftype_b": pair_value(args.scftype, "scftype", 1)}
    if args.software == "xtb":
        params["xtb_model"] = args.xtb_model
    result = run_scan_bond(args.inputs[0], args.inputs[1], args.atoms, args.orientations,
                         params, args.output, args.scan_end, args.scan_step, args.scan_points,
                         args.through, reaction_options)
    if isinstance(result, dict) and "status" in result:
        print(f"scan-bond: {result['status']} through {args.through}; summaries: {result['output_dir']}")
        if result["status"] != "complete":
            parser.exit(1, "One or more orientations failed; inspect summary.json for the failed stage and gate.\n")
    return result


if __name__ == "__main__":
    main()
