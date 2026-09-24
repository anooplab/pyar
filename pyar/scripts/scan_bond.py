"""Command-line interface for the ORCA-only bond scan."""

import argparse
from pathlib import Path

from pyar.workflows.scan_bond import run_scan_bond


def main(argv=None):
    parser = argparse.ArgumentParser(description="Run an ORCA relaxed scan along a fragment-local bond coordinate.")
    parser.add_argument("inputs", nargs=2, metavar="XYZ", help="fragment A and fragment B XYZ files")
    parser.add_argument("--atoms", nargs=2, type=int, required=True, metavar=("I", "J"),
                        help="0-based atom indices local to fragments A and B")
    parser.add_argument("-N", "--number-of-orientations", type=int, required=True, dest="orientations")
    parser.add_argument("--software", default="orca")
    parser.add_argument("--method", default="BP86")
    parser.add_argument("--basis", default="def2-SVP")
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
    args = parser.parse_args(argv)
    if args.orientations < 1:
        parser.error("orientation count must be at least 1")
    if args.software.lower() != "orca":
        parser.error("scan-bond currently supports only the ORCA backend")
    if not all(Path(path).is_file() for path in args.inputs):
        parser.error("both input XYZ files must exist")
    def pair_value(values, label, index):
        if len(values) == 1: return values[0]
        if len(values) == 2: return values[index]
        parser.error(f"{label} accepts one value or one value per fragment")
    params = {"software": "orca", "method": args.method, "basis": args.basis,
              "nprocs": args.nprocs, "scf_cycles": args.scf_cycles,
              "opt_cycles": args.opt_cycles, "opt_threshold": args.opt_threshold,
              "charge_a": pair_value(args.charge, "charge", 0),
              "charge_b": pair_value(args.charge, "charge", 1),
              "multiplicity_a": pair_value(args.multiplicity, "multiplicity", 0),
              "multiplicity_b": pair_value(args.multiplicity, "multiplicity", 1),
              "scftype_a": pair_value(args.scftype, "scftype", 0),
              "scftype_b": pair_value(args.scftype, "scftype", 1)}
    return run_scan_bond(args.inputs[0], args.inputs[1], args.atoms, args.orientations,
                         params, args.output, args.scan_end, args.scan_step, args.scan_points)


if __name__ == "__main__":
    main()
