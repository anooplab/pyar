"""Dedicated CLI for fixed-seed sequential growth."""
import argparse

from pyar.cli import (_normalize_parameter_list, _infer_default_multiplicities,
                      _preflight_cli_requirements,
                      _validate_backend_qc_options)
from pyar.core.molecule import Molecule, parse_xyz
from pyar.growth.request import GrowRequest
from pyar.state.grow import GrowStateError
from pyar.workflows.grow import grow


def build_parser():
    parser = argparse.ArgumentParser(prog="pyar-cli grow",
                                    description="Grow a specific seed by sequential additions of one species; retain every selected stage.")
    parser.add_argument("seed", help="seed XYZ geometry")
    parser.add_argument("monomer", help="repeated addend XYZ geometry")
    parser.add_argument("--count", type=int, required=True)
    parser.add_argument("-N", "--number-of-orientations", type=int, default=8, dest="orientations")
    parser.add_argument("--maximum-number-of-seeds", type=int, default=12)
    parser.add_argument("--software", help="optimization backend; omit for bounded geometry-only growth")
    parser.add_argument("--geometry-optimizer", choices=("native", "geometric"), default="native")
    parser.add_argument("--xtb-model", choices=("gfn2", "gxtb"))
    parser.add_argument("--method")
    parser.add_argument("--basis")
    parser.add_argument("--model")
    parser.add_argument("--gxtb-wrapper")
    parser.add_argument("--nprocs", type=int, default=1)
    parser.add_argument("--opt-cycles", type=int, default=100)
    parser.add_argument("--scf-cycles", type=int, default=1000)
    parser.add_argument("--opt-threshold", choices=("loose", "normal", "tight"), default="normal")
    parser.add_argument("--scf-threshold", default="tight")
    parser.add_argument("--custom-keywords")
    parser.add_argument("-c", "--charge", nargs="+", type=int)
    parser.add_argument("-m", "--multiplicity", nargs="+", type=int)
    parser.add_argument("--scftype", nargs="+")
    parser.add_argument("--site", nargs=2, type=int, metavar=("SEED_INDEX", "MONOMER_INDEX"),
                        help="0-based fragment-local placement sites; seed index refers to the original seed")
    parser.add_argument("--connectivity-policy", choices=("auto", "off", "prefer", "strict"), default="auto")
    parser.add_argument("--selection-feature", default="auto")
    parser.add_argument("--selection-algorithm", default="auto")
    parser.add_argument("--selection-distance", default="euclidean")
    parser.add_argument("--selection-system-type", default="auto")
    parser.add_argument("--output", default="grow")
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        molecules = []
        for path in (args.seed, args.monomer):
            atoms, coordinates, name, title, energy = parse_xyz(path)
            molecules.append(Molecule(atoms, coordinates, name=name, title=title, energy=energy))
        charges = _normalize_parameter_list(args.charge, 0, 2, "Charges")
        multiplicities = (_infer_default_multiplicities(molecules, charges) if args.multiplicity is None
                          else _normalize_parameter_list(args.multiplicity, 1, 2, "Multiplicities"))
        scftypes = _normalize_parameter_list(args.scftype, "rhf", 2, "SCF types")
        for molecule, charge, multiplicity, scftype in zip(molecules, charges, multiplicities, scftypes):
            molecule.charge, molecule.multiplicity, molecule.scftype = charge, multiplicity, scftype
        qc = {key: value for key, value in vars(args).items()
              if key in {"software", "geometry_optimizer", "xtb_model", "method", "basis", "model",
                         "gxtb_wrapper", "nprocs", "opt_cycles", "scf_cycles", "opt_threshold", "scf_threshold", "custom_keywords"}
              and value is not None}
        request = GrowRequest(molecules[0], molecules[1], args.count, args.orientations, qc,
                              args.maximum_number_of_seeds, None if args.site is None else tuple(args.site),
                              args.connectivity_policy, args.selection_feature, args.selection_algorithm,
                              args.selection_distance, args.selection_system_type)
        if args.software:
            supplied = {key for key in ("method", "basis", "model", "custom_keywords") if getattr(args, key) is not None}
            _, unsupported = _validate_backend_qc_options(request.backend_parameters["software"], supplied,
                                                         args.geometry_optimizer)
            if unsupported:
                raise ValueError(f"Unsupported backend options: {', '.join(sorted(unsupported))}")
            _preflight_cli_requirements("grow", request.backend_parameters["software"], args.geometry_optimizer)
            if request.backend_parameters["software"] == "aimnet_2":
                from pyar.backends.aimnet2_assets import validate_aimnet2_runtime_assets
                validate_aimnet2_runtime_assets()
        elif any((args.method, args.basis, args.model, args.xtb_model)):
            raise ValueError("Backend options require --software")
        result = grow(request, output=args.output)
    except (ValueError, OSError, GrowStateError) as exc:
        parser.error(str(exc))
    print(f"Growth {result.status}: {result.metadata['completed_additions']}/{args.count} additions")
    print(f"Run directory: {result.run_directory}")
    for path in result.selected_paths:
        print(path)
    if result.status not in {"completed", "stopped"}:
        raise SystemExit(1)
    return result
