#!/usr/bin/env python3
"""Command-line entrypoint for RDKit conformational search."""

from __future__ import annotations

import argparse
import logging
import sys

from pyar.data import defualt_parameters
from pyar.workflows.conformer import ConformerWorkflowError, conformer_search

logger = logging.getLogger("pyar-conformer")


def argument_parse(argv=None):
    """Parse conformer-search CLI arguments."""
    parser = argparse.ArgumentParser(
        prog="pyar-cli conformer",
        description="Generate and optionally refine RDKit conformers.",
    )
    parser.add_argument("input", help="SMILES string, SDF/MOL file, or XYZ file")
    parser.add_argument(
        "--input-format",
        choices=["auto", "smiles", "sdf", "mol", "xyz"],
        default="auto",
        help="Input format; auto treats existing .xyz/.sdf/.mol paths as files and other input as SMILES.",
    )
    parser.add_argument("--num-conformers", type=int, default=100)
    parser.add_argument("--top-n", type=int, default=10)
    parser.add_argument("--backend-top-n", type=int)
    parser.add_argument("--num-seeds", type=int, default=1)
    parser.add_argument("--diversity-fraction", type=float, default=0.5)
    parser.add_argument(
        "--rms-threshold",
        "--prune-rms-threshold",
        dest="rms_threshold",
        type=float,
        default=0.5,
        help="RDKit greedy prune RMS threshold; lower values keep more embedded conformers.",
    )
    parser.add_argument(
        "--use-random-coords",
        action="store_true",
        help="Start RDKit embedding from random coordinates instead of distance geometry eigenvectors.",
    )
    parser.add_argument("--force-field", choices=["auto", "mmff", "uff"], default="auto")
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--num-threads", type=int, default=0)
    parser.add_argument("--max-iterations", type=int, default=200)

    molecule_group = parser.add_argument_group("molecule")
    molecule_group.add_argument("-c", "--charge", type=int)
    molecule_group.add_argument("-m", "--multiplicity", type=int, default=1)
    molecule_group.add_argument("--scftype", type=str, default="rhf")

    backend_group = parser.add_argument_group("backend refinement")
    backend_group.add_argument(
        "--software",
        type=str,
        choices=[
            "gaussian",
            "mopac",
            "obabel",
            "orca",
            "psi4",
            "turbomole",
            "xtb",
            "xtb_turbo",
            "mlatom_aiqm1",
            "aimnet_2",
            "aiqm1_mlatom",
            "xtb-aimnet2",
            "xtb-aiqm1",
        ],
    )
    backend_group.add_argument(
        "--geometry-optimizer",
        type=str,
        default=defualt_parameters.values["geometry_optimizer"],
        choices=["native", "geometric"],
    )
    backend_group.add_argument(
        "--opt-target",
        type=str,
        default=defualt_parameters.values["opt_target"],
        choices=["minimum", "ts"],
    )
    backend_group.add_argument("-basis", "--basis", type=str, default=defualt_parameters.values["basis"])
    backend_group.add_argument("-method", "--method", type=str, default=defualt_parameters.values["method"])
    backend_group.add_argument("--opt-threshold", type=str, default=defualt_parameters.values["opt_threshold"])
    backend_group.add_argument("--opt-cycles", type=int, default=defualt_parameters.values["opt_cycles"])
    backend_group.add_argument("--scf-threshold", type=str, default=defualt_parameters.values["scf_threshold"])
    backend_group.add_argument("--scf-cycles", type=int, default=defualt_parameters.values["scf_cycles"])
    backend_group.add_argument("-nprocs", "--nprocs", type=int, default=defualt_parameters.values["nprocs"])
    backend_group.add_argument("--custom-keywords", type=str)
    backend_group.add_argument("-model", "--model", type=str, default=defualt_parameters.values["model"])
    return parser.parse_args(argv)


def _backend_qc_params(args):
    """Return backend parameters when backend refinement was requested."""
    if args.software is None:
        return None
    return {
        "basis": args.basis,
        "method": args.method,
        "software": args.software,
        "geometry_optimizer": args.geometry_optimizer,
        "opt_target": args.opt_target,
        "opt_cycles": args.opt_cycles,
        "opt_threshold": args.opt_threshold,
        "scf_cycles": args.scf_cycles,
        "scf_threshold": args.scf_threshold,
        "nprocs": args.nprocs,
        "gamma": None,
        "custom_keywords": args.custom_keywords,
        "custom_keyword": args.custom_keywords,
        "model": args.model,
    }


def main(argv=None):
    """Run the RDKit conformer-search workflow."""
    args = argument_parse(argv)
    qc_params = _backend_qc_params(args)
    if qc_params is not None:
        from pyar.cli import _preflight_cli_requirements

        _preflight_cli_requirements(
            "conformer",
            qc_params["software"],
            qc_params["geometry_optimizer"],
        )

    try:
        result = conformer_search(
            args.input,
            input_format=args.input_format,
            num_conformers=args.num_conformers,
            top_n=args.top_n,
            backend_top_n=args.backend_top_n,
            num_seeds=args.num_seeds,
            diversity_fraction=args.diversity_fraction,
            rms_threshold=args.rms_threshold,
            use_random_coords=args.use_random_coords,
            force_field=args.force_field,
            seed=args.seed,
            num_threads=args.num_threads,
            max_iterations=args.max_iterations,
            charge=args.charge,
            multiplicity=args.multiplicity,
            scftype=args.scftype,
            qc_params=qc_params,
        )
    except (ConformerWorkflowError, FileNotFoundError, ImportError, ValueError) as exc:
        raise SystemExit(str(exc)) from exc

    print(f"Conformer workflow: {result.status}")
    print(f"Run directory: {result.run_directory}")
    print(f"Selected conformers: {len(result.selected_paths)}")
    return result


if __name__ == "__main__":
    main(sys.argv[1:])
