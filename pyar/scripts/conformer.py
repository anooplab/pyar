#!/usr/bin/env python3
"""Command-line entrypoint for RDKit conformational search."""

from __future__ import annotations

import argparse
import logging
import sys

from pyar.data import defualt_parameters
from pyar.workflows.conformer import ConformerWorkflowError, conformer_search

logger = logging.getLogger("pyar-conformer")


def build_parser(prog=None, *, modern=False):
    """Parse conformer-search CLI arguments."""
    parser = argparse.ArgumentParser(
        prog=prog or "pyar-cli conformer",
        description="Generate and optionally refine RDKit conformers.",
    )
    parser.add_argument("input", help="SMILES string, SDF/MOL file, or XYZ file")
    parser.add_argument(
        "--input-format",
        choices=["auto", "smiles", "sdf", "mol", "xyz"],
        default="auto",
        help="Input format; auto treats existing .xyz/.sdf/.mol paths as files and other input as SMILES.",
    )
    parser.add_argument("--num-conformers", type=int, default=150)
    parser.add_argument("--top-n", type=int, default=10)
    parser.add_argument("--backend-top-n", type=int)
    advanced = parser.add_argument_group("Advanced generation controls") if modern else parser
    advanced.add_argument("--num-seeds", type=int, default=5)
    advanced.add_argument("--diversity-fraction", type=float, default=0.2)
    advanced.add_argument(
        "--compactness-fraction",
        type=float,
        default=0.2,
        help="Protected contact-rich folded-basin quota; matched by an open-basin quota for diversity.",
    )
    advanced.add_argument(
        "--rms-threshold",
        "--prune-rms-threshold",
        dest="rms_threshold",
        type=float,
        default=0.25,
        help="RDKit greedy prune RMS threshold; lower values keep more embedded conformers.",
    )
    advanced.add_argument(
        "--use-random-coords",
        dest="use_random_coords",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Start RDKit embedding from random coordinates instead of distance geometry eigenvectors.",
    )
    advanced.add_argument(
        "--torsion-kicks",
        dest="torsion_kicks",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Generate local torsion-perturbed conformers before backend refinement.",
    )
    advanced.add_argument(
        "--torsion-mode",
        choices=["random"],
        default="random",
        help="Use the stratified random torsion-kick sampler.",
    )
    advanced.add_argument("--torsion-rounds", type=int, default=2)
    advanced.add_argument("--torsion-kicks-per-conformer", type=int, default=6)
    advanced.add_argument("--torsion-max-bonds", type=int, default=3)
    advanced.add_argument("--torsion-dedup-rms", type=float, default=0.5)
    advanced.add_argument(
        "--dedup-atom-mode", choices=["heavy", "all"], default="heavy",
        help="Atoms scored by shared graph RMSD; full connectivity is always checked.",
    )
    parser.add_argument("--force-field", choices=["auto", "mmff", "uff"], default="auto")
    advanced.add_argument("--seed", type=int, default=1)
    advanced.add_argument("--num-threads", type=int, default=0)
    advanced.add_argument("--max-iterations", type=int, default=200)

    molecule_group = parser.add_argument_group("molecule")
    molecule_group.add_argument("-c", "--charge", type=int)
    molecule_group.add_argument("-m", "--multiplicity", type=int, default=None if modern else 1)
    molecule_group.add_argument("--scftype", type=str, default=None if modern else "rhf")

    if modern:
        from pyar.modern_workflow import add_backend_arguments
        add_backend_arguments(parser)
        parser.add_argument('--check', action='store_true', help='Validate without embedding or creating a run directory')
        return parser
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
    return parser


def argument_parse(argv=None):
    return build_parser().parse_args(argv)


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


def workflow_options(args, qc_params):
    """Shared mapping for legacy and modern interfaces to the canonical engine."""
    names = (
        'input_format', 'num_conformers', 'top_n', 'backend_top_n', 'num_seeds',
        'diversity_fraction', 'compactness_fraction', 'rms_threshold',
        'use_random_coords', 'torsion_kicks', 'torsion_mode', 'torsion_rounds',
        'torsion_kicks_per_conformer', 'torsion_max_bonds', 'torsion_dedup_rms',
        'dedup_atom_mode', 'force_field', 'seed', 'num_threads', 'max_iterations',
        'charge', 'multiplicity', 'scftype',
    )
    return {**{name: getattr(args, name) for name in names}, "qc_params": qc_params}


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
        result = conformer_search(args.input, **workflow_options(args, qc_params))
    except (ConformerWorkflowError, FileNotFoundError, ImportError, ValueError) as exc:
        raise SystemExit(str(exc)) from exc

    print(f"Conformer workflow: {result.status}")
    print(f"Run directory: {result.run_directory}")
    print(f"Selected conformers: {len(result.selected_paths)}")
    return None


if __name__ == "__main__":
    main(sys.argv[1:])
