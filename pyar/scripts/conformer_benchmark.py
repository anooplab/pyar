#!/usr/bin/env python3
"""Command-line entrypoint for conformer benchmark diagnosis."""

from __future__ import annotations

import argparse
import sys

from pyar.data import defualt_parameters
from pyar.benchmarks.conformer import ConformerBenchmarkError, run_conformer_benchmark


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


def argument_parse(argv=None):
    """Parse conformer benchmark CLI arguments."""
    parser = argparse.ArgumentParser(
        prog="pyar-conformer-benchmark",
        description="Benchmark pyar-conformer and classify missed-reference causes.",
    )
    parser.add_argument("benchmark", help="JSON conformer benchmark specification")
    parser.add_argument("--num-conformers", type=int, default=100)
    parser.add_argument("--num-seeds", type=int, default=3)
    parser.add_argument("--top-n", type=int, default=10)
    parser.add_argument("--backend-top-n", type=int)
    parser.add_argument("--rms-hit-threshold", type=float, default=0.50)
    parser.add_argument(
        "--energy-window",
        type=float,
        help=(
            "Classify reference-like generated conformers outside this RDKit "
            "native force-field energy window as wrong_ranking."
        ),
    )
    parser.add_argument("--output", default="conformer_benchmark_results")

    conformer_group = parser.add_argument_group("conformer search")
    conformer_group.add_argument("--diversity-fraction", type=float, default=0.2)
    conformer_group.add_argument(
        "--compactness-fraction",
        type=float,
        default=0.2,
        help="Protected contact-rich folded-basin quota; matched by an open-basin quota for diversity.",
    )
    conformer_group.add_argument("--rms-threshold", type=float, default=0.25)
    conformer_group.add_argument(
        "--use-random-coords",
        action=argparse.BooleanOptionalAction,
        default=True,
    )
    conformer_group.add_argument(
        "--torsion-kicks",
        action=argparse.BooleanOptionalAction,
        default=True,
    )
    conformer_group.add_argument("--torsion-rounds", type=int, default=2)
    conformer_group.add_argument("--torsion-mode", choices=["evolve", "mc", "grid", "random"], default="evolve")
    conformer_group.add_argument("--torsion-kicks-per-conformer", type=int, default=8)
    conformer_group.add_argument("--torsion-max-bonds", type=int, default=3)
    conformer_group.add_argument("--torsion-dedup-rms", type=float, default=0.5)
    conformer_group.add_argument("--force-field", choices=["auto", "mmff", "uff"], default="auto")
    conformer_group.add_argument("--seed", type=int, default=1)
    conformer_group.add_argument("--num-threads", type=int, default=0)
    conformer_group.add_argument("--max-iterations", type=int, default=200)

    backend_group = parser.add_argument_group("backend refinement")
    backend_group.add_argument(
        "--software",
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
        default=defualt_parameters.values["geometry_optimizer"],
        choices=["native", "geometric"],
    )
    backend_group.add_argument(
        "--opt-target",
        default=defualt_parameters.values["opt_target"],
        choices=["minimum", "ts"],
    )
    backend_group.add_argument("-basis", "--basis", default=defualt_parameters.values["basis"])
    backend_group.add_argument("-method", "--method", default=defualt_parameters.values["method"])
    backend_group.add_argument("--opt-threshold", default=defualt_parameters.values["opt_threshold"])
    backend_group.add_argument("--opt-cycles", type=int, default=defualt_parameters.values["opt_cycles"])
    backend_group.add_argument("--scf-threshold", default=defualt_parameters.values["scf_threshold"])
    backend_group.add_argument("--scf-cycles", type=int, default=defualt_parameters.values["scf_cycles"])
    backend_group.add_argument("-nprocs", "--nprocs", type=int, default=defualt_parameters.values["nprocs"])
    backend_group.add_argument("--custom-keywords")
    backend_group.add_argument("-model", "--model", default=defualt_parameters.values["model"])
    return parser.parse_args(argv)


def _conformer_options(args):
    return {
        "diversity_fraction": args.diversity_fraction,
        "compactness_fraction": args.compactness_fraction,
        "rms_threshold": args.rms_threshold,
        "use_random_coords": args.use_random_coords,
        "torsion_kicks": args.torsion_kicks,
        "torsion_rounds": args.torsion_rounds,
        "torsion_mode": args.torsion_mode,
        "torsion_kicks_per_conformer": args.torsion_kicks_per_conformer,
        "torsion_max_bonds": args.torsion_max_bonds,
        "torsion_dedup_rms": args.torsion_dedup_rms,
        "force_field": args.force_field,
        "seed": args.seed,
        "num_threads": args.num_threads,
        "max_iterations": args.max_iterations,
    }


def main(argv=None):
    """Run conformer benchmark diagnosis from the command line."""
    args = argument_parse(argv)
    qc_params = _backend_qc_params(args)
    if qc_params is not None:
        from pyar.cli import _preflight_cli_requirements

        _preflight_cli_requirements(
            "conformer-benchmark",
            qc_params["software"],
            qc_params["geometry_optimizer"],
        )
    try:
        result = run_conformer_benchmark(
            args.benchmark,
            output=args.output,
            num_conformers=args.num_conformers,
            num_seeds=args.num_seeds,
            top_n=args.top_n,
            backend_top_n=args.backend_top_n,
            software=args.software,
            rms_hit_threshold=args.rms_hit_threshold,
            energy_window=args.energy_window,
            conformer_options=_conformer_options(args),
            qc_params=qc_params,
        )
    except (ConformerBenchmarkError, FileNotFoundError, ValueError) as exc:
        raise SystemExit(str(exc)) from exc

    print(f"Conformer benchmark: {result['name']}")
    print(f"Cases: {len(result['cases'])}")
    print(f"Output: {args.output}")
    return None


if __name__ == "__main__":
    main(sys.argv[1:])
