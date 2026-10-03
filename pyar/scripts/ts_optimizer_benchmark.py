"""CLI for paired transition-state optimizer benchmark runs and collection."""

from __future__ import annotations

import argparse
import json

from pyar.benchmarks.ts_optimizer import (
    OPTIMIZERS,
    TSOptimizerBenchmarkError,
    collect_ts_optimizer_benchmark,
    run_ts_optimizer_case,
)
from pyar.benchmarks.rgd1_tsopt import TIERS, prepare_rgd1_tsopt_manifest


def argument_parse(argv=None):
    parser = argparse.ArgumentParser(
        prog="pyar-benchmark-ts",
        description="Run or collect paired PyAR transition-state optimizer benchmarks.",
    )
    commands = parser.add_subparsers(dest="command", required=True)
    run = commands.add_parser("run", help="run one independently scheduled case/optimizer job")
    run.add_argument("benchmark", help="JSON benchmark manifest")
    run.add_argument("--case", required=True, help="case ID from the manifest")
    run.add_argument("--optimizer", required=True, choices=OPTIMIZERS)
    run.add_argument("--output", required=True, help="shared benchmark run directory")
    collect = commands.add_parser("collect", help="collect complete and incomplete paired jobs")
    collect.add_argument("output", help="benchmark run directory containing run_manifest.json")
    prepare = commands.add_parser(
        "prepare-rgd1",
        help="create a small manifest from an extracted RGD1-TSopt-GFN2 dataset",
    )
    prepare.add_argument("--dataset-dir", required=True, help="extracted dataset directory")
    prepare.add_argument("--output", required=True, help="new output JSON manifest path")
    prepare.add_argument(
        "--reactions", type=int, default=15,
        help="number of distinct reactions to sample (default: 15)",
    )
    prepare.add_argument(
        "--seed", type=int, default=20261002,
        help="deterministic atom-count-stratified sampling seed",
    )
    prepare.add_argument("--tiers", nargs="+", choices=TIERS, default=list(TIERS),
                         help="distortion tiers to include (default: easy med hard)")
    prepare.add_argument(
        "--ts-fmax", type=float, default=0.02,
        help="shared optimizer force threshold in eV/angstrom",
    )
    prepare.add_argument(
        "--ts-max-cycles", type=int, default=200,
        help="shared optimizer step ceiling",
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = argument_parse(argv)
    try:
        if args.command == "run":
            result = run_ts_optimizer_case(
                args.benchmark, case_id=args.case, optimizer=args.optimizer, output=args.output,
            )
            print(json.dumps({
                "case_id": result["case_id"], "optimizer": result["optimizer"],
                "status": result["status"], "outcome": result["outcome"],
                "result": f"{result['run_directory']}/benchmark_result.json",
            }, indent=2))
        elif args.command == "collect":
            result = collect_ts_optimizer_benchmark(args.output)
            print(json.dumps(result, indent=2))
        else:
            result = prepare_rgd1_tsopt_manifest(
                args.dataset_dir, args.output, reactions=args.reactions, seed=args.seed,
                tiers=args.tiers, ts_fmax=args.ts_fmax,
                ts_max_cycles=args.ts_max_cycles,
            )
            print(result)
    except (TSOptimizerBenchmarkError, OSError, ValueError) as exc:
        raise SystemExit(f"pyar-benchmark-ts: {exc}") from exc
    return None


if __name__ == "__main__":
    main()
