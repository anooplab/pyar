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
        else:
            result = collect_ts_optimizer_benchmark(args.output)
            print(json.dumps(result, indent=2))
    except (TSOptimizerBenchmarkError, OSError, ValueError) as exc:
        raise SystemExit(f"pyar-benchmark-ts: {exc}") from exc
    return None


if __name__ == "__main__":
    main()
