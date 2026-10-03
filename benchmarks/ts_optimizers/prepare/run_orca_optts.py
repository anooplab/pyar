"""Run ORCA OptTS on cases from a prepared TS benchmark manifest."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parents[3]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from pyar.benchmarks.orca_optts import run_orca_optts_case
from pyar.benchmarks.ts_optimizer import load_ts_optimizer_benchmark


def argument_parse(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", help="prepared TS benchmark JSON manifest")
    parser.add_argument("--output", required=True, help="ORCA run directory")
    parser.add_argument("--case", action="append", help="case ID; repeat to select cases")
    parser.add_argument("--orca-executable", help="ORCA executable path (default: PATH lookup)")
    parser.add_argument("--xtb-executable", help="xTB executable path (default: PATH lookup)")
    parser.add_argument(
        "--match-ts-fmax", action="store_true",
        help="set ORCA gradient thresholds from the manifest ts_fmax value",
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = argument_parse(argv)
    spec = load_ts_optimizer_benchmark(args.manifest)
    selected = set(args.case or (case.id for case in spec.cases))
    unknown = selected - {case.id for case in spec.cases}
    if unknown:
        raise SystemExit("unknown case ID(s): " + ", ".join(sorted(unknown)))
    output = Path(args.output).expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    for case in spec.cases:
        if case.id not in selected:
            continue
        try:
            result = run_orca_optts_case(
                spec, case_id=case.id, output=output,
                orca_executable=args.orca_executable,
                xtb_executable=args.xtb_executable,
                match_ts_fmax=args.match_ts_fmax,
            )
            print(json.dumps({
                "case_id": case.id,
                "outcome": result["outcome"],
                "steps": result["ts_optimizer_steps"],
                "gradient_events": result["ts_orca_gradient_events"],
                "wall_seconds": result["ts_wall_seconds"],
            }), flush=True)
        except Exception as exc:
            print(json.dumps({
                "case_id": case.id,
                "error": f"{type(exc).__name__}: {exc}",
            }), flush=True)
    return None


if __name__ == "__main__":
    main()
