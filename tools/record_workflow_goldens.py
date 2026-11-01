"""Record or check golden evidence using only the approved original checkout.

Run this developer tool with --baseline-root pointing to a clean checkout of
9f1e567. It starts a subprocess importing that checkout, while using the current
characterization harness. Never regenerate fixtures during a pure refactor.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

BASELINE = "9f1e567ce22b5d382d1facc9ff1da624de5d55a9"
CASES = (
    "conformer",
    "aggregate",
    "grow",
    "microsolvation",
    "microsolvation_resampling",
    "reaction",
    "scan_inter",
    "scan_intra",
)


def worker(args):
    import pytest

    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tests"))
    import pyar
    from workflow_golden_support import run_scenario

    assert Path(pyar.__file__).resolve().is_relative_to(args.baseline_root), (
        pyar.__file__
    )
    for iteration in range(args.repeat):
        for case in CASES:
            with tempfile.TemporaryDirectory(prefix="pyar-golden-") as temporary:
                with pytest.MonkeyPatch.context() as patch:
                    evidence = run_scenario(case, Path(temporary) / "run", patch)
            encoded = json.dumps(evidence, indent=2, sort_keys=True) + "\n"
            target = args.output / case / "evidence.json"
            if args.check or iteration:
                assert target.read_text() == encoded, (
                    f"{case}: original evidence changed on repeat {iteration + 1}"
                )
            else:
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_text(encoded)
            print(f"{case}: original evidence verified ({iteration + 1})", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-root", type=Path, required=True)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "tests/fixtures/golden",
    )
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--repeat", type=int, default=5)
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    args.baseline_root, args.output = (
        args.baseline_root.resolve(),
        args.output.resolve(),
    )
    if args.repeat < 1:
        parser.error("--repeat must be positive")
    revision = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=args.baseline_root, text=True
    ).strip()
    dirty = subprocess.check_output(
        ["git", "status", "--porcelain", "--", "pyar"],
        cwd=args.baseline_root,
        text=True,
    )
    if revision != BASELINE or dirty:
        parser.error("Baseline must be the clean approved 9f1e567 checkout")
    if args.worker:
        worker(args)
    else:
        environment = dict(os.environ, PYTHONPATH=str(args.baseline_root))
        subprocess.run(
            [sys.executable, str(Path(__file__).resolve()), *sys.argv[1:], "--worker"],
            env=environment,
            cwd=args.baseline_root,
            check=True,
        )


if __name__ == "__main__":
    main()
