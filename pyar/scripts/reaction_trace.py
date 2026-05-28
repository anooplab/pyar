#!/usr/bin/env python3
"""Standalone reaction-trace analysis and plotting entrypoint."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from pyar.reaction_analysis import analyse_reaction_trace, plot_reaction_trace


def _build_parser():
    """Parse reaction-trace CLI arguments."""
    parser = argparse.ArgumentParser(
        prog="pyar-reaction-trace",
        description=(
            "Analyze a PyAR reaction trace, write the summary artifacts, and "
            "optionally generate PNG plots."
        ),
    )
    parser.add_argument(
        "path",
        nargs="?",
        default=".",
        help="Job directory or reaction_trace directory to analyze",
    )
    parser.add_argument(
        "--plot",
        action="store_true",
        help="Generate PNG plots from the trace after analysis",
    )
    parser.add_argument(
        "--plot-only",
        action="store_true",
        help="Generate PNG plots without rewriting summary artifacts",
    )
    parser.add_argument(
        "--plot-directory",
        type=str,
        default=None,
        help="Directory for generated plots; defaults to <job>/trace_plots",
    )
    return parser


def argument_parse(argv=None):
    """Parse reaction-trace CLI arguments."""
    return _build_parser().parse_args(argv)


def main(argv=None):
    """Run trace analysis and optionally emit plot artifacts."""
    args = argument_parse(argv)
    job_path = Path(args.path)
    plot_only = bool(getattr(args, "plot_only", False))
    do_plot = bool(getattr(args, "plot", False)) or plot_only
    plot_directory = getattr(args, "plot_directory", None)
    summary = None
    if not plot_only:
        summary = analyse_reaction_trace(job_path)
        if summary is None:
            raise SystemExit(f"No reaction trace records found in {job_path}")

    result = {}
    if summary is not None:
        result["analysis"] = summary
    if do_plot:
        plot_summary = plot_reaction_trace(job_path, plot_directory)
        if plot_summary is None:
            raise SystemExit(f"No reaction trace records found in {job_path}")
        result["plots"] = plot_summary

    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
