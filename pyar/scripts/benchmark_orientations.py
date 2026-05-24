#!/usr/bin/env python3
"""Benchmark direction samplers used for trial geometry placement."""

from __future__ import annotations

import argparse
import time

import numpy as np

from pyar.orientation_sampling import generate_directions, sphere_coverage_metrics


DEFAULT_METHODS = [
    "random",
    "lhs",
    "lhs_maximin",
    "fibonacci",
    "fibonacci_maximin",
]


def _benchmark(method, count, repeat, evaluation_points, oversample_factor):
    metrics = []
    runtimes = []
    for index in range(repeat):
        seed = index if method in {
            "random",
            "lhs",
            "lhs_maximin",
        } else None
        start = time.perf_counter()
        points = generate_directions(
            method,
            count,
            seed=seed,
            oversample_factor=oversample_factor,
        )
        metrics.append(
            sphere_coverage_metrics(points, evaluation_points=evaluation_points)
        )
        runtimes.append(time.perf_counter() - start)
    values = np.array([
        [
            metric.minimum_separation_degrees,
            metric.mean_nearest_neighbor_degrees,
            metric.covering_radius_degrees,
            metric.mean_covering_distance_degrees,
            metric.centroid_norm,
        ]
        for metric in metrics
    ])
    return np.mean(values, axis=0), np.std(values, axis=0), np.mean(runtimes)


def main():
    """Run reproducible unit-sphere coverage comparisons."""
    parser = argparse.ArgumentParser(
        description="Benchmark direction samplers for trial geometry generation."
    )
    parser.add_argument(
        "-N",
        "--number-of-points",
        nargs="+",
        type=int,
        default=[8, 12, 20, 28],
        help="Numbers of placement directions to benchmark.",
    )
    parser.add_argument(
        "-m",
        "--methods",
        nargs="+",
        default=DEFAULT_METHODS,
        choices=DEFAULT_METHODS,
        help="Sampling methods to compare.",
    )
    parser.add_argument(
        "-r",
        "--repeat",
        type=int,
        default=10,
        help="Repeats for randomized methods; deterministic methods repeat identically.",
    )
    parser.add_argument(
        "--evaluation-points",
        type=int,
        default=10000,
        help="Dense reference directions used to approximate covering radius.",
    )
    parser.add_argument(
        "--oversample-factor",
        type=int,
        default=8,
        help="Candidate-pool multiplier for maximin methods.",
    )
    args = parser.parse_args()

    print(
        "N\tmethod\tmin_sep_deg\tmean_nn_deg\tcover_radius_deg"
        "\tmean_cover_deg\tcentroid_norm\truntime_s"
    )
    for count in args.number_of_points:
        for method in args.methods:
            mean, std, runtime = _benchmark(
                method,
                count,
                args.repeat,
                args.evaluation_points,
                args.oversample_factor,
            )
            variation = "" if method.startswith("fibonacci") else (
                f" (+/- {std[0]:.3f}, {std[2]:.3f})"
            )
            print(
                f"{count}\t{method}\t{mean[0]:.3f}\t{mean[1]:.3f}\t"
                f"{mean[2]:.3f}\t{mean[3]:.3f}\t{mean[4]:.6f}\t"
                f"{runtime:.6f}{variation}"
            )


if __name__ == "__main__":
    main()
