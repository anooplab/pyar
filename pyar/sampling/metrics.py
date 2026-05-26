"""Coverage metrics for sphere and rotation sampling."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from pyar.sampling.rotation import canonicalize_quaternion, halton_quaternions
from pyar.sampling.sphere import _normalize_points, _require_count, fibonacci_sphere


@dataclass(frozen=True)
class SphereCoverageMetrics:
    """Quality measures for directional sphere samples."""
    minimum_separation_degrees: float
    mean_nearest_neighbor_degrees: float
    covering_radius_degrees: float
    mean_covering_distance_degrees: float
    centroid_norm: float


@dataclass(frozen=True)
class RotationCoverageMetrics:
    """Quality measures for sampled monomer rotations."""
    minimum_separation_degrees: float
    mean_nearest_neighbor_degrees: float
    covering_radius_degrees: float
    mean_covering_distance_degrees: float


def sphere_coverage_metrics(
    points: np.ndarray,
    *,
    evaluation_points: int = 10000,
) -> SphereCoverageMetrics:
    """Measure angular spread and covering quality for sphere directions."""
    samples = _normalize_points(points)
    if len(samples) < 2:
        minimum_separation = 180.0
        mean_nearest_neighbor = 180.0
    else:
        similarities = np.clip(samples @ samples.T, -1.0, 1.0)
        np.fill_diagonal(similarities, -1.0)
        nearest_angles = np.degrees(np.arccos(np.max(similarities, axis=1)))
        minimum_separation = float(np.min(nearest_angles))
        mean_nearest_neighbor = float(np.mean(nearest_angles))

    reference = fibonacci_sphere(_require_count(evaluation_points))
    nearest_similarity = np.max(np.clip(reference @ samples.T, -1.0, 1.0), axis=1)
    covering_distances = np.degrees(np.arccos(nearest_similarity))
    return SphereCoverageMetrics(
        minimum_separation_degrees=minimum_separation,
        mean_nearest_neighbor_degrees=mean_nearest_neighbor,
        covering_radius_degrees=float(np.max(covering_distances)),
        mean_covering_distance_degrees=float(np.mean(covering_distances)),
        centroid_norm=float(np.linalg.norm(np.mean(samples, axis=0))),
    )


def _quaternion_distance_degrees(left: np.ndarray, right: np.ndarray) -> float:
    dot = float(abs(np.dot(left, right)))
    dot = float(np.clip(dot, -1.0, 1.0))
    return float(np.degrees(2.0 * np.arccos(dot)))


def quaternion_coverage_metrics(
    quaternions: np.ndarray,
    *,
    evaluation_points: int = 10000,
) -> RotationCoverageMetrics:
    """Measure angular spread and covering quality for quaternions."""
    samples = np.asarray(quaternions, dtype=float)
    if samples.ndim != 2 or samples.shape[1] != 4 or len(samples) == 0:
        raise ValueError("quaternions must have shape (N, 4) with N > 0")
    samples = np.array([canonicalize_quaternion(q) for q in samples], dtype=float)

    if len(samples) < 2:
        minimum_separation = 180.0
        mean_nearest_neighbor = 180.0
    else:
        similarities = np.clip(np.abs(samples @ samples.T), -1.0, 1.0)
        np.fill_diagonal(similarities, -1.0)
        nearest_angles = np.degrees(2.0 * np.arccos(np.max(similarities, axis=1)))
        minimum_separation = float(np.min(nearest_angles))
        mean_nearest_neighbor = float(np.mean(nearest_angles))

    reference = halton_quaternions(_require_count(evaluation_points))
    reference = np.array([canonicalize_quaternion(q) for q in reference], dtype=float)
    nearest_distances = []
    for ref in reference:
        nearest_distances.append(min(_quaternion_distance_degrees(ref, sample) for sample in samples))
    return RotationCoverageMetrics(
        minimum_separation_degrees=minimum_separation,
        mean_nearest_neighbor_degrees=mean_nearest_neighbor,
        covering_radius_degrees=float(np.max(nearest_distances)),
        mean_covering_distance_degrees=float(np.mean(nearest_distances)),
    )
