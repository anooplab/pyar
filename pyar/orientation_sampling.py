"""Sampling and quality metrics for trial placement and fragment rotation.

Spatial approach directions and rotational orientations are generated and
measured separately. Rotation must not enter the sphere-coverage metrics.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
from scipy.stats import qmc


_GOLDEN_ANGLE = math.pi * (3.0 - math.sqrt(5.0))


@dataclass(frozen=True)
class SphereCoverageMetrics:
    """Quality measurements for a unit-sphere point set.

    Attributes
    ----------
    minimum_separation_degrees
        Smallest pairwise angular distance. Higher is better.
    mean_nearest_neighbor_degrees
        Mean distance to each point's closest neighbour. Higher indicates
        better separation, when considered with covering radius.
    covering_radius_degrees
        Largest distance from a dense evaluation grid to its closest sampled
        point. Lower is better.
    mean_covering_distance_degrees
        Mean distance from the evaluation grid to its closest sampled point.
        Lower is better.
    centroid_norm
        Norm of the point-set centroid. Lower indicates less directional bias.
    """

    minimum_separation_degrees: float
    mean_nearest_neighbor_degrees: float
    covering_radius_degrees: float
    mean_covering_distance_degrees: float
    centroid_norm: float


def _require_count(number_of_points: int) -> int:
    """Validate and normalize a requested point count."""
    value = int(number_of_points)
    if value < 1:
        raise ValueError("number_of_points must be positive")
    return value


def _normalize_points(points: np.ndarray) -> np.ndarray:
    """Return validated rows normalized to unit length."""
    values = np.asarray(points, dtype=float)
    if values.ndim != 2 or values.shape[1] != 3 or len(values) == 0:
        raise ValueError("points must have shape (N, 3) with N > 0")
    norms = np.linalg.norm(values, axis=1)
    if np.any(norms == 0.0):
        raise ValueError("sphere points cannot contain zero vectors")
    return values / norms[:, None]


def fibonacci_sphere(
    number_of_points: int,
    offset: float = 0.5,
    azimuth_offset: float = 0.0,
) -> np.ndarray:
    """Generate deterministic equal-area spherical Fibonacci points.

    Parameters
    ----------
    number_of_points
        Number of unit vectors to generate.
    offset
        Half-step offset along the vertical equal-area coordinate. The
        default avoids placing points exactly at the poles.
    azimuth_offset
        Global azimuthal phase in radians. The default preserves the canonical
        deterministic design.
    """
    count = _require_count(number_of_points)
    indices = np.arange(count, dtype=float)
    z = 1.0 - (2.0 * (indices + float(offset)) / count)
    radius = np.sqrt(np.maximum(0.0, 1.0 - z * z))
    azimuth = indices * _GOLDEN_ANGLE + float(azimuth_offset)
    return np.column_stack((radius * np.cos(azimuth),
                            radius * np.sin(azimuth),
                            z))


def latin_hypercube_sphere(number_of_points: int, seed: int | None = None) -> np.ndarray:
    """Generate equal-area sphere points from a Latin hypercube design."""
    count = _require_count(number_of_points)
    samples = qmc.LatinHypercube(d=2, seed=seed).random(n=count)
    z = 1.0 - 2.0 * samples[:, 0]
    azimuth = 2.0 * np.pi * samples[:, 1]
    radius = np.sqrt(np.maximum(0.0, 1.0 - z * z))
    return np.column_stack((radius * np.cos(azimuth),
                            radius * np.sin(azimuth),
                            z))


def random_sphere(number_of_points: int, seed: int | None = None) -> np.ndarray:
    """Generate uniformly random sphere points for a baseline comparison."""
    count = _require_count(number_of_points)
    random = np.random.default_rng(seed)
    values = random.normal(size=(count, 3))
    return _normalize_points(values)


def _halton_value(index: int, base: int) -> float:
    """Return one radical-inverse Halton value."""
    value = 0.0
    factor = 1.0
    current = int(index)
    while current > 0:
        factor /= base
        value += factor * (current % base)
        current //= base
    return value


def _halton_sequence(number_of_points: int, bases=(2, 3, 5), offset: int = 0) -> np.ndarray:
    """Return a Halton sequence in the unit cube."""
    count = _require_count(number_of_points)
    sequence = np.empty((count, len(bases)), dtype=float)
    for row in range(count):
        index = row + 1 + int(offset)
        for column, base in enumerate(bases):
            sequence[row, column] = _halton_value(index, int(base))
    return sequence


def random_quaternions(number_of_points: int, seed: int | None = None) -> np.ndarray:
    """Generate uniform random unit quaternions with scalar part last."""
    count = _require_count(number_of_points)
    random = np.random.default_rng(seed)
    values = random.normal(size=(count, 4))
    norms = np.linalg.norm(values, axis=1)
    if np.any(norms == 0.0):
        raise ValueError("random quaternion sampling produced a zero vector")
    return values / norms[:, None]


def halton_quaternions(number_of_points: int, seed: int | None = None) -> np.ndarray:
    """Generate deterministic low-discrepancy unit quaternions.

    The implementation uses a Halton sequence for the three uniform inputs in
    Shoemake's quaternion sampling map, which gives a simple and clean
    rotational default for trial geometries.
    """
    count = _require_count(number_of_points)
    offset = 0 if seed is None else int(seed)
    cube = _halton_sequence(count, bases=(2, 3, 5), offset=offset)
    u1 = cube[:, 0]
    u2 = cube[:, 1]
    u3 = cube[:, 2]
    sqrt_one_minus_u1 = np.sqrt(1.0 - u1)
    sqrt_u1 = np.sqrt(u1)
    qx = sqrt_one_minus_u1 * np.sin(2.0 * np.pi * u2)
    qy = sqrt_one_minus_u1 * np.cos(2.0 * np.pi * u2)
    qz = sqrt_u1 * np.sin(2.0 * np.pi * u3)
    qw = sqrt_u1 * np.cos(2.0 * np.pi * u3)
    return np.column_stack((qx, qy, qz, qw))


def canonicalize_quaternion(quaternion: np.ndarray) -> np.ndarray:
    """Return a quaternion with a unique sign convention."""
    values = np.asarray(quaternion, dtype=float).reshape(4)
    if values[3] < 0.0:
        values = -values
    return values / np.linalg.norm(values)


def quaternions_to_euler_zxz(quaternions: np.ndarray) -> np.ndarray:
    """Convert quaternions to ZXZ Euler angles compatible with ``rotate_3d``."""
    from scipy.spatial.transform import Rotation

    values = np.asarray(quaternions, dtype=float)
    if values.ndim == 1:
        values = values.reshape(1, 4)
    rotation = Rotation.from_quat(values)
    return rotation.as_euler("zxz", degrees=False)


def maximin_select(points: np.ndarray, number_of_points: int) -> np.ndarray:
    """Choose directions greedily to maximize nearest angular separation.

    The first selected candidate is the one closest to the north pole. After
    that, each selection maximizes its minimum chord distance to the already
    selected set. This is deterministic and appropriate for reducing an
    oversized candidate set.
    """
    candidates = _normalize_points(points)
    count = _require_count(number_of_points)
    if count > len(candidates):
        raise ValueError("cannot select more points than candidate directions")
    if count == len(candidates):
        return candidates.copy()

    selected_indices = [int(np.argmax(candidates[:, 2]))]
    min_squared_distances = np.sum(
        (candidates - candidates[selected_indices[0]]) ** 2,
        axis=1,
    )
    min_squared_distances[selected_indices[0]] = -np.inf
    while len(selected_indices) < count:
        next_index = int(np.argmax(min_squared_distances))
        selected_indices.append(next_index)
        distance_to_new = np.sum(
            (candidates - candidates[next_index]) ** 2,
            axis=1,
        )
        min_squared_distances = np.minimum(min_squared_distances, distance_to_new)
        min_squared_distances[selected_indices] = -np.inf
    return candidates[selected_indices]


def _sphere_to_trial_vectors(directions: np.ndarray, quaternions: np.ndarray | None = None, use_angles: bool = True) -> np.ndarray:
    """Combine spatial directions with rotation parameters."""
    normalized_directions = _normalize_points(directions)
    if not use_angles:
        angles = np.zeros((len(normalized_directions), 3), dtype=float)
    else:
        if quaternions is None:
            raise ValueError("quaternions are required when use_angles is True")
        canonical_quaternions = np.array([canonicalize_quaternion(q) for q in quaternions], dtype=float)
        angles = quaternions_to_euler_zxz(canonical_quaternions)
    return np.column_stack((normalized_directions, angles))


def generate_trial_vectors(
    number_of_points: int,
    *,
    direction_method: str = "fibonacci",
    rotation_method: str = "halton",
    seed: int | None = None,
    sequence_offset: int = 0,
    oversample_factor: int = 8,
    use_angles: bool = True,
) -> np.ndarray:
    """Generate direction-plus-rotation trial vectors.

    The default pair is spherical Fibonacci directions plus Halton-based
    quaternions. This keeps the sphere coverage clean and the monomer
    orientation sampling low-discrepancy.
    """
    count = _require_count(number_of_points)
    directions = generate_directions(
        direction_method,
        count,
        seed=seed,
        sequence_offset=sequence_offset,
        oversample_factor=oversample_factor,
    )
    if use_angles:
        rotation_offset = None
        if sequence_offset or seed is not None:
            rotation_offset = int(sequence_offset) * count
            if seed is not None:
                rotation_offset += int(seed)
        if rotation_method == "halton":
            quaternions = halton_quaternions(count, seed=rotation_offset)
        elif rotation_method == "random":
            quaternions = random_quaternions(count, seed=rotation_offset)
        else:
            raise ValueError(f"Unknown rotation method: {rotation_method}")
    else:
        quaternions = None

    return _sphere_to_trial_vectors(directions, quaternions, use_angles=use_angles)


def generate_directions(
    method: str,
    number_of_points: int,
    *,
    seed: int | None = None,
    sequence_offset: int = 0,
    oversample_factor: int = 8,
) -> np.ndarray:
    """Generate unit-sphere directions by a named benchmark method."""
    count = _require_count(number_of_points)
    if method == "fibonacci":
        if sequence_offset:
            phase = (int(sequence_offset) * ((math.sqrt(5.0) - 1.0) / 2.0)) % 1.0
            return fibonacci_sphere(
                count,
                offset=0.25 + 0.5 * phase,
                azimuth_offset=2.0 * math.pi * phase,
            )
        return fibonacci_sphere(count)
    if method == "lhs":
        return latin_hypercube_sphere(count, seed=seed)
    if method == "random":
        return random_sphere(count, seed=seed)
    if method == "lhs_maximin":
        candidates = latin_hypercube_sphere(count * oversample_factor, seed=seed)
        return maximin_select(candidates, count)
    if method == "fibonacci_maximin":
        candidates = fibonacci_sphere(count * oversample_factor)
        return maximin_select(candidates, count)
    raise ValueError(f"Unknown sphere sampling method: {method}")


def sphere_coverage_metrics(
    points: np.ndarray,
    *,
    evaluation_points: int = 10000,
) -> SphereCoverageMetrics:
    """Measure angular separation and approximate full-sphere coverage.

    Coverage is approximated against a dense deterministic Fibonacci grid.
    This gives reproducible relative comparisons between sampling methods.
    """
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


@dataclass(frozen=True)
class RotationCoverageMetrics:
    """Quality measurements for a unit-quaternion sample."""

    minimum_separation_degrees: float
    mean_nearest_neighbor_degrees: float
    covering_radius_degrees: float
    mean_covering_distance_degrees: float


def _quaternion_distance_degrees(left: np.ndarray, right: np.ndarray) -> float:
    dot = float(abs(np.dot(left, right)))
    dot = float(np.clip(dot, -1.0, 1.0))
    return float(np.degrees(2.0 * np.arccos(dot)))


def quaternion_coverage_metrics(
    quaternions: np.ndarray,
    *,
    evaluation_points: int = 10000,
) -> RotationCoverageMetrics:
    """Measure separation and approximate coverage on ``SO(3)``."""
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
