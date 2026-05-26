"""Sphere sampling helpers for trial placement."""

from __future__ import annotations

import math

import numpy as np
from scipy.stats import qmc


_GOLDEN_ANGLE = math.pi * (3.0 - math.sqrt(5.0))


def _require_count(number_of_points: int) -> int:
    value = int(number_of_points)
    if value < 1:
        raise ValueError("number_of_points must be positive")
    return value


def _normalize_points(points: np.ndarray) -> np.ndarray:
    values = np.asarray(points, dtype=float)
    if values.ndim != 2 or values.shape[1] != 3 or len(values) == 0:
        raise ValueError("points must have shape (N, 3) with N > 0")
    norms = np.linalg.norm(values, axis=1)
    if np.any(norms == 0.0):
        raise ValueError("sphere points cannot contain zero vectors")
    return values / norms[:, None]


def fibonacci_sphere(number_of_points: int, offset: float = 0.5, azimuth_offset: float = 0.0) -> np.ndarray:
    """Generate deterministic approximately uniform directions on a sphere."""
    count = _require_count(number_of_points)
    indices = np.arange(count, dtype=float)
    z = 1.0 - (2.0 * (indices + float(offset)) / count)
    radius = np.sqrt(np.maximum(0.0, 1.0 - z * z))
    azimuth = indices * _GOLDEN_ANGLE + float(azimuth_offset)
    return np.column_stack((radius * np.cos(azimuth), radius * np.sin(azimuth), z))


def latin_hypercube_sphere(number_of_points: int, seed: int | None = None) -> np.ndarray:
    """Generate sphere directions from a Latin-hypercube sample."""
    count = _require_count(number_of_points)
    samples = qmc.LatinHypercube(d=2, seed=seed).random(n=count)
    z = 1.0 - 2.0 * samples[:, 0]
    azimuth = 2.0 * np.pi * samples[:, 1]
    radius = np.sqrt(np.maximum(0.0, 1.0 - z * z))
    return np.column_stack((radius * np.cos(azimuth), radius * np.sin(azimuth), z))


def random_sphere(number_of_points: int, seed: int | None = None) -> np.ndarray:
    """Generate uniformly distributed random sphere directions."""
    count = _require_count(number_of_points)
    random = np.random.default_rng(seed)
    values = random.normal(size=(count, 3))
    return _normalize_points(values)


def maximin_select(points: np.ndarray, number_of_points: int) -> np.ndarray:
    """Select a spread-out subset of candidate directions greedily."""
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
        distance_to_new = np.sum((candidates - candidates[next_index]) ** 2, axis=1)
        min_squared_distances = np.minimum(min_squared_distances, distance_to_new)
        min_squared_distances[selected_indices] = -np.inf
    return candidates[selected_indices]


def generate_directions(
    method: str,
    number_of_points: int,
    *,
    seed: int | None = None,
    sequence_offset: int = 0,
    oversample_factor: int = 8,
) -> np.ndarray:
    """Generate trial-placement directions using a named sampling method."""
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
