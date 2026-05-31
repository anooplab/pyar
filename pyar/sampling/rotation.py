"""Rotation and orientation sampling helpers for trial placement."""

from __future__ import annotations

import numpy as np

from pyar.sampling.sphere import _require_count


def _halton_value(index: int, base: int) -> float:
    value = 0.0
    factor = 1.0
    current = int(index)
    while current > 0:
        factor /= base
        value += factor * (current % base)
        current //= base
    return value


def _halton_sequence(number_of_points: int, bases=(2, 3, 5), offset: int = 0) -> np.ndarray:
    count = _require_count(number_of_points)
    sequence = np.empty((count, len(bases)), dtype=float)
    for row in range(count):
        index = row + 1 + int(offset)
        for column, base in enumerate(bases):
            sequence[row, column] = _halton_value(index, int(base))
    return sequence


def random_quaternions(number_of_points: int, seed: int | None = None) -> np.ndarray:
    """Generate random unit quaternions for monomer rotation sampling."""
    count = _require_count(number_of_points)
    random = np.random.default_rng(seed)
    values = random.normal(size=(count, 4))
    norms = np.linalg.norm(values, axis=1)
    if np.any(norms == 0.0):
        raise ValueError("random quaternion sampling produced a zero vector")
    return values / norms[:, None]


def halton_quaternions(number_of_points: int, seed: int | None = None) -> np.ndarray:
    """Generate deterministic low-discrepancy unit quaternions.

    ``seed`` advances the Halton sequence deterministically so repeated
    calls with the same offset return the same rotation set.
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
    """Normalize a quaternion and choose a consistent sign representative."""
    values = np.asarray(quaternion, dtype=float).reshape(4)
    if values[3] < 0.0:
        values = -values
    return values / np.linalg.norm(values)


def quaternions_to_euler_zxz(quaternions: np.ndarray) -> np.ndarray:
    """Convert unit quaternions to Z-X-Z Euler angles in radians."""
    from scipy.spatial.transform import Rotation

    values = np.asarray(quaternions, dtype=float)
    if values.ndim == 1:
        values = values.reshape(1, 4)
    rotation = Rotation.from_quat(values)
    return rotation.as_euler("zxz", degrees=False)
