"""Validated pairwise distance backends for structural feature matrices."""

from __future__ import annotations

import numpy as np
from scipy.spatial.distance import cdist

FEATURE_DISTANCE_METRICS = ("euclidean", "manhattan", "cosine")
STRUCTURAL_DISTANCE_METRICS = ("graph-rmsd", "fragment-rmsd", "soap-rematch")
DISTANCE_METRICS = FEATURE_DISTANCE_METRICS + STRUCTURAL_DISTANCE_METRICS


def pairwise_distances(values, metric="euclidean"):
    """Return a finite, symmetric distance matrix for a feature matrix."""
    values = np.asarray(values, dtype=float)
    metric = str(metric).strip().lower()
    if metric not in FEATURE_DISTANCE_METRICS:
        raise ValueError(f"Unknown distance metric {metric!r} for feature vectors. Choose one of: {', '.join(FEATURE_DISTANCE_METRICS)}")
    if values.ndim != 2 or not np.isfinite(values).all():
        raise ValueError("Distance input must be a finite two-dimensional matrix")
    if metric == "euclidean":
        result = cdist(values, values, metric="euclidean")
    elif metric == "manhattan":
        result = cdist(values, values, metric="cityblock")
    else:
        norms = np.linalg.norm(values, axis=1)
        denominator = norms[:, None] * norms[None, :]
        similarity = np.divide(
            values @ values.T,
            denominator,
            out=np.zeros((len(values), len(values)), dtype=float),
            where=denominator > 0,
        )
        zero_rows = norms == 0
        zero_pairs = denominator == 0
        both_zero = zero_rows[:, None] & zero_rows[None, :]
        similarity[zero_pairs] = both_zero[zero_pairs].astype(float)
        result = 1.0 - np.clip(similarity, -1.0, 1.0)
    result = np.asarray(result, dtype=float)
    np.fill_diagonal(result, 0.0)
    if not np.isfinite(result).all():
        raise ValueError("Distance backend produced non-finite values")
    return result


def validate_distance_matrix(values, size=None):
    """Validate a complete dissimilarity matrix without inventing missing pairs."""
    values = np.asarray(values, dtype=float)
    if values.ndim != 2 or values.shape[0] != values.shape[1]:
        raise ValueError("Distance matrix must be square")
    if size is not None and values.shape != (size, size):
        raise ValueError("Distance matrix size must match the structure pool")
    if not np.isfinite(values).all() or np.any(values < 0):
        raise ValueError("Distance matrix must be finite and nonnegative")
    if not np.allclose(values, values.T, rtol=1e-10, atol=1e-12):
        raise ValueError("Distance matrix must be symmetric")
    if not np.allclose(np.diag(values), 0, rtol=0, atol=1e-12):
        raise ValueError("Distance matrix must have a zero diagonal")
    return values


__all__ = ["DISTANCE_METRICS", "FEATURE_DISTANCE_METRICS", "STRUCTURAL_DISTANCE_METRICS",
           "pairwise_distances", "validate_distance_matrix"]
