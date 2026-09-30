"""Validated pairwise distance backends for structural feature matrices."""

from __future__ import annotations

import numpy as np

DISTANCE_METRICS = ("euclidean", "manhattan", "cosine")


def pairwise_distances(values, metric="euclidean"):
    """Return a finite, symmetric distance matrix for a feature matrix."""
    values = np.asarray(values, dtype=float)
    metric = str(metric).strip().lower()
    if metric not in DISTANCE_METRICS:
        raise ValueError(f"Unknown distance metric {metric!r}. Choose one of: {', '.join(DISTANCE_METRICS)}")
    if values.ndim != 2 or not np.isfinite(values).all():
        raise ValueError("Distance input must be a finite two-dimensional matrix")
    if metric == "euclidean":
        differences = values[:, None, :] - values[None, :, :]
        result = np.linalg.norm(differences, axis=-1)
    elif metric == "manhattan":
        result = np.sum(np.abs(values[:, None, :] - values[None, :, :]), axis=-1)
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


__all__ = ["DISTANCE_METRICS", "pairwise_distances"]
