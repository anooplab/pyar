"""Validated cluster-label algorithms and deterministic fallback handling."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from pyar.selection.distances import DISTANCE_METRICS, pairwise_distances

CLUSTERING_ALGORITHMS = ("hybrid", "hdbscan", "agglomerative", "dbscan", "optics")
_ALGORITHM_ALIASES = {"hybrid": "hdbscan", "ward": "agglomerative"}


@dataclass(frozen=True)
class ClusteringResult:
    """Cluster labels plus the exact feature and algorithm path used."""

    labels: np.ndarray
    feature_requested: str
    feature_used: str
    algorithm_requested: str
    algorithm_used: str
    distance_metric: str
    feature_values: np.ndarray
    feature_fallbacks: tuple[dict[str, str], ...] = ()
    algorithm_fallbacks: tuple[dict[str, str], ...] = ()

    @property
    def number_of_clusters(self):
        return len({int(label) for label in self.labels if int(label) >= 0})

    @property
    def number_of_noise_points(self):
        return int(np.count_nonzero(np.asarray(self.labels) < 0))

    def to_dict(self):
        return {
            "feature_requested": self.feature_requested,
            "feature_used": self.feature_used,
            "feature_fallbacks": list(self.feature_fallbacks),
            "algorithm_requested": self.algorithm_requested,
            "algorithm_used": self.algorithm_used,
            "distance_metric": self.distance_metric,
            "algorithm_fallbacks": list(self.algorithm_fallbacks),
            "number_of_structures": int(len(self.labels)),
            "number_of_clusters": self.number_of_clusters,
            "number_of_noise_points": self.number_of_noise_points,
            "labels": [int(label) for label in self.labels],
        }


def determine_dbscan_params(values, min_samples=2, metric="euclidean"):
    """Estimate DBSCAN epsilon from k-neighbour distances in Euclidean units."""
    values = np.asarray(values, dtype=float)
    if values.ndim != 2 or not np.isfinite(values).all():
        raise ValueError("DBSCAN requires a finite two-dimensional feature matrix")
    if metric not in DISTANCE_METRICS:
        raise ValueError(f"Unknown distance metric {metric!r}")
    if len(values) < 2:
        return 0.0, 1
    min_samples = max(2, min(int(min_samples), len(values)))
    distances = pairwise_distances(values, metric=metric)
    np.fill_diagonal(distances, np.inf)
    kth_nearest = np.sort(distances, axis=1)[:, min_samples - 2]
    finite = kth_nearest[np.isfinite(kth_nearest)]
    if not len(finite):
        return 0.0, min_samples
    epsilon = float(np.quantile(finite, 0.75))
    if epsilon <= 0.0:
        positive = finite[finite > 0.0]
        epsilon = float(np.min(positive)) if len(positive) else 0.0
    return epsilon, min_samples


def _cluster_agglomerative(values, maximum_number_of_clusters, metric="euclidean"):
    """Average-linkage hierarchy using SciPy, the required SciPy dependency."""
    from scipy.cluster.hierarchy import cut_tree, linkage
    from scipy.spatial.distance import squareform

    values = np.asarray(values, dtype=float)
    if len(values) < 2:
        return np.zeros(len(values), dtype=int)
    if metric not in DISTANCE_METRICS:
        raise ValueError(f"Unknown distance metric {metric!r}")
    condensed = squareform(pairwise_distances(values, metric=metric), checks=False)
    if not np.isfinite(condensed).all():
        raise ValueError("Agglomerative clustering received invalid distances")
    if np.all(condensed <= 1e-14):
        return np.zeros(len(values), dtype=int)
    tree = linkage(condensed, method="average", optimal_ordering=False)
    labels = cut_tree(
        tree,
        n_clusters=[max(1, min(int(maximum_number_of_clusters), len(values)))],
    ).reshape(-1)
    return labels.astype(int)


def _run_algorithm(values, algorithm, maximum_number_of_clusters, options):
    algorithm = _ALGORITHM_ALIASES.get(algorithm, algorithm)
    values = np.asarray(values, dtype=float)
    if len(values) < 2:
        return np.zeros(len(values), dtype=int)
    metric = str(options.get("metric", "euclidean")).lower()
    if metric not in DISTANCE_METRICS:
        raise ValueError(f"Unknown distance metric {metric!r}. Choose one of: {', '.join(DISTANCE_METRICS)}")
    if algorithm == "agglomerative":
        return _cluster_agglomerative(values, maximum_number_of_clusters, metric=metric)
    if algorithm == "hdbscan":
        import hdbscan

        min_cluster_size = int(options.get("min_cluster_size", 2))
        min_samples = options.get("min_samples")
        clusterer = hdbscan.HDBSCAN(
            min_cluster_size=max(2, min(min_cluster_size, len(values))),
            min_samples=min_samples,
            metric=metric,
        )
        return clusterer.fit_predict(values)
    if algorithm == "dbscan":
        from sklearn.cluster import DBSCAN

        min_samples = int(options.get("min_samples", 2))
        epsilon = options.get("eps")
        if epsilon is None:
            epsilon, min_samples = determine_dbscan_params(values, min_samples, metric=metric)
        return DBSCAN(eps=float(epsilon), min_samples=min_samples, metric=metric).fit_predict(values)
    if algorithm == "optics":
        from sklearn.cluster import OPTICS

        return OPTICS(
            min_samples=max(2, min(int(options.get("min_samples", 2)), len(values))),
            xi=float(options.get("xi", 0.05)),
            min_cluster_size=options.get("min_cluster_size", 2),
            metric=metric,
        ).fit_predict(values)
    raise ValueError(
        f"Unknown clustering algorithm {algorithm!r}. "
        f"Choose one of: {', '.join(CLUSTERING_ALGORITHMS)}"
    )


def cluster_molecules(
    molecules,
    *,
    feature="mbtr",
    algorithm="hybrid",
    maximum_number_of_clusters=12,
    feature_fallbacks=True,
    algorithm_fallbacks=True,
    algorithm_options=None,
    distance_metric="euclidean",
):
    """Cluster geometries and report every feature or algorithm fallback.

    HDBSCAN is the requested method for the backwards-compatible ``hybrid``
    name. If it is unavailable or returns only noise, average-linkage
    agglomerative clustering supplies an explicit operational fallback.
    """
    from pyar.selection.features import compute_feature_matrix, standardize_features

    molecules = list(molecules)
    if maximum_number_of_clusters < 1:
        raise ValueError("maximum_number_of_clusters must be at least 1")
    requested_algorithm = str(algorithm).strip().lower()
    canonical_algorithm = _ALGORITHM_ALIASES.get(requested_algorithm, requested_algorithm)
    if requested_algorithm not in CLUSTERING_ALGORITHMS and canonical_algorithm not in CLUSTERING_ALGORITHMS:
        raise ValueError(
            f"Unknown clustering algorithm {requested_algorithm!r}. "
            f"Choose one of: {', '.join(CLUSTERING_ALGORITHMS)}"
        )
    if not molecules:
        raise ValueError("Cannot cluster an empty structure pool")

    feature_result = compute_feature_matrix(
        molecules,
        feature,
        allow_fallbacks=feature_fallbacks,
    )
    values = standardize_features(feature_result.values)
    distance_metric = str(distance_metric).strip().lower()
    if distance_metric not in DISTANCE_METRICS:
        raise ValueError(f"Unknown distance metric {distance_metric!r}. Choose one of: {', '.join(DISTANCE_METRICS)}")
    options = dict(algorithm_options or {})
    options["metric"] = distance_metric
    attempts = [canonical_algorithm]
    if algorithm_fallbacks and canonical_algorithm != "agglomerative":
        attempts.append("agglomerative")
    failures = []
    for candidate in attempts:
        try:
            labels = np.asarray(
                _run_algorithm(values, candidate, maximum_number_of_clusters, options),
                dtype=int,
            )
            if labels.shape != (len(molecules),):
                raise ValueError(f"Clusterer returned labels with shape {labels.shape}")
            if len(molecules) > 1 and np.all(labels < 0):
                raise ValueError("Clusterer classified every structure as noise")
            return ClusteringResult(
                labels=labels,
                feature_requested=str(feature),
                feature_used=feature_result.name,
                algorithm_requested=requested_algorithm,
                algorithm_used=candidate,
                distance_metric=distance_metric,
                feature_values=feature_result.values,
                feature_fallbacks=feature_result.fallbacks,
                algorithm_fallbacks=tuple(failures),
            )
        except Exception as exc:
            failures.append({"algorithm": candidate, "reason": f"{type(exc).__name__}: {exc}"})
    # Preserve every geometry as its own cluster if both requested and
    # hierarchical labelers fail. Selection can then apply its explicit
    # maximum-budget trimming rule without silently treating structures as
    # equivalent or dropping all candidates.
    failures.append({
        "algorithm": "singleton-preservation",
        "reason": "all cluster-label algorithms failed; each structure kept in its own cluster",
    })
    return ClusteringResult(
        labels=np.arange(len(molecules), dtype=int),
        feature_requested=str(feature),
        feature_used=feature_result.name,
        algorithm_requested=requested_algorithm,
        algorithm_used="singleton-preservation",
        distance_metric=distance_metric,
        feature_values=feature_result.values,
        feature_fallbacks=feature_result.fallbacks,
        algorithm_fallbacks=tuple(failures),
    )


__all__ = [
    "CLUSTERING_ALGORITHMS",
    "DISTANCE_METRICS",
    "ClusteringResult",
    "cluster_molecules",
    "determine_dbscan_params",
]
