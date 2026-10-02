"""Validated cluster-label algorithms and deterministic fallback handling."""

from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral

import numpy as np
from pyar.selection.distances import (
    DISTANCE_METRICS, FEATURE_DISTANCE_METRICS, pairwise_distances, validate_distance_matrix,
)
from pyar.selection.structural_distances import compute_distance_matrix, validate_distance_options

# ``hybrid`` remains a hidden compatibility alias for old API/state files.
CLUSTERING_ALGORITHMS = ("auto", "hdbscan", "agglomerative", "dbscan", "optics")
_ALGORITHM_ALIASES = {"auto": "hdbscan", "hybrid": "hdbscan", "ward": "agglomerative"}


def validate_clustering_options(algorithm, maximum_number_of_clusters, distance_metric, options=None):
    """Reject invalid user settings before operational fallbacks can hide them."""
    if (isinstance(maximum_number_of_clusters, bool)
            or not isinstance(maximum_number_of_clusters, Integral)
            or maximum_number_of_clusters < 1):
        raise ValueError("maximum_number_of_clusters must be a positive integer")
    algorithm = str(algorithm).strip().lower()
    if algorithm not in (*CLUSTERING_ALGORITHMS, "hybrid", "ward"):
        raise ValueError(f"Unknown clustering algorithm {algorithm!r}")
    metric = str(distance_metric).strip().lower()
    if metric not in DISTANCE_METRICS:
        raise ValueError(f"Unknown distance metric {metric!r}")
    options = dict(options or {})
    unknown = set(options) - {"min_samples", "min_cluster_size", "eps", "xi", "metric"}
    if unknown:
        raise ValueError(f"Unknown clustering options: {', '.join(sorted(unknown))}")
    for name, minimum in (("min_samples", 1), ("min_cluster_size", 2)):
        value = options.get(name)
        if value is not None and (
            isinstance(value, bool) or not isinstance(value, Integral) or value < minimum
        ):
            raise ValueError(f"{name} must be an integer >= {minimum}")
    epsilon = options.get("eps")
    if epsilon is not None and (not np.isfinite(epsilon) or epsilon <= 0):
        raise ValueError("eps must be finite and positive")
    xi = options.get("xi")
    if xi is not None and (not np.isfinite(xi) or not 0 < xi < 1):
        raise ValueError("xi must be finite and between 0 and 1")
    options = {key: value for key, value in options.items() if value is not None}
    options["metric"] = metric
    return options


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
    system_type: str = "unknown"
    system_confidence: str = "low"
    policy_reason: str = ""
    topology_group_ids: tuple[int, ...] = ()
    distance_matrix: np.ndarray | None = None
    distance_used: str = ""
    distance_units: str = "standardized-feature-units"
    distance_parameters: dict | None = None
    distance_fallbacks: tuple[dict, ...] = ()

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
            "distance_requested": self.distance_metric,
            "distance_used": self.distance_used or self.distance_metric,
            "distance_units": self.distance_units,
            "distance_parameters": self.distance_parameters or {},
            "distance_fallbacks": list(self.distance_fallbacks),
            "algorithm_fallbacks": list(self.algorithm_fallbacks),
            "system_type": self.system_type,
            "system_confidence": self.system_confidence,
            "policy_reason": self.policy_reason,
            "topology_group_ids": list(self.topology_group_ids),
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
    if metric not in (*FEATURE_DISTANCE_METRICS, "precomputed"):
        raise ValueError(f"Unknown distance metric {metric!r}")
    if len(values) < 2:
        return 0.0, 1
    min_samples = max(2, min(int(min_samples), len(values)))
    distances = (validate_distance_matrix(values).copy() if metric == "precomputed"
                 else pairwise_distances(values, metric=metric))
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
    if metric not in (*FEATURE_DISTANCE_METRICS, "precomputed"):
        raise ValueError(f"Unknown distance metric {metric!r}")
    distances = validate_distance_matrix(values) if metric == "precomputed" else pairwise_distances(values, metric=metric)
    condensed = squareform(distances, checks=False)
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
    if metric not in (*FEATURE_DISTANCE_METRICS, "precomputed"):
        raise ValueError(f"Unknown distance metric {metric!r}")
    if metric == "precomputed":
        values = validate_distance_matrix(values)
    if algorithm == "agglomerative":
        return _cluster_agglomerative(values, maximum_number_of_clusters, metric=metric)
    if algorithm == "hdbscan":
        import hdbscan

        min_cluster_size = int(options.get("min_cluster_size", 2))
        min_samples = options.get("min_samples")
        # HDBSCAN's tree backend cannot evaluate cosine. A shared distance
        # matrix also preserves PyAR's explicit convention for zero vectors.
        data = pairwise_distances(values, metric=metric) if metric == "cosine" else values
        clusterer = hdbscan.HDBSCAN(
            min_cluster_size=max(2, min(min_cluster_size, len(values))),
            min_samples=min_samples,
            metric="precomputed" if metric == "cosine" else metric,
        )
        return clusterer.fit_predict(data)
    if algorithm == "dbscan":
        from sklearn.cluster import DBSCAN

        min_samples = int(options.get("min_samples", 2))
        epsilon = options.get("eps")
        if epsilon is None:
            epsilon, min_samples = determine_dbscan_params(values, min_samples, metric=metric)
        data = pairwise_distances(values, metric=metric) if metric == "cosine" else values
        # sklearn requires eps > 0 even when every point has zero distance.
        if epsilon == 0.0:
            epsilon = np.finfo(float).eps
        return DBSCAN(
            eps=float(epsilon), min_samples=min_samples,
            metric="precomputed" if metric == "cosine" else metric,
        ).fit_predict(data)
    if algorithm == "optics":
        from sklearn.cluster import OPTICS

        data = pairwise_distances(values, metric=metric) if metric == "cosine" else values
        return OPTICS(
            min_samples=max(2, min(int(options.get("min_samples", 2)), len(values))),
            xi=float(options.get("xi", 0.05)),
            min_cluster_size=options.get("min_cluster_size", 2),
            metric="precomputed" if metric == "cosine" else metric,
        ).fit_predict(data)
    raise ValueError(
        f"Unknown clustering algorithm {algorithm!r}. "
        f"Choose one of: {', '.join(CLUSTERING_ALGORITHMS)}"
    )


def _run_grouped_algorithm(values, algorithm, maximum_number_of_clusters, options, groups):
    """Cluster each inferred molecular/aggregate topology independently."""
    labels = np.full(len(values), -1, dtype=int)
    next_label = 0
    for group_id in sorted(set(int(group) for group in groups)):
        indices = np.flatnonzero(np.asarray(groups) == group_id)
        local = (
            np.zeros(len(indices), dtype=int)
            if len(indices) == 1
            else np.asarray(
                _run_algorithm(
                    values[np.ix_(indices, indices)] if options.get("metric") == "precomputed" else values[indices],
                    algorithm, maximum_number_of_clusters, options
                ),
                dtype=int,
            )
        )
        if local.shape != (len(indices),):
            raise ValueError(f"Clusterer returned labels with shape {local.shape}")
        for local_label in sorted(int(label) for label in np.unique(local) if label >= 0):
            labels[indices[local == local_label]] = next_label
            next_label += 1
    return labels


def cluster_molecules(
    molecules,
    *,
    feature="auto",
    algorithm="auto",
    maximum_number_of_clusters=12,
    feature_fallbacks=True,
    algorithm_fallbacks=True,
    algorithm_options=None,
    distance_metric="euclidean",
    system_type="auto",
    distance_options=None,
    distance_fallbacks=True,
):
    """Cluster geometries and report every feature or algorithm fallback.

    ``auto`` tries HDBSCAN first. If it is unavailable or returns only noise,
    average-linkage agglomerative clustering supplies an operational fallback.
    ``hybrid`` is accepted only as a compatibility alias for older callers.
    """
    from pyar.selection.features import compute_feature_matrix, standardize_features

    molecules = list(molecules)
    requested_algorithm = str(algorithm).strip().lower()
    options = validate_clustering_options(
        requested_algorithm, maximum_number_of_clusters, distance_metric, algorithm_options
    )
    canonical_algorithm = _ALGORITHM_ALIASES.get(requested_algorithm, requested_algorithm)
    if requested_algorithm not in CLUSTERING_ALGORITHMS and canonical_algorithm not in CLUSTERING_ALGORITHMS:
        raise ValueError(
            f"Unknown clustering algorithm {requested_algorithm!r}. "
            f"Choose one of: {', '.join(CLUSTERING_ALGORITHMS)}"
        )
    if not molecules:
        raise ValueError("Cannot cluster an empty structure pool")

    distance_options = validate_distance_options(distance_options)
    distance_metric = str(distance_metric).strip().lower()
    distance_result = None
    if distance_metric in FEATURE_DISTANCE_METRICS:
        feature_result = compute_feature_matrix(
            molecules, feature, allow_fallbacks=feature_fallbacks,
            system_type=system_type, algorithm=requested_algorithm,
        )
        values = standardize_features(feature_result.values)
    else:
        from pyar.selection.features import FeatureMatrix
        from pyar.selection.policy import classify_system_pool, resolve_clustering_policy

        resolve_clustering_policy(system_type, feature, requested_algorithm)
        classification = classify_system_pool(molecules, system_type)
        feature_result = None

        def fallback_features():
            nonlocal feature_result
            feature_result = compute_feature_matrix(
                molecules, feature, allow_fallbacks=feature_fallbacks,
                system_type=system_type, algorithm=requested_algorithm,
            )
            return standardize_features(feature_result.values)

        distance_result = compute_distance_matrix(
            molecules, distance_metric, feature_values=fallback_features,
            allow_fallbacks=distance_fallbacks, options=distance_options,
        )
        if feature_result is None:
            feature_result = FeatureMatrix(
                "local-soap" if distance_result.used == "soap-rematch" else "coordinates",
                np.empty((len(molecules), 0)),
                tuple(sorted({atom for molecule in molecules for atom in molecule.atoms_list})),
                system_type=classification.system_type,
                system_confidence=classification.confidence,
                policy_reason=classification.reason,
                topology_group_ids=classification.topology_group_ids,
            )
        values = distance_result.values
        options["metric"] = "precomputed"
        if distance_result.used != distance_result.requested and "eps" in options:
            # A radius in Angstrom cannot silently become a kernel/descriptor
            # radius. Re-estimate DBSCAN epsilon in the successful backend's scale.
            distance_result.parameters["requested_eps_discarded"] = options.pop("eps")
            distance_result.parameters["epsilon_policy"] = "reestimate-after-distance-fallback"
    distance_record = {
        "distance_matrix": None if distance_result is None else distance_result.values,
        "distance_used": distance_metric if distance_result is None else distance_result.used,
        "distance_units": ("dimensionless" if distance_metric == "cosine" else "standardized-feature-units")
                          if distance_result is None else distance_result.units,
        "distance_parameters": {} if distance_result is None else distance_result.parameters,
        "distance_fallbacks": () if distance_result is None else distance_result.fallbacks,
    }
    attempts = [canonical_algorithm]
    if algorithm_fallbacks and canonical_algorithm != "agglomerative":
        attempts.append("agglomerative")
    failures = []
    if (len(molecules) > 1 and not np.any(values)
            and (distance_result is None or distance_result.used not in {"graph-rmsd", "fragment-rmsd"})):
        failures.append({
            "algorithm": canonical_algorithm,
            "reason": "descriptor has no variation; cluster similarity cannot be established",
        })
        attempts = []
    for candidate in attempts:
        try:
            labels = np.asarray(
                _run_grouped_algorithm(
                    values,
                    candidate,
                    maximum_number_of_clusters,
                    options,
                    feature_result.topology_group_ids,
                )
                if feature_result.topology_group_ids
                else _run_algorithm(values, candidate, maximum_number_of_clusters, options),
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
                system_type=feature_result.system_type,
                system_confidence=feature_result.system_confidence,
                policy_reason=feature_result.policy_reason,
                topology_group_ids=feature_result.topology_group_ids,
                **distance_record,
            )
        except Exception as exc:
            failures.append({"algorithm": candidate, "reason": f"{type(exc).__name__}: {exc}"})
    # Preserve every geometry as its own cluster if both requested and
    # hierarchical labelers fail. Selection can then apply its explicit
    # maximum-budget trimming rule without silently treating structures as
    # equivalent or dropping all candidates.
    failures.append({
        "algorithm": "singleton-preservation",
        "reason": "cluster labels could not be established; each structure kept in its own cluster",
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
        system_type=feature_result.system_type,
        system_confidence=feature_result.system_confidence,
        policy_reason=feature_result.policy_reason,
        topology_group_ids=feature_result.topology_group_ids,
        **distance_record,
    )


__all__ = [
    "CLUSTERING_ALGORITHMS",
    "DISTANCE_METRICS",
    "ClusteringResult",
    "cluster_molecules",
    "determine_dbscan_params",
]
