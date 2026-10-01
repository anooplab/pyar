"""Selection orchestration and clustering algorithms.

This module owns the main ``choose_geometries`` entrypoint and the clustering
algorithm implementations. Shared helpers are provided by the specialized
modules under :mod:`pyar.selection`.
"""

import logging
import os
from numbers import Integral
import numpy as np
from pyar.selection.clusterers import cluster_molecules, validate_clustering_options
from pyar.selection.features import compute_feature_matrix, standardize_features
from pyar.selection.policy import resolve_clustering_policy
from pyar.selection.structural_distances import validate_distance_options
from pyar.optional_dependencies import optional_dependency_error

cluster_logger = logging.getLogger('pyar.cluster')

__all__ = [
    "affinity_propagation_clustering",
    "agglomerative_clustering",
    "choose_geometries",
    "cluster_molecules",
    "dbscan_clustering",
    "determine_dbscan_params",
    "gaussian_mixture_clustering",
    "generate_labels",
    "get_the_best_molecule",
    "hdbscan_clustering",
    "kmeans_clustering",
    "mean_shift_clustering",
    "optics_clustering",
    "plot_energy_histogram",
    "print_energy_table",
    "rbf_kernel_clustering",
    "read_energy_from_xyz_file",
    "remove_similar",
    "select_best_from_each_cluster",
    "spectral_clustering",
]


def _require_hdbscan():
    try:
        import hdbscan
    except ImportError as exc:
        raise optional_dependency_error("hdbscan", feature="HDBSCAN selection") from exc
    return hdbscan


def _require_pandas():
    try:
        import pandas as pd
    except ImportError as exc:
        raise optional_dependency_error("pandas", feature="selection CSV output") from exc
    return pd


def _require_sklearn_cluster(name):
    try:
        import sklearn.cluster as cluster
    except ImportError as exc:
        raise optional_dependency_error("sklearn", feature="selection clustering") from exc
    return getattr(cluster, name)


def _require_sklearn_mixture(name):
    try:
        import sklearn.mixture as mixture
    except ImportError as exc:
        raise optional_dependency_error("sklearn", feature="selection clustering") from exc
    return getattr(mixture, name)


def _require_sklearn_pairwise(name):
    try:
        import sklearn.metrics.pairwise as pairwise
    except ImportError as exc:
        raise optional_dependency_error("sklearn", feature="selection clustering") from exc
    return getattr(pairwise, name)


def choose_geometries(
    list_of_molecules,
    maximum_number_of_seeds=12,
    persist_basin_memory=True,
    apply_basin_memory=True,
    algorithm=None,
    group_basin_by_stoichiometry=True,
    connectivity_policy="off",
    feature="auto",
    distance_metric="euclidean",
    algorithm_options=None,
    system_type="auto",
    diagnostics=None,
    distance_options=None,
):
    list_of_molecules = list(list_of_molecules)
    if (isinstance(maximum_number_of_seeds, bool)
            or not isinstance(maximum_number_of_seeds, Integral)
            or maximum_number_of_seeds < 1):
        raise ValueError("maximum_number_of_seeds must be a positive integer")
    normalized_connectivity_policy = connectivity_policy.lower()
    if normalized_connectivity_policy not in {"off", "prefer", "strict"}:
        raise ValueError(
            f"Unknown connectivity policy: {connectivity_policy!r}. "
            "Expected one of 'off', 'prefer', or 'strict'."
        )

    # ``maxmin`` is a cluster-minima trimming request. Clustering still
    # defines the eligible candidate set, as required by the seed policy.
    if algorithm is None:
        algorithm = os.environ.get('PYAR_CLUSTERING_ALGORITHM', 'auto')
    algorithm = algorithm.strip().lower()
    cluster_algorithm = "auto" if algorithm in {"maxmin", "max-min", "max_min"} else algorithm
    validate_clustering_options(
        cluster_algorithm, maximum_number_of_seeds, distance_metric, algorithm_options
    )
    resolve_clustering_policy(system_type, feature, cluster_algorithm)
    validate_distance_options(distance_options)
    if diagnostics is not None:
        diagnostics.clear()
        diagnostics.update({
            "selection_algorithm_requested": algorithm,
            "algorithm_used": "not-run", "feature_requested": feature,
            "feature_used": None, "distance_metric": distance_metric,
            "feature_fallbacks": [], "algorithm_fallbacks": [],
            "distance_requested": distance_metric, "distance_used": "not-run",
            "distance_units": None, "distance_fallbacks": [],
            "input_count": len(list_of_molecules),
        })
    cluster_logger.info(f'Seed selection on {len(list_of_molecules)} geometries using {algorithm}')

    pruned_molecules = remove_similar(list_of_molecules)
    basin_registry_path = _basin_registry_path(
        pruned_molecules[0],
        group_by_stoichiometry=group_basin_by_stoichiometry,
    ) if pruned_molecules and os.path.isdir('selected') else None
    write_path = basin_registry_path if persist_basin_memory else None
    basin_entries = []
    basin_memory_diagnostics = {"status": "disabled", "schema_version": 2}
    if not apply_basin_memory and not persist_basin_memory:
        basin_memory_diagnostics["status"] = "disabled-by-request"
    elif basin_registry_path is None:
        basin_memory_diagnostics["status"] = "no-registry-path"
    if basin_registry_path and (apply_basin_memory or persist_basin_memory):
        from pyar.selection.basin_memory import BasinMemoryError

        try:
            basin_entries = _load_basin_registry(basin_registry_path)
        except BasinMemoryError as exc:
            # A damaged or future-version registry must never be overwritten
            # or used to remove candidates. Continue with memory disabled.
            cluster_logger.warning(
                "Basin registry unavailable; retaining all candidates and leaving it untouched: %s",
                exc,
            )
            write_path = None
            basin_entries = []
            basin_memory_diagnostics.update({
                "status": "registry-unavailable-keep-all",
                "registry_path": basin_registry_path,
                "reason": str(exc),
            })
        else:
            basin_memory_diagnostics.update({
                "status": "loaded" if basin_entries else "empty-registry",
                "registry_path": basin_registry_path,
                "schema_version": 2,
            })
    if apply_basin_memory and basin_entries:
        pruned_molecules = _apply_basin_memory(
            pruned_molecules,
            maximum_number_of_seeds,
            basin_entries,
            feature=feature,
            distance_metric=distance_metric,
            system_type=system_type,
            algorithm=cluster_algorithm,
            distance_options=distance_options,
            diagnostics=basin_memory_diagnostics,
        )
    if normalized_connectivity_policy in {"prefer", "strict"}:
        pruned_molecules = _prefer_connected_structures(
            pruned_molecules,
            policy=normalized_connectivity_policy,
        )

    if diagnostics is not None:
        if basin_memory_diagnostics.get("status") in {"loaded", "empty-registry"}:
            basin_memory_diagnostics.setdefault("output_candidates", len(pruned_molecules))
            basin_memory_diagnostics.setdefault("input_candidates", len(list_of_molecules))
        diagnostics["basin_memory"] = basin_memory_diagnostics
        diagnostics["candidate_names"] = [molecule.name for molecule in pruned_molecules]
        diagnostics["candidate_paths"] = [
            str(getattr(molecule, "relative_path", molecule.name)) for molecule in pruned_molecules
        ]

    def finish(selected, stage):
        if diagnostics is not None:
            diagnostics.update({
                "selection_stage": stage,
                "selected_names": [molecule.name for molecule in selected],
                "selected_count": len(selected),
            })
        return _finalize_selection(selected, write_path, existing_entries=basin_entries)

    if len(pruned_molecules) <= maximum_number_of_seeds:
        cluster_logger.info(
            "Similarity pruning reduced seed pool to %d; skipping diversity selection.",
            len(pruned_molecules),
        )
        if len(pruned_molecules) < maximum_number_of_seeds:
            _log_seed_shortfall(maximum_number_of_seeds, len(pruned_molecules), "similarity/connectivity filtering")
        return finish(pruned_molecules, "within-budget-after-filtering")

    clustering_result = cluster_molecules(
        pruned_molecules,
        feature=feature,
        algorithm=cluster_algorithm,
        maximum_number_of_clusters=maximum_number_of_seeds,
        distance_metric=distance_metric,
        algorithm_options=algorithm_options,
        system_type=system_type,
        distance_options=distance_options,
    )
    _log_clustering_result(clustering_result)
    for fallback in clustering_result.distance_fallbacks:
        cluster_logger.warning(
            "Distance %s failed; using %s (%s).", fallback["distance"],
            clustering_result.distance_used, fallback["reason"],
        )
    if diagnostics is not None:
        diagnostics.update(clustering_result.to_dict())
    labels = clustering_result.labels
    dt_scaled = standardize_features(clustering_result.feature_values)
    if clustering_result.feature_values.shape[1]:
        _save_features_if_requested(clustering_result.feature_used, dt_scaled)

    best_from_each_cluster = select_best_from_each_cluster(labels, pruned_molecules)

    if len(best_from_each_cluster) > maximum_number_of_seeds:
        cluster_logger.info(
            "Cluster selection returned %d minima; trimming to %d with max-min.",
            len(best_from_each_cluster),
            maximum_number_of_seeds,
        )
        selected_ids = {id(molecule) for molecule in best_from_each_cluster}
        cluster_feature_indices = [
            index for index, molecule in enumerate(pruned_molecules)
            if id(molecule) in selected_ids
        ]
        cluster_subset_features = dt_scaled[cluster_feature_indices]
        cluster_subset_molecules = [pruned_molecules[index] for index in cluster_feature_indices]
        trimmed = _max_min_diversity_select(
            cluster_subset_features,
            cluster_subset_molecules,
            maximum_number_of_seeds,
            distance_metric=distance_metric,
            distance_matrix=(None if clustering_result.distance_matrix is None else
                             clustering_result.distance_matrix[np.ix_(cluster_feature_indices, cluster_feature_indices)]),
        )
        selected = _limit_seed_count(
            trimmed,
            maximum_number_of_seeds,
            reason="cluster minima trimming",
        )
        return finish(selected, "cluster-minima-max-min-trimming")

    # Keep only representatives justified by the clustering result. Max-min
    # is used above to trim an overfull set of cluster minima; it must not add
    # geometries from clusters that did not contribute a representative.
    selected = _limit_seed_count(
        best_from_each_cluster,
        maximum_number_of_seeds,
        reason="cluster selection",
    )
    return finish(selected, "cluster-minima")


def _log_feature_result(feature_result):
    cluster_logger.info(
        "Structural policy: system=%s confidence=%s; feature=%s (dimensions=%d)",
        feature_result.system_type,
        feature_result.system_confidence,
        feature_result.name,
        feature_result.values.shape[1],
    )
    for fallback in feature_result.fallbacks:
        cluster_logger.warning(
            "Feature %s failed; using fallback %s (%s).",
            fallback["feature"],
            feature_result.name,
            fallback["reason"],
        )


def _save_features_if_requested(feature_name, values):
    if os.environ.get("PYAR_SAVE_MBTR_FEATURES") != "1":
        return
    pd = _require_pandas()
    filename = "mbtr_features.csv" if feature_name == "mbtr" else f"{feature_name}_features.csv"
    pd.DataFrame(values).to_csv(filename)


def _log_clustering_result(result):
    cluster_logger.info(
        "Structural policy: system=%s confidence=%s; feature=%s; reason=%s",
        result.system_type,
        result.system_confidence,
        result.feature_used,
        result.policy_reason,
    )
    for fallback in result.feature_fallbacks:
        cluster_logger.warning(
            "Feature %s failed; using fallback %s (%s).",
            fallback["feature"],
            result.feature_used,
            fallback["reason"],
        )
    cluster_logger.info(
        "Clusterer: requested=%s used=%s distance=%s (requested=%s; units=%s) clusters=%d noise=%d",
        result.algorithm_requested,
        result.algorithm_used,
        result.distance_used,
        result.distance_metric,
        result.distance_units,
        result.number_of_clusters,
        result.number_of_noise_points,
    )
    for fallback in result.algorithm_fallbacks:
        cluster_logger.warning(
            "Clusterer %s failed; using fallback %s (%s).",
            fallback["algorithm"],
            result.algorithm_used,
            fallback["reason"],
        )

def generate_labels(dt, algorithm='hdbscan', maximum_number_of_seeds=8):
    """Compatibility wrapper for clustering precomputed feature vectors."""
    algorithm = str(algorithm).strip().lower()
    if algorithm in {'auto', 'hybrid', 'hdbscan'}:
        return hdbscan_clustering(dt)
    if algorithm in {'maxmin', 'max-min', 'max_min'}:
        raise ValueError("max-min is a selector, not a cluster-label algorithm")
    if algorithm == 'kmeans':
        return kmeans_clustering(dt, maximum_number_of_seeds)
    elif algorithm == 'dbscan':
        return dbscan_clustering(dt)
    elif algorithm == 'optics':
        return optics_clustering(dt)
    elif algorithm in {'affinity', 'affinity_propagation'}:
        return affinity_propagation_clustering(dt)
    elif algorithm in {'meanshift', 'mean_shift'}:
        return mean_shift_clustering(dt)
    elif algorithm in {'agglomerative', 'ward'}:
        return agglomerative_clustering(dt, maximum_number_of_seeds)
    elif algorithm == 'spectral':
        return spectral_clustering(dt, maximum_number_of_seeds)
    elif algorithm == 'hdbscan':
        return hdbscan_clustering(dt)
    elif algorithm == 'gaussian_mixture':
        return gaussian_mixture_clustering(dt, maximum_number_of_seeds)
    elif algorithm == 'rbf_kernel':
        return rbf_kernel_clustering(dt)
    else:
        from pyar.selection.clusterers import CLUSTERING_ALGORITHMS
        raise ValueError(
            f"Unknown algorithm: {algorithm!r}. Choose one of: "
            f"{', '.join(CLUSTERING_ALGORITHMS)}"
        )

def kmeans_clustering(dt, n_clusters):
    KMeans = _require_sklearn_cluster("KMeans")
    kmeans = KMeans(n_clusters=n_clusters, random_state=42)
    return kmeans.fit_predict(dt)

def dbscan_clustering(dt):
    from pyar.selection.clusterers import _run_algorithm
    return _run_algorithm(dt, "dbscan", 12, {})


def optics_clustering(dt):
    OPTICS = _require_sklearn_cluster("OPTICS")
    clusterer = OPTICS(min_samples=2, xi=0.05, min_cluster_size=2)
    return clusterer.fit_predict(dt)

def hdbscan_clustering(dt):
    hdbscan = _require_hdbscan()
    clusterer = hdbscan.HDBSCAN(min_cluster_size=2, min_samples=1)
    return clusterer.fit_predict(dt)


def affinity_propagation_clustering(dt):
    AffinityPropagation = _require_sklearn_cluster("AffinityPropagation")
    clusterer = AffinityPropagation(random_state=42)
    return clusterer.fit_predict(dt)


def mean_shift_clustering(dt):
    MeanShift = _require_sklearn_cluster("MeanShift")
    estimate_bandwidth = _require_sklearn_cluster("estimate_bandwidth")
    bandwidth = estimate_bandwidth(dt, quantile=0.2, n_samples=min(len(dt), 500))
    if not np.isfinite(bandwidth) or bandwidth <= 0:
        bandwidth = None
    clusterer = MeanShift(bandwidth=bandwidth, bin_seeding=True)
    return clusterer.fit_predict(dt)


def spectral_clustering(dt, n_clusters):
    SpectralClustering = _require_sklearn_cluster("SpectralClustering")
    n_clusters = max(2, min(n_clusters, len(dt)))
    n_neighbors = max(1, min(10, len(dt) - 1))
    clusterer = SpectralClustering(
        n_clusters=n_clusters,
        random_state=42,
        assign_labels='kmeans',
        affinity='nearest_neighbors',
        n_neighbors=n_neighbors,
    )
    return clusterer.fit_predict(dt)


def agglomerative_clustering(dt, n_clusters):
    from pyar.selection.clusterers import _cluster_agglomerative
    return _cluster_agglomerative(dt, n_clusters)

def gaussian_mixture_clustering(dt, n_components):
    GaussianMixture = _require_sklearn_mixture("GaussianMixture")
    gm = GaussianMixture(n_components=n_components, random_state=42)
    return gm.fit_predict(dt)

def rbf_kernel_clustering(dt, threshold=0.99):
    """Cluster connected components of the RBF similarity threshold graph."""
    if not 0.0 <= float(threshold) <= 1.0:
        raise ValueError("RBF similarity threshold must be between 0 and 1")
    rbf_kernel = _require_sklearn_pairwise("rbf_kernel")
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import connected_components

    similarities = rbf_kernel(dt)
    if not np.isfinite(similarities).all():
        raise ValueError("RBF kernel produced non-finite similarities")
    adjacency = similarities >= float(threshold)
    np.fill_diagonal(adjacency, True)
    _, labels = connected_components(csr_matrix(adjacency), directed=False)
    return labels.astype(int)

def determine_dbscan_params(dt):
    from pyar.selection.clusterers import determine_dbscan_params as estimate
    return estimate(dt)

# def select_best_from_each_cluster(labels, list_of_molecules):
#     unique_labels = np.unique(labels)
#     cluster_logger.info(f"The distribution of file in each cluster: {np.bincount(labels)}")
#     best_from_each_cluster = []
#     for this_label in unique_labels:
#         if this_label != -1:  # -1 is the noise label in some clustering algorithms
#             molecules_in_this_group = [m for label, m in zip(labels, list_of_molecules) if label == this_label]
#             best_from_each_cluster.append(get_the_best_molecule(molecules_in_this_group))
#     cluster_logger.info("Lowest energy structures from each cluster")
#     print_energy_table(best_from_each_cluster)
#     return best_from_each_cluster

def select_best_from_each_cluster(labels, list_of_molecules):
    labels = np.array(labels)  # Ensure labels is a numpy array
    unique_labels = np.unique(labels)
    noise_molecules = []

    # Handle the case where there are negative labels (noise points)
    if np.any(labels < 0):
        cluster_logger.info("Clustering algorithm identified noise points.")
        noise_molecules = [m for l, m in zip(labels, list_of_molecules) if l == -1]
        positive_labels = labels[labels >= 0]
        if len(positive_labels) > 0:
            cluster_logger.info(f"The distribution of files in each cluster (excluding noise): {np.bincount(positive_labels)}")
        else:
            cluster_logger.info("No valid clusters found.")
    else:
        cluster_logger.info(f"The distribution of files in each cluster: {np.bincount(labels)}")

    best_from_each_cluster = []
    for label in unique_labels:
        if label != -1:  # Exclude noise points (label -1)
            molecules_in_this_group = [m for l, m in zip(labels, list_of_molecules) if l == label]
            if molecules_in_this_group:
                best_from_each_cluster.append(get_the_best_molecule(molecules_in_this_group))

    if noise_molecules:
        cluster_logger.info(
            "Keeping all %d noise points as candidates; they may represent rare structures.",
            len(noise_molecules),
        )
        best_from_each_cluster.extend(noise_molecules)

    cluster_logger.info("Lowest energy structures from each cluster:")
    print_energy_table(best_from_each_cluster)
    return best_from_each_cluster

def get_the_best_molecule(list_of_molecules):
    return min(list_of_molecules, key=lambda m: m.energy)

# Shared selection helpers live in the focused service modules.
from pyar.selection.basin_memory import (  # noqa: E402
    BasinMemoryError,
    _apply_basin_memory,
    _basin_novelty_scores,
    _basin_registry_path,
    _entry_fingerprint,
    _fingerprint_signature,
    _load_basin_registry,
    _persist_basin_registry,
    _stoichiometry_label,
    migrate_basin_registry,
    record_selected_basins,
)
from pyar.selection.deduplication import (  # noqa: E402
    _adaptive_duplicate_rmsd_threshold,
    _assigned_element_order,
    _equivalent_atom_groups,
    _exact_element_orders,
    _iterative_assigned_rmsd,
    _kabsch_rotation,
    _kabsch_rmsd,
    _prefer_connected_structures,
    _rmsd_after_alignment,
    _structure_is_similar,
    calc_fingerprint_distance,
    remove_similar,
)
from pyar.selection.diversity import (  # noqa: E402
    _finalize_selection,
    _limit_seed_count,
    _log_seed_shortfall,
    _max_min_diversity_select,
)
from pyar.selection.reports import (  # noqa: E402
    plot_energy_histogram,
    print_energy_table,
    read_energy_from_xyz_file,
)


def main():
    pass

if __name__ == "__main__":
    main()
