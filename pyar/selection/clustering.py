"""Selection orchestration and clustering algorithms.

This module owns the main ``choose_geometries`` entrypoint and the clustering
algorithm implementations. Shared helpers are provided by the specialized
modules under :mod:`pyar.selection`.
"""

import logging
import os
import numpy as np
from pyar.selection.clusterers import cluster_molecules
from pyar.selection.features import compute_feature_matrix, standardize_features
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
    feature="mbtr",
    distance_metric="euclidean",
    algorithm_options=None,
):
    normalized_connectivity_policy = connectivity_policy.lower()
    if normalized_connectivity_policy not in {"off", "prefer", "strict"}:
        raise ValueError(
            f"Unknown connectivity policy: {connectivity_policy!r}. "
            "Expected one of 'off', 'prefer', or 'strict'."
        )

    if len(list_of_molecules) < 2:
        _log_seed_shortfall(maximum_number_of_seeds, len(list_of_molecules), "selection")
        basin_registry_path = _basin_registry_path(
            list_of_molecules[0],
            group_by_stoichiometry=group_basin_by_stoichiometry,
        ) if list_of_molecules and os.path.isdir('selected') else None
        write_path = basin_registry_path if persist_basin_memory else None
        return _finalize_selection(list_of_molecules, write_path)

    if len(list_of_molecules) <= maximum_number_of_seeds:
        cluster_logger.info('Not enough data for clustering. Removing similar geometries from the list')
        basin_registry_path = _basin_registry_path(
            list_of_molecules[0],
            group_by_stoichiometry=group_basin_by_stoichiometry,
        ) if os.path.isdir('selected') else None
        write_path = basin_registry_path if persist_basin_memory else None
        selected = _limit_seed_count(
            remove_similar(list_of_molecules),
            maximum_number_of_seeds,
            reason="similarity pruning",
        )
        return _finalize_selection(selected, write_path)

    # Preserve the existing selector policy while allowing the representation
    # and clustering algorithm to be selected independently.
    if algorithm is None:
        algorithm = os.environ.get('PYAR_CLUSTERING_ALGORITHM', 'hybrid')
    algorithm = algorithm.lower()
    cluster_logger.info(f'Seed selection on {len(list_of_molecules)} geometries using {algorithm}')

    pruned_molecules = remove_similar(list_of_molecules)
    basin_registry_path = _basin_registry_path(
        pruned_molecules[0],
        group_by_stoichiometry=group_basin_by_stoichiometry,
    ) if pruned_molecules and os.path.isdir('selected') else None
    write_path = basin_registry_path if persist_basin_memory else None
    basin_entries = _load_basin_registry(basin_registry_path) if apply_basin_memory else []
    if basin_entries:
        pruned_molecules = _apply_basin_memory(pruned_molecules, maximum_number_of_seeds, basin_entries)
    if normalized_connectivity_policy in {"prefer", "strict"}:
        pruned_molecules = _prefer_connected_structures(
            pruned_molecules,
            policy=normalized_connectivity_policy,
        )

    if len(pruned_molecules) <= maximum_number_of_seeds:
        cluster_logger.info(
            "Similarity pruning reduced seed pool to %d; skipping diversity selection.",
            len(pruned_molecules),
        )
        return _finalize_selection(pruned_molecules, write_path, existing_entries=basin_entries)

    if algorithm in {"maxmin", "max-min", "max_min"}:
        feature_result = compute_feature_matrix(pruned_molecules, feature)
        dt_scaled = standardize_features(feature_result.values)
        _log_feature_result(feature_result)
        _save_features_if_requested(feature_result.name, dt_scaled)
        selected = _limit_seed_count(
            _max_min_diversity_select(
                dt_scaled,
                pruned_molecules,
                maximum_number_of_seeds,
                distance_metric=distance_metric,
            ),
            maximum_number_of_seeds,
            reason="max-min selection",
        )
        return _finalize_selection(selected, write_path, existing_entries=basin_entries)

    clustering_result = cluster_molecules(
        pruned_molecules,
        feature=feature,
        algorithm=algorithm,
        maximum_number_of_clusters=maximum_number_of_seeds,
        distance_metric=distance_metric,
        algorithm_options=algorithm_options,
    )
    _log_clustering_result(clustering_result)
    labels = clustering_result.labels
    dt_scaled = standardize_features(clustering_result.feature_values)
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
        )
        selected = _limit_seed_count(
            trimmed,
            maximum_number_of_seeds,
            reason="cluster minima trimming",
        )
        return _finalize_selection(selected, write_path, existing_entries=basin_entries)

    # Keep only representatives justified by the clustering result. Max-min
    # is used above to trim an overfull set of cluster minima; it must not add
    # geometries from clusters that did not contribute a representative.
    selected = _limit_seed_count(
        best_from_each_cluster,
        maximum_number_of_seeds,
        reason="cluster selection",
    )
    return _finalize_selection(selected, write_path, existing_entries=basin_entries)


def _log_feature_result(feature_result):
    cluster_logger.info(
        "Structural feature: %s (dimensions=%d)",
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
    cluster_logger.info("Structural feature: %s", result.feature_used)
    for fallback in result.feature_fallbacks:
        cluster_logger.warning(
            "Feature %s failed; using fallback %s (%s).",
            fallback["feature"],
            result.feature_used,
            fallback["reason"],
        )
    cluster_logger.info(
        "Clusterer: requested=%s used=%s distance=%s clusters=%d noise=%d",
        result.algorithm_requested,
        result.algorithm_used,
        result.distance_metric,
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
    if algorithm in {'hybrid', 'hdbscan'}:
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
        best_noise = get_the_best_molecule(noise_molecules)
        cluster_logger.info(
            "Including best noise-point representative: %s (energy %.6f)",
            best_noise.name,
            float(best_noise.energy),
        )
        best_from_each_cluster.append(best_noise)

    cluster_logger.info("Lowest energy structures from each cluster:")
    print_energy_table(best_from_each_cluster)
    return best_from_each_cluster

def get_the_best_molecule(list_of_molecules):
    return min(list_of_molecules, key=lambda m: m.energy)

# Shared selection helpers live in the focused service modules.
from pyar.selection.basin_memory import (  # noqa: E402
    _apply_basin_memory,
    _basin_novelty_scores,
    _basin_registry_path,
    _entry_fingerprint,
    _fingerprint_signature,
    _load_basin_registry,
    _persist_basin_registry,
    _stoichiometry_label,
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
