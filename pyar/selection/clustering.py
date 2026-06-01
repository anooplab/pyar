"""Selection orchestration and clustering algorithms.

This module owns the main ``choose_geometries`` entrypoint and the clustering
algorithm implementations. Shared helpers are provided by the specialized
modules under :mod:`pyar.selection`.
"""

import logging
import os
import numpy as np
import pyar.representations
from pyar.optional_dependencies import optional_dependency_error

cluster_logger = logging.getLogger('pyar.cluster')
_MBTR_RUNTIME_DISABLED = False
_MBTR_DISABLE_REASON = None

__all__ = [
    "affinity_propagation_clustering",
    "agglomerative_clustering",
    "choose_geometries",
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


def _require_sklearn_preprocessing(name):
    try:
        import sklearn.preprocessing as preprocessing
    except ImportError as exc:
        raise optional_dependency_error("sklearn", feature="selection clustering") from exc
    return getattr(preprocessing, name)


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
):
    global _MBTR_RUNTIME_DISABLED, _MBTR_DISABLE_REASON
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

    # Read selection algorithm from environment variable, default to the
    # hybrid cluster-plus-fill behavior when no explicit algorithm is provided.
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

    if _MBTR_RUNTIME_DISABLED:
        cluster_logger.info(
            "MBTR clustering disabled for this run (%s). Falling back to similarity pruning.",
            _MBTR_DISABLE_REASON,
        )
        selected = _limit_seed_count(
            pruned_molecules,
            maximum_number_of_seeds,
            reason="MBTR-disabled fallback",
        )
        return _finalize_selection(selected, write_path, existing_entries=basin_entries)

    # Use MBTR for feature representation. Some DScribe/ASE combinations can
    # fail during conversion, so disable MBTR for the current run after first
    # failure and use similarity pruning thereafter.
    try:
        dt = np.array([pyar.representations.mbtr_descriptor(m.atoms_list, m.coordinates) for m in pruned_molecules])
    except Exception as exc:
        _MBTR_RUNTIME_DISABLED = True
        _MBTR_DISABLE_REASON = str(exc)
        cluster_logger.warning(
            "MBTR clustering unavailable, disabling MBTR for this run and falling back to similarity pruning: %s",
            exc,
        )
        selected = _limit_seed_count(
            pruned_molecules,
            maximum_number_of_seeds,
            reason="MBTR failure fallback",
        )
        return _finalize_selection(selected, write_path, existing_entries=basin_entries)

    # Scale descriptor features before diversity/clustering selection.
    StandardScaler = _require_sklearn_preprocessing("StandardScaler")
    scaler = StandardScaler()
    dt_scaled = scaler.fit_transform(dt)

    if os.environ.get("PYAR_SAVE_MBTR_FEATURES") == "1":
        pd = _require_pandas()
        pd.DataFrame(dt_scaled).to_csv("mbtr_features.csv")

    if algorithm in {"maxmin", "max-min", "max_min"}:
        selected = _limit_seed_count(
            _max_min_diversity_select(dt_scaled, pruned_molecules, maximum_number_of_seeds),
            maximum_number_of_seeds,
            reason="max-min selection",
        )
        return _finalize_selection(selected, write_path, existing_entries=basin_entries)

    try:
        cluster_algorithm = "hdbscan" if algorithm == "hybrid" else algorithm
        labels = generate_labels(dt_scaled, cluster_algorithm, maximum_number_of_seeds)
    except Exception as exc:
        cluster_logger.warning(
            "Clustering algorithm %s failed (%s); falling back to max-min selection.",
            algorithm,
            exc,
        )
        cluster_logger.debug("Clustering failure traceback", exc_info=True)
        selected = _limit_seed_count(
            _max_min_diversity_select(
                dt_scaled,
                pruned_molecules,
                maximum_number_of_seeds,
            ),
            maximum_number_of_seeds,
            reason=f"{algorithm} failure fallback",
        )
        return _finalize_selection(selected, write_path, existing_entries=basin_entries)

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
        )
        selected = _limit_seed_count(
            trimmed,
            maximum_number_of_seeds,
            reason="cluster minima trimming",
        )
        return _finalize_selection(selected, write_path, existing_entries=basin_entries)

    if len(best_from_each_cluster) == maximum_number_of_seeds:
        selected = _limit_seed_count(
            best_from_each_cluster,
            maximum_number_of_seeds,
            reason="cluster selection",
        )
        return _finalize_selection(selected, write_path, existing_entries=basin_entries)

    selected_ids = {id(m) for m in best_from_each_cluster}
    if len(selected_ids) == len(pruned_molecules):
        return _finalize_selection(best_from_each_cluster, write_path, existing_entries=basin_entries)

    cluster_logger.info(
        "Cluster selection returned %d seeds; filling remaining %d with max-min.",
        len(best_from_each_cluster),
        maximum_number_of_seeds - len(best_from_each_cluster),
    )
    selected_indices = [
        index for index, molecule in enumerate(pruned_molecules)
        if id(molecule) in selected_ids
    ]
    filler = _max_min_diversity_select(
        dt_scaled,
        pruned_molecules,
        maximum_number_of_seeds,
        initial_selected_indices=selected_indices,
    )
    combined = best_from_each_cluster + filler
    selected = _limit_seed_count(
        combined,
        maximum_number_of_seeds,
        reason="hybrid cluster-fill selection",
    )
    return _finalize_selection(selected, write_path, existing_entries=basin_entries)

def generate_labels(dt, algorithm='hdbscan', maximum_number_of_seeds=8):
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
        cluster_logger.warning(f"Unknown algorithm: {algorithm}. Using HDBSCAN.")
        return hdbscan_clustering(dt)

def kmeans_clustering(dt, n_clusters):
    KMeans = _require_sklearn_cluster("KMeans")
    kmeans = KMeans(n_clusters=n_clusters, random_state=42)
    return kmeans.fit_predict(dt)

def dbscan_clustering(dt):
    DBSCAN = _require_sklearn_cluster("DBSCAN")
    eps, min_samples = determine_dbscan_params(dt)
    dbscan = DBSCAN(eps=eps, min_samples=min_samples)
    return dbscan.fit_predict(dt)


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
    AgglomerativeClustering = _require_sklearn_cluster("AgglomerativeClustering")
    n_clusters = max(2, min(n_clusters, len(dt)))
    clusterer = AgglomerativeClustering(n_clusters=n_clusters, linkage='ward')
    return clusterer.fit_predict(dt)

def gaussian_mixture_clustering(dt, n_components):
    GaussianMixture = _require_sklearn_mixture("GaussianMixture")
    gm = GaussianMixture(n_components=n_components, random_state=42)
    return gm.fit_predict(dt)

def rbf_kernel_clustering(dt, threshold=0.99):
    rbf_kernel = _require_sklearn_pairwise("rbf_kernel")
    similarities = rbf_kernel(dt)
    n_samples = similarities.shape[0]
    labels = np.zeros(n_samples, dtype=int)
    current_label = 0
    for i in range(n_samples):
        if labels[i] == 0:
            current_label += 1
            labels[i] = current_label
            for j in range(i+1, n_samples):
                if similarities[i, j] > threshold:
                    labels[j] = current_label
    return labels

def determine_dbscan_params(dt):
    # Simple heuristic for DBSCAN parameters
    distances = np.sort(np.sum((dt[:, None, :] - dt[None, :, :]) ** 2, axis=-1), axis=1)
    eps = np.median(distances[:, 1])
    min_samples = 2
    return eps, min_samples

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
