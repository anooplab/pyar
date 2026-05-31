"""Compatibility shim for the legacy ``pyar.data_analysis.clustering`` path.

Use :mod:`pyar.selection.clustering` instead.
"""

from pyar.selection.clustering import (
    affinity_propagation_clustering,
    agglomerative_clustering,
    calc_fingerprint_distance,
    choose_geometries,
    cluster_logger,
    dbscan_clustering,
    determine_dbscan_params,
    gaussian_mixture_clustering,
    generate_labels,
    get_the_best_molecule,
    hdbscan_clustering,
    kmeans_clustering,
    mean_shift_clustering,
    optics_clustering,
    plot_energy_histogram,
    print_energy_table,
    rbf_kernel_clustering,
    read_energy_from_xyz_file,
    remove_similar,
    select_best_from_each_cluster,
    spectral_clustering,
)

__all__ = [
    "affinity_propagation_clustering",
    "agglomerative_clustering",
    "calc_fingerprint_distance",
    "choose_geometries",
    "cluster_logger",
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
