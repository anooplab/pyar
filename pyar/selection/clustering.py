"""Clustering algorithms and cluster representatives for selection."""

from __future__ import annotations

import logging

import hdbscan
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from sklearn.cluster import (
    DBSCAN,
    OPTICS,
    AffinityPropagation,
    AgglomerativeClustering,
    KMeans,
    MeanShift,
    SpectralClustering,
    estimate_bandwidth,
)
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import StandardScaler
from sklearn.metrics.pairwise import rbf_kernel

import pyar.representations


cluster_logger = logging.getLogger('pyar.cluster')


def determine_dbscan_params(dt):
    # Simple heuristic for DBSCAN parameters
    distances = np.sort(np.sum((dt[:, None, :] - dt[None, :, :]) ** 2, axis=-1), axis=1)
    eps = np.median(distances[:, 1])
    min_samples = 2
    return eps, min_samples


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
    kmeans = KMeans(n_clusters=n_clusters, random_state=42)
    return kmeans.fit_predict(dt)


def dbscan_clustering(dt):
    eps, min_samples = determine_dbscan_params(dt)
    dbscan = DBSCAN(eps=eps, min_samples=min_samples)
    return dbscan.fit_predict(dt)


def optics_clustering(dt):
    clusterer = OPTICS(min_samples=2, xi=0.05, min_cluster_size=2)
    return clusterer.fit_predict(dt)


def hdbscan_clustering(dt):
    clusterer = hdbscan.HDBSCAN(min_cluster_size=2, min_samples=1)
    return clusterer.fit_predict(dt)


def affinity_propagation_clustering(dt):
    clusterer = AffinityPropagation(random_state=42)
    return clusterer.fit_predict(dt)


def mean_shift_clustering(dt):
    bandwidth = estimate_bandwidth(dt, quantile=0.2, n_samples=min(len(dt), 500))
    if not np.isfinite(bandwidth) or bandwidth <= 0:
        bandwidth = None
    clusterer = MeanShift(bandwidth=bandwidth, bin_seeding=True)
    return clusterer.fit_predict(dt)


def spectral_clustering(dt, n_clusters):
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
    n_clusters = max(2, min(n_clusters, len(dt)))
    clusterer = AgglomerativeClustering(n_clusters=n_clusters, linkage='ward')
    return clusterer.fit_predict(dt)


def gaussian_mixture_clustering(dt, n_components):
    gm = GaussianMixture(n_components=n_components, random_state=42)
    return gm.fit_predict(dt)


def rbf_kernel_clustering(dt, threshold=0.99):
    similarities = rbf_kernel(dt)
    n_samples = similarities.shape[0]
    labels = np.zeros(n_samples, dtype=int)
    current_label = 0
    for i in range(n_samples):
        if labels[i] == 0:
            current_label += 1
            labels[i] = current_label
            for j in range(i + 1, n_samples):
                if similarities[i, j] > threshold:
                    labels[j] = current_label
    return labels


def get_the_best_molecule(list_of_molecules):
    return min(list_of_molecules, key=lambda m: m.energy)


def print_energy_table(molecules, stream=None, title=None):
    """Report energies with relative values against the global minimum."""
    e_dict = {i.name: float(i.energy) for i in molecules}
    if not e_dict:
        return

    ref = min(e_dict.values())
    lines = []
    if title:
        lines.append(title)
    lines.append(f"{'Name':>35}  {'Energy':>12}  {'R. E. (kcal/mol)':>18}")
    for name, energy in sorted(e_dict.items(), key=lambda item: item[1]):
        lines.append(f"{name:>35}  {energy:12.6f}  {(energy - ref) * 627.51:18.2f}")
    lines.append(f"Global minimum: {min(e_dict, key=e_dict.get)} ({ref:12.6f} Eh)")
    lines.append("")

    if stream is not None:
        for line in lines:
            print(line, file=stream)
    else:
        for line in lines:
            cluster_logger.info(line)


def select_best_from_each_cluster(labels, list_of_molecules):
    labels = np.array(labels)
    unique_labels = np.unique(labels)
    noise_molecules = []

    if np.any(labels < 0):
        cluster_logger.info("Clustering algorithm identified noise points.")
        noise_molecules = [m for l, m in zip(labels, list_of_molecules) if l == -1]
        positive_labels = labels[labels >= 0]
        if len(positive_labels) > 0:
            cluster_logger.info(
                f"The distribution of files in each cluster (excluding noise): {np.bincount(positive_labels)}"
            )
        else:
            cluster_logger.info("No valid clusters found.")
    else:
        cluster_logger.info(f"The distribution of files in each cluster: {np.bincount(labels)}")

    best_from_each_cluster = []
    for label in unique_labels:
        if label != -1:
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
