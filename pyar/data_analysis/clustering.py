import itertools
import logging
import operator
import os
import numpy as np
import pandas as pd
from sklearn.cluster import DBSCAN, KMeans
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import StandardScaler
from sklearn.metrics.pairwise import rbf_kernel
import hdbscan
import pyar.property
import pyar.representations

cluster_logger = logging.getLogger('pyar.cluster')

def remove_similar(list_of_molecules):
    final_list = list_of_molecules[:]
    cluster_logger.debug('Number of molecules before similarity elimination,  {}'.format(len(final_list)))
    for a, b in itertools.combinations(list_of_molecules, 2):
        energy_difference = calc_energy_difference(a, b)
        fingerprint_distance = calc_fingerprint_distance(a, b)
        if abs(energy_difference) < 1e-5 and abs(fingerprint_distance) < 1.0:
            if energy_difference < 0:
                if a in final_list:
                    cluster_logger.debug('Removing {}'.format(a.name))
                    final_list.remove(a)
            else:
                if b in final_list:
                    cluster_logger.debug('Removing {}'.format(b.name))
                    final_list.remove(b)
    cluster_logger.debug('Number of molecules after similarity elimination,  {}'.format(len(final_list)))
    print_energy_table(final_list)
    return final_list

def calc_energy_difference(a, b):
    return float(a.energy) - float(b.energy)

def calc_fingerprint_distance(a, b):
    """Calculate the distance between two fingerprints"""
    return np.linalg.norm(
        pyar.representations.fingerprint(a.atoms_list, a.coordinates)
        - pyar.representations.fingerprint(b.atoms_list, b.coordinates)
    )


def choose_geometries(list_of_molecules, maximum_number_of_seeds=12):
    if len(list_of_molecules) < 2:
        cluster_logger.info("Not enough data to cluster (only %d), returning original" % len(list_of_molecules))
        return list_of_molecules

    if len(list_of_molecules) <= maximum_number_of_seeds:
        cluster_logger.info('Not enough data for clustering. Removing similar geometries from the list')
        return remove_similar(list_of_molecules)

    # Read clustering algorithm from environment variable, default to 'hdbscan'
    algorithm = os.environ.get('PYAR_CLUSTERING_ALGORITHM', 'hdbscan').lower()
    cluster_logger.info(f'Clustering on {len(list_of_molecules)} geometries using {algorithm}')

    # Use MBTR for feature representation
    dt = np.array([pyar.representations.mbtr_descriptor(m.atoms_list, m.coordinates) for m in list_of_molecules])

    # Scale the data
    scaler = StandardScaler()
    dt_scaled = scaler.fit_transform(dt)

    # Save features to CSV
    pd.DataFrame(dt_scaled).to_csv("mbtr_features.csv")

    try:
        labels = generate_labels(dt_scaled, algorithm, maximum_number_of_seeds)
    except Exception as e:
        cluster_logger.exception(f"Clustering algorithm {algorithm} failed")
        cluster_logger.exception(e)
        return list_of_molecules

    best_from_each_cluster = select_best_from_each_cluster(labels, list_of_molecules)

    if len(best_from_each_cluster) == 1:
        return best_from_each_cluster
    else:
        cluster_logger.info("Removing similar molecules after clustering.")
        reduced_best_from_each_cluster = remove_similar(best_from_each_cluster)

    if len(reduced_best_from_each_cluster) > maximum_number_of_seeds:
        return choose_geometries(reduced_best_from_each_cluster, maximum_number_of_seeds)
    else:
        return reduced_best_from_each_cluster

def generate_labels(dt, algorithm='hdbscan', maximum_number_of_seeds=8):
    if algorithm == 'kmeans':
        return kmeans_clustering(dt, maximum_number_of_seeds)
    elif algorithm == 'dbscan':
        return dbscan_clustering(dt)
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

def hdbscan_clustering(dt):
    clusterer = hdbscan.HDBSCAN(min_cluster_size=2, min_samples=1)
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
    
    # Handle the case where there are negative labels (noise points)
    if np.any(labels < 0):
        cluster_logger.info("Clustering algorithm identified noise points.")
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

    cluster_logger.info("Lowest energy structures from each cluster:")
    print_energy_table(best_from_each_cluster)
    return best_from_each_cluster

def get_the_best_molecule(list_of_molecules):
    return min(list_of_molecules, key=lambda m: m.energy)

def print_energy_table(molecules):
    e_dict = {i.name: float(i.energy) for i in molecules}
    if len(e_dict) > 1:
        ref = min(e_dict.values())
        cluster_logger.info(f"{'Name':>35}:{'Energy':>12}{'R. E. (kcal/mol)':>18}")
        for name, energy in sorted(e_dict.items(), key=operator.itemgetter(1), reverse=True):
            cluster_logger.info(f"{name:>35}:{energy:12.6f}{(energy - ref) * 627.51:12.2f}")
        cluster_logger.info("")



def read_energy_from_xyz_file(xyz_file):
    import re
    try:
        with open(xyz_file, 'r') as fr:
            lines = fr.readlines()
            second_line = lines[1].strip()
            energy_str = re.search(r'-?\d+\.?\d*', second_line).group()
            energy = float(energy_str)
    except (IndexError, ValueError):
        with open(xyz_file, 'r') as fr:
            comments_line = fr.readlines()[1].rstrip()
            energy = float(re.split(':|=|\s+', comments_line)[1])

    return energy


def plot_energy_histogram(molecules):
    energies = [i.energy for i in molecules]
    ref = min(energies)
    relative_energies = [(float(energy) - ref) for energy in energies]
    histogram_bin = np.linspace(0, max(relative_energies), 10)
    import matplotlib.pyplot as plt
    plt.hist(relative_energies, histogram_bin)
    plt.xlabel('Energy')
    plt.ylabel('Population')
    plt.title('Histogram of energies')

    
def main():
    pass

if __name__ == "__main__":
    main()

