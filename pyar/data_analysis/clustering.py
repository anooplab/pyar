import itertools
import logging
import operator

import numpy as np
import pandas as pd
from sklearn.cluster import DBSCAN, AffinityPropagation, KMeans, \
    MiniBatchKMeans, MeanShift, estimate_bandwidth
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import RobustScaler

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


def memoize(f):
    """ Memoization decorator for functions taking one or more arguments.
    https://code.activestate.com/recipes/578231-probably-the-fastest-memoization-decorator-in-the-/
    """

    class MemoDict(dict):
        # noinspection PyCompatibility
        def __init__(self, fun):
            super().__init__()
            self.f = fun

        def __call__(self, *args):
            return self[args]

        def __missing__(self, key):
            ret = self[key] = self.f(*key)
            return ret

    return MemoDict(f)


@memoize
def calc_energy_difference(a, b):
    return a.energy - b.energy


@memoize
def calc_fingerprint_distance(a, b):
    """Calculate the distance between two fingerprints"""
    return np.linalg.norm(
        pyar.representations.fingerprint(a.atoms_list, a.coordinates)
        - pyar.representations.fingerprint(b.atoms_list, b.coordinates)
    )


def choose_geometries(list_of_molecules, features='fingerprint', maximum_number_of_seeds=8):
    if len(list_of_molecules) < 2:
        cluster_logger.info("    Not enough data to cluster (only %d), returning original" % len(list_of_molecules))
        return list_of_molecules

    if len(list_of_molecules) <= maximum_number_of_seeds:
        cluster_logger.info('    Not enough data for clustering. '
                            '    Removing similar geometries from the list')
        return remove_similar(list_of_molecules)

    cluster_logger.info('Clustering on {} geometries'.format(len(list_of_molecules)))

    if features == 'fingerprint':
        dt = [pyar.representations.fingerprint(i.atoms_list, i.coordinates) for i in list_of_molecules]
    elif features == 'scm':
        dt = [
            pyar.representations.sorted_coulomb_matrix(pyar.representations.coulomb_matrix(i.atoms_list, i.coordinates))
            for i in list_of_molecules]
    elif features == 'moi':
        dt = [pyar.property.get_principal_axes(i.moments_of_inertia_tensor) for i in list_of_molecules]
    elif features == 'rsmd':
        dt = [pyar.representations.get_rsmd(i.moments_of_inertia_tensor) for i in list_of_molecules]
    else:
        cluster_logger.error('This feature is not implemented')
        return list_of_molecules

    dt = np.around(dt, decimals=5)

    df = pd.DataFrame(dt)
    df.to_csv("features.csv")

    scale_it = RobustScaler()
    dt = scale_it.fit_transform(dt)

    try:
        labels = generate_labels(dt)
    except Exception as e:
        cluster_logger.exception("All Clustering algorithms failed")
        cluster_logger.exception(e)
        return list_of_molecules

    best_from_each_cluster = select_best_from_each_cluster(labels,
                                                           list_of_molecules)

    if len(best_from_each_cluster) == 1:
        return best_from_each_cluster
    else:
        cluster_logger.info("    Removing similar molecules after clustering.")
        reduced_best_from_each_cluster = remove_similar(best_from_each_cluster)

    if len(reduced_best_from_each_cluster) > maximum_number_of_seeds:
        return choose_geometries(reduced_best_from_each_cluster,
                                 maximum_number_of_seeds=maximum_number_of_seeds)
    else:
        return reduced_best_from_each_cluster


def print_energy_table(molecules):
    e_dict = {i.name: i.energy for i in molecules}
    if len(e_dict) > 1:
        ref = min(e_dict.values())
        cluster_logger.info(f"     {'Name':>35}:{'Energy':>12}{'R. E. (kcal/mol)':>18}")
        for name, energy in sorted(e_dict.items(), key=operator.itemgetter(1), reverse=True):
            cluster_logger.info(f"     {name:>35}:{energy:12.6f}{(energy - ref) * 627.51:12.2f}")
        cluster_logger.info(f"")


def get_labels(data_as_list, algorithm='combo'):
    dt = np.array(data_as_list)
    labels = []

    cluster_logger.debug(' Algorithm = {}'.format(algorithm))

    if algorithm == 'combo':
        kmeans_labels, centres = n_clusters_optimized_with_kmeans(dt)
        cluster_logger.info('    Clustering with MeahShift algorithm using seeds '
                            'from K-Means')
        ms = MeanShift(bandwidth=None, bin_seeding=True, seeds=centres)
        ms.fit(dt)
        meanshift_labels = ms.labels_
        n_labels = len(np.unique(meanshift_labels))
        labels = meanshift_labels if 1 < n_labels < 8 else kmeans_labels
    if algorithm == 'meanshift':

        try:
            # The following bandwidth can be automatically detected using
            bandwidth = estimate_bandwidth(dt, quantile=0.2, n_samples=len(dt))
            cluster_logger.debug(' Estimated bandwidth: {}'.format(bandwidth))
        except Exception as e:
            # If it fails to estimate bandwidth, give a value (it is arbitrary
            # now,
            cluster_logger.exception(' Estimation of bandwidth failed. Using '
                                     '0.5 as bandwidth')
            cluster_logger.exception(e)

            bandwidth = 0.5
        try:
            ms = MeanShift(bandwidth=bandwidth, bin_seeding=True)
            ms.fit(dt)
            labels = ms.labels_
        except Exception as e:
            cluster_logger.error('MeanShift failed')

    if algorithm == 'dbscan':
        dbs = DBSCAN(eps=0.1)
        try:
            dbs.fit(dt)
            labels = dbs.labels_
        except Exception:
            cluster_logger.exception('DBSCAN Failed')

    if algorithm == 'kmeans':
        kmeans = KMeans(n_clusters=8)
        try:
            kmeans.fit(dt)
            labels = kmeans.labels_
        except Exception:
            cluster_logger.exception('K-Means  Failed')

    if algorithm == 'affinitypropagation':
        af = AffinityPropagation()
        try:
            af.fit(dt)
            labels = af.labels_
        except Exception:
            cluster_logger.exception('Affinity Propagation Failed')

    return labels


def n_clusters_optimized_with_kmeans(dt):
    if len(dt) > 10000:
        kmeans = MiniBatchKMeans(n_clusters=8)
        kmeans.fit(dt)
        labels = kmeans.labels_
        centres = kmeans.cluster_centers_
        cluster_logger.info('    For more than 10k samples, KMeans with '
                            'n_clusters=8 is used as memory becomes a problem.')
        return labels, centres

    labels = {}
    centres = {}
    scores = {}
    for i in range(2, min(len(dt), 9)):
        kmeans = KMeans(n_clusters=i)
        try:
            kmeans.fit(dt)
            labels[i] = kmeans.labels_
            centres[i] = kmeans.cluster_centers_
            scores[i] = silhouette_score(dt, labels[i])
            cluster_logger.debug('n_clusters: {}; score: {}'.format(i, scores[i]))
        except Exception as e:
            cluster_logger.error('K-Means failed')
            cluster_logger.error(e)
    if scores:
        best = max(scores, key=scores.get)
        cluster_logger.info('Best was {} clusters with Silhouette score of {'
                            '}'.format(best, scores[best]))
        return labels[best], centres[best]
    else:
        return [0 for _ in range(len(dt))], [1.e+00]


def generate_labels(dt):
    # type: (np.array) -> list
    """
    Try many algorithms for clustering one after the other.
    :param dt: data for clustering as np.array
    :return: list of labels
    """
    methods = ['combo', 'meahshift', 'affinitypropagation', 'dbscan']
    for method in methods:
        try:
            return get_labels(dt, algorithm=method)
        except Exception as e:
            cluster_logger.exception(f'{method} failed\n{e}')


def select_best_from_each_cluster(labels, list_of_molecules):
    unique_labels = np.unique(labels)
    cluster_logger.info("    The distribution of file in each cluster {}:".format(np.bincount(labels)))
    best_from_each_cluster = []
    for this_label in unique_labels:
        molecules_in_this_group = [m for label, m in
                                   zip(labels, list_of_molecules)
                                   if label == this_label]
        best_from_each_cluster.append(get_the_best_molecule(molecules_in_this_group))
    cluster_logger.info("    Lowest energy structures from each cluster")
    print_energy_table(best_from_each_cluster)
    return best_from_each_cluster


def get_the_best_molecule(list_of_molecules):
    """Give a list of molecule objects with name and energy, and this function
    returns the molecule with lower energy
    :type list_of_molecules: list of Molecule objects
    """
    energy_dict = {each_molecule.name: each_molecule.energy
                   for each_molecule in list_of_molecules}
    key_of_the_least_value = min(energy_dict, key=energy_dict.get)
    for i in list_of_molecules:
        if i.name == key_of_the_least_value:
            return i


def read_energy_from_xyz_file(xyz_file):
    import re
    with open(xyz_file, 'r') as fr:
        comments_line = fr.readlines()[1].rstrip()
    return float(re.split(':|=|\s+', comments_line)[1])


def plot_energy_histogram(molecules):
    energies = [i.energy for i in molecules]
    ref = min(energies)
    relative_energies = [(energy - ref) * 627.51 for energy in energies]
    histogram_bin = np.linspace(0, max(relative_energies), 10)
    import matplotlib.pyplot as plt
    plt.hist(relative_energies, histogram_bin)
    plt.xlabel('Energy')
    plt.ylabel('Population')
    plt.title('Histogram of energies')
    plt.show()


# main program
def main():
    pass


if __name__ == "__main__":
    main()
