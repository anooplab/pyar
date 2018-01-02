import itertools
import collections
import operator
import sys

import numpy as np
from sklearn.cluster import DBSCAN, AffinityPropagation, KMeans, MeanShift, estimate_bandwidth


def get_labels(data_as_list, algorithm='meanshift'):
    dt = np.array(data_as_list)
    labels = []

    print('    Algorithm =', algorithm)

    if algorithm == 'dbscan':
        dbs = DBSCAN(eps=0.1)
        dbs.fit(dt)
        labels = dbs.labels_

    if algorithm == 'kmeans':
        kmeans = KMeans(n_clusters=10)
        kmeans.fit(dt)
        labels = kmeans.labels_

    if algorithm == 'meanshift':
        # The following bandwidth can be automatically detected using
        try:
            bandwidth = estimate_bandwidth(dt, quantile=0.2, n_samples=len(dt))
        except:
            bandwidth = 0.5
        ms = MeanShift(bandwidth=bandwidth, bin_seeding=True)
        ms.fit(dt)
        labels = ms.labels_

    if algorithm == 'affinitypropagation':
        af = AffinityPropagation()
        af.fit(dt)
        labels = af.labels_

    return labels


def remove_similar(dict_of_molecules):
    bfec = dict_of_molecules.copy()

    for a, b in itertools.combinations(dict_of_molecules, 2):
        energy_difference = dict_of_molecules[a].energy - dict_of_molecules[b].energy
        fingerprint_distance = np.linalg.norm(dict_of_molecules[a].fingerprint - dict_of_molecules[b].fingerprint)
        if abs(energy_difference) < 1e-5 and abs(fingerprint_distance) < 1.0:
            if energy_difference < 1:
                if a in bfec:
                    del bfec[a]
            else:
                if b in bfec:
                    del bfec[b]
    return bfec


def choose_geometries(t_molecules):
    if len(t_molecules) < 2:
        print("    Not enough data to cluster, returning original")
        return t_molecules

    print('   Clustering on ', len(t_molecules), 'geometries')
    # print_energy_table(t_molecules)

    molecules_sorted = sorted(t_molecules.keys())
    #    list_of_coulomb_matrices = [t_molecules[i].sorted_coulomb_matrix for i in sorted(t_molecules.keys())]
    # list_of_fingerprints = [t_molecules[i].fingerprint for i in molecules_sorted]
    # dt = list_of_fingerprints

    dt = [(t_molecules[i]. energy, np.sum(t_molecules[i].fingerprint)) for i in molecules_sorted]

    dt = np.around(dt, decimals=5)

    import pandas as pd
    df = pd.DataFrame(dt)
    df.to_csv('dataframe.csv')

    try:
        labels = generate_labels(dt)
    except:
        print("   All Clustering algorithms failed")
        return t_molecules.values()

    best_from_each_cluster = select_best_from_each_cluster(labels, t_molecules)

    if len(best_from_each_cluster) == 1:
        return best_from_each_cluster
    else:
        reduced_best_from_each_cluster = remove_similar(best_from_each_cluster)

    print("After removing similar ones, the lowest energy structures from each cluster")
    print_energy_table(reduced_best_from_each_cluster)

    if len(reduced_best_from_each_cluster) > 8:
        return choose_geometries(reduced_best_from_each_cluster)
    else:
        return reduced_best_from_each_cluster


def print_energy_table(molecules):
    e_dict = {i: molecules[i].energy for i in molecules}
    ref = min(e_dict.values())
    for name, energy in sorted(e_dict.items(), key=operator.itemgetter(1), reverse=True):
        print("      {:>15}:{:12.6f}{:12.2f}".format(name, energy, (energy - ref) * 627.51))


def generate_labels(dt):
    # type: (np.array) -> list
    """
    Try many algorithms for clutering one after the other.
    :param dt: data for clustering as np.array
    :return: list of labels
    """
    try:
        return get_labels(dt, algorithm='meanshift')
    except:
        print('    MeanShift failed')
    try:
        return get_labels(dt, algorithm='affinitypropagation')
    except:
        print('    Affinity Propagation failed')
    try:
        return get_labels(dt, algorithm='dbscan')
    except:
        print('    DBSCAN failed')

    return get_labels(dt, algorithm='kmeans')


def select_best_from_each_cluster(labels, t_molecules):
    # type: (list, dict) -> dict
    """
    :rtype: dict
    :type t_molecules dict of molecules
    :type labels list
    """
    molecules_sorted = sorted(t_molecules.keys())
    unique_labels = np.unique(labels)
    print("   The distribution of file in each cluster:", np.bincount(labels[labels >= 0]))
    print("   The distribution of file in each cluster:", collections.Counter(labels))

    best_from_each_cluster = {}
    for n in unique_labels:
        best_in_this_cluster = {}
        for label, k in zip(labels, molecules_sorted):
            if label == n:
                best_in_this_cluster[k] = t_molecules[k].energy
        best = min(best_in_this_cluster, key=best_in_this_cluster.get)
        best_from_each_cluster[best] = t_molecules[best]
    print("    Lowest energy structures from each cluster")
    print_energy_table(best_from_each_cluster)
    return best_from_each_cluster


def read_energy_from_xyz_file(xyzfile):
    with open(xyzfile, 'r') as fr:
        comments_line = fr.readlines()[1]
    energy = float(comments_line.split(':')[1])
    return energy


# main program
def main():
    from Molecule import Molecule
    try:
        input_files = sys.argv[1:]
    except:
        print('usage: cluster.;y <xyz-file(s)>')
        sys.exit(1)

    if len(input_files) < 2:
        print('Not enough files to cluster')
        sys.exit(0)
    mols = {}
    for each_file in input_files:
        mol = Molecule.from_xyz(each_file)
        mol.energy = read_energy_from_xyz_file(each_file)
        mols[mol.name] = mol
    choose_geometries(mols)
    return


if __name__ == "__main__":
    main()
