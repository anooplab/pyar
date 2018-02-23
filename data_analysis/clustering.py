import itertools
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
        try:
            # The following bandwidth can be automatically detected using
            bandwidth = estimate_bandwidth(dt, quantile=0.2, n_samples=len(dt))
        except:
            # If it fails to estimate bandwidth, give a value (it is arbitrary now,
            # TODO can we find a reasonable value?)
            bandwidth = 0.5
        ms = MeanShift(bandwidth=bandwidth, bin_seeding=True)
        ms.fit(dt)
        labels = ms.labels_

    if algorithm == 'affinitypropagation':
        af = AffinityPropagation()
        af.fit(dt)
        labels = af.labels_

    return labels


def remove_similar(list_of_molecules):
    final_list = list_of_molecules[:]

    for a, b in itertools.combinations(list_of_molecules, 2):
        energy_difference = a.energy - b.energy
        fingerprint_distance = np.linalg.norm(a.fingerprint - b.fingerprint)
        if abs(energy_difference) < 1e-5 and abs(fingerprint_distance) < 1.0:
            if energy_difference < 1:
                if a in final_list:
                    final_list.remove(a)
            else:
                if b in final_list:
                    final_list.remove(b)
    return final_list


def choose_geometries(list_of_molecules):
    if len(list_of_molecules) < 2:
        print("    Not enough data to cluster (only ", len(list_of_molecules), ") , returning original")
        return list_of_molecules

    print('   Clustering on ', len(list_of_molecules), 'geometries')

    # dt = [i.fingerprint for i in list_of_molecules]
    dt = [i.sorted_coulomb_matrix for i in list_of_molecules]

    dt = np.around(dt, decimals=5)

    import pandas as pd
    df = pd.DataFrame(dt)
    df.to_csv("features.csv")

    try:
        labels = generate_labels(dt)
    except:
        print("   All Clustering algorithms failed")
        return list_of_molecules

    best_from_each_cluster = select_best_from_each_cluster(labels, list_of_molecules)

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
    e_dict = {i.name: i.energy for i in molecules}
    ref = min(e_dict.values())
    for name, energy in sorted(e_dict.items(), key=operator.itemgetter(1), reverse=True):
        print("      {:>15}:{:12.6f}{:12.2f}".format(name, energy, (energy - ref) * 627.51))


def generate_labels(dt):
    # type: (np.array) -> list
    """
    Try many algorithms for clustering one after the other.
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


def select_best_from_each_cluster(labels, list_of_molecules):
    unique_labels = np.unique(labels)
    print("   The distribution of file in each cluster:", np.bincount(labels[labels >= 0]))

    best_from_each_cluster = []
    for this_label in unique_labels:
        mols_in_this_grp = [m for label, m in zip(labels, list_of_molecules) if label == this_label]
        best_from_each_cluster.append(best_molecule(mols_in_this_grp))
    print("    Lowest energy structures from each cluster")
    print_energy_table(best_from_each_cluster)
    return best_from_each_cluster


def best_molecule(list_of_molecules):
    energy_dict = {each_molecule.name: each_molecule.energy for each_molecule in list_of_molecules}
    key_of_the_least_value = min(energy_dict, key=energy_dict.get)
    return list_of_molecules[key_of_the_least_value]


def read_energy_from_xyz_file(xyz_file):
    with open(xyz_file, 'r') as fr:
        comments_line = fr.readlines()[1]
    energy = float(comments_line.split(':')[-1])
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
    return


if __name__ == "__main__":
    main()
