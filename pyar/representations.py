import itertools
from itertools import product

import numpy as np

import pyar.property
from pyar.Molecule import Molecule


def get_rsmd(mol):
    """Molecular Descriptor written by Rajat Shubhro Majumdar"""
    mol.move_to_centre_of_mass()
    centre_of_mass = pyar.property.get_centre_of_mass(mol.coordinates, mol.atomic_mass)
    number_of_atoms = mol.number_of_atoms

    radii = np.array([pyar.property.distance(c, centre_of_mass) for c in mol.coordinates])
    one_s = np.sum(np.exp(-1 * radii / mol.atomic_number))
    two_s = np.sum(np.exp(-1 * radii / mol.atomic_number) * (2 - radii / mol.atomic_number))
    p_ex = 0.0
    p_yi = 0.0
    p_zd = 0.0
    zc = [[z, c] for z, c in zip(mol.atomic_number, mol.coordinates)]
    for i, j in itertools.permutations(zc, 2):
        dist = pyar.property.distance(i[1], j[1])
        z1 = i[0]
        z2 = j[0]
        if dist < 4.0:
            a = i[1]
            b = centre_of_mass
            c = j[1]
            v1 = a - b
            v2 = c - b
            cos_theta = (np.dot(v1, v2) / (
                    np.linalg.norm(v1) *
                    np.linalg.norm(v2)))
            cos_theta = min(cos_theta, 1.)
            cos_theta = max(cos_theta, -1.)
            sin_theta = np.sqrt(1 - (cos_theta ** 2))
            r = np.linalg.norm(a)
            p_ex += (r / (z1 * z2)) * (np.exp(-r / z1 * z2)) * cos_theta
            p_yi += (r / (z1 * z2)) * (np.exp(-r / z1 * z2)) * sin_theta
            p_zd += ((r / z1 * z2) ** 2) * (3 * (cos_theta ** 2) - 1)

    return np.array([number_of_atoms, one_s, two_s, p_ex, p_yi, p_zd])


def make_internal_coordinates(mol):
    """


    :rtype: list of internal coordinates
    """
    coordinates = mol.coordinates
    covalent_radius = mol.covalent_radius
    dm = pyar.property.get_distance_matrix(coordinates)
    bm = pyar.property.get_bond_matrix(coordinates, covalent_radius)
    bl = []
    for a_i in range(len(coordinates)):
        for a_j in range(len(coordinates)):
            if bm[a_i, a_j]:
                seq = [a_i, a_j]
                if seq[::-1] not in [i[0] for i in bl]:
                    bl.append([seq, dm[a_i, a_j]])
    al = []
    for a_i, a_j, k in product(range(len(coordinates)), repeat=3):
        if a_i != a_j and a_j != k and a_i != k and bm[a_i, a_j] and bm[a_j, k]:
            seq = [a_i, a_j, k]
            if seq[::-1] not in [i[0] for i in al]:
                al.append([seq, pyar.property.calculate_angle(coordinates[a_i], coordinates[a_j],
                                                              coordinates[k])])

    dl = []
    for a_i, a_j, k, l in product(range(len(coordinates)), repeat=4):
        if (
                a_i != a_j
                and a_i != k
                and a_i != l
                and a_j != k
                and a_j != l
                and k != l
                and bm[a_i, a_j]
                and bm[a_j, k]
                and bm[k, l]
        ):
            seq = [a_i, a_j, k, l]
            if seq[::-1] not in [i[0] for i in dl]:
                dl.append([seq, pyar.property.calculate_dihedral(a_i, a_j, k, l)])

    return [bl, al, dl]


def sorted_coulomb_matrix(cm):
    """
    From: https://github.com/pythonpanda/coulomb_matrix/
    Takes in a Coulomb matrix of (mxn) dimension and performs a
    row-wise sorting such that ||C(j,:)|| > ||C(j+1,:)||, J=
    0,1,.......,(m-1) Finally returns a vectorized (m*n,1) column
    matrix.

    :param cm: Coulomb matrix
    :type cm: ndarray
    """
    summation = np.array([sum(x ** 2) for x in cm])
    sorted_matrix = cm[np.argsort(summation)[::-1, ], :]
    return sorted_matrix.ravel()


def coulomb_matrix(atoms_list, coordinates):
    """

    :return: Coulomb Matrix

    """
    from pyar.data import new_atomic_data as atomic_data
    charges = [atomic_data.atomic_number[c.capitalize()] for c in atoms_list]
    number_of_atoms = len(atoms_list)
    coords = coordinates
    c_matrix = np.zeros((number_of_atoms, number_of_atoms))
    for i in range(number_of_atoms):
        for j in range(number_of_atoms):
            if i == j:
                c_matrix[i, j] = 0.5 * (charges[i] ** 2.4)
            else:
                r_ij = np.linalg.norm(coords[i, :] - coords[j, :])
                c_matrix[i, j] = (charges[i] * charges[j]) / r_ij
    return c_matrix


def fingerprint(atoms_list, coordinates):
    eigenvalues = np.linalg.eigvals(coulomb_matrix(atoms_list, coordinates))
    eigenvalues[::-1].sort()
    return eigenvalues


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('file', metavar='file',
                        help='Coordinate files')
    parser.add_argument('-d', '--descriptor',
                        choices=['rsm', 'internal'], default='rsm',
                        help='Choose the descriptor')
    args = parser.parse_args()
    file = args.file
    mol = Molecule.from_xyz(file)

    if args.descriptor == 'rsm':
        d = get_rsmd(mol)
        print(d)
        return d

    if args.descriptor == 'internal':
        d = make_internal_coordinates(mol)
        print(d)
        return d


if __name__ == '__main__':
    main()
