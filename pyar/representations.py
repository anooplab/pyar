import itertools
from itertools import product

import numpy as np

import pyar.property
from pyar.Molecule import Molecule
# import torch
# import torchani
# from DBCV import DBCV
from dscribe.descriptors import MBTR, SOAP, LMBTR, ACSF, SineMatrix, ValleOganov
# from ase.io import read
from ase import Atoms
# import glob

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
    for a_i, a_j in product(range(len(coordinates)), range(len(coordinates))):
        if bm[a_i, a_j]:
            seq = [a_i, a_j]
            if seq[::-1] not in [i[0] for i in bl]:
                bl.append([seq, dm[a_i, a_j]])
    al = []
    for a_i, a_j, k in product(range(len(coordinates)), repeat=3):
        if a_i != a_j and a_j != k and a_i != k and bm[a_i, a_j] and bm[a_j, k]:
            seq = [a_i, a_j, k]
            if seq[::-1] not in [i[0] for i in al]:
                al.append([seq, pyar.property.calculate_angle(coordinates[a_i], coordinates[a_j], coordinates[k])])

    dl = []
    for a_i, a_j, k, l in product(range(len(coordinates)), repeat=4):
        if a_i != a_j and a_i != k and a_i != l and a_j != k and a_j != l and k != l and bm[a_i, a_j] and bm[a_j, k] and bm[k, l]:
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
    # charges = [atomic_data.atomic_number[c.capitalize()] for c in atoms_list]
    charges = [atomic_data.atomic_number[c.capitalize()] if isinstance(c, str) else c for c in atoms_list]
    number_of_atoms = len(atoms_list)
    coords = coordinates
    c_matrix = np.zeros((number_of_atoms, number_of_atoms))
    for i, j in product(range(number_of_atoms), range(number_of_atoms)):
        if i == j:
            c_matrix[i, j] = 0.5 * charges[i]**2.4
        else:
            r_ij = np.linalg.norm(coords[i, :] - coords[j, :])
            c_matrix[i, j] = charges[i] * charges[j] / r_ij
    return c_matrix

#last old version
# def coulomb_matrix(atoms_list, coordinates):
#     """
#     :return: Coulomb Matrix
#     """
#     from pyar.data import new_atomic_data as atomic_data

#     charges = [atomic_data.atomic_number[c.capitalize()] if isinstance(c, str) else c for c in atoms_list]
#     number_of_atoms = len(atoms_list)
#     coords = coordinates
#     c_matrix = np.zeros((number_of_atoms, number_of_atoms))

#     for i, j in product(range(number_of_atoms), range(number_of_atoms)):
#         charge_i = charges[i]
#         charge_j = charges[j]

#         if isinstance(charge_i, np.ndarray):
#             charge_i = charge_i.item()
#         if isinstance(charge_j, np.ndarray):
#             charge_j = charge_j.item()

#         if i == j:
#             c_matrix[i, j] = 0.5 * charge_i**2.4
#         else:
#             r_ij = np.linalg.norm(coords[i, :] - coords[j, :])
#             c_matrix[i, j] = charge_i * charge_j / r_ij

#     return c_matrix

# def coulomb_matrix(atoms_list, coordinates):
#     """
#     :return: Coulomb Matrix
#     """
#     from pyar.data import new_atomic_data as atomic_data

#     charges = [atomic_data.atomic_number[c.capitalize()] if isinstance(c, str) else c for c in atoms_list]
#     number_of_atoms = len(atoms_list)
#     coords = coordinates
#     c_matrix = np.zeros((number_of_atoms, number_of_atoms))

#     for i, j in product(range(number_of_atoms), range(number_of_atoms)):
#         charge_i = charges[i]
#         charge_j = charges[j]

#         if isinstance(charge_i, np.ndarray):
#             charge_i = charge_i[0]  # Take the first element of the array
#         if isinstance(charge_j, np.ndarray):
#             charge_j = charge_j[0]  # Take the first element of the array

#         if i == j:
#             c_matrix[i, j] = 0.5 * charge_i**2.4
#         else:
#             r_ij = np.linalg.norm(coords[i, :] - coords[j, :])
#             if r_ij < 1e-8:  # Check if the distance is very close to zero
#                 c_matrix[i, j] = 0.0  # Set the value to zero to avoid division by near-zero
#             else:
#                 c_matrix[i, j] = charge_i * charge_j / r_ij

#     # Check if the Coulomb matrix contains any NaN or infinite values
#     if np.isnan(c_matrix).any() or np.isinf(c_matrix).any():
#         raise ValueError("Coulomb matrix contains NaN or infinite values")

#     return c_matrix



def fingerprint(atoms_list, coordinates):
    eigenvalues = np.linalg.eigvals(coulomb_matrix(atoms_list, coordinates))
    # eigenvalues[::-1].sort()
    eigenvalues = np.sort(eigenvalues)[::-1]
    return eigenvalues

def cutoff_func(r, Rc):
    return 0.5 * (np.cos(np.pi * r / Rc) + 1) if r < Rc else 0

def radial_symmetry_func(r, eta, Rs, Rc):
    return np.exp(-eta * (r - Rs)**2) * cutoff_func(r, Rc)

def angular_symmetry_func(r_ij, r_ik, r_jk, theta_ijk, eta, zeta, lambda_, Rc):
    return 2**(1 - zeta) * (1 + lambda_ * np.cos(theta_ijk))**zeta * np.exp(-eta * (r_ij**2 + r_ik**2 + r_jk**2)) * cutoff_func(r_ij, Rc) * cutoff_func(r_ik, Rc) * cutoff_func(r_jk, Rc)

def generate_aev(atom_list, coordinates):
    # Define the Behler-Parrinello constants
    eta = 0.5
    Rs = [1.0, 2.0, 3.0]  # Example radial symmetry function distances
    zeta = 1.0
    lambda_ = 1.0
    Rc = 6.0
    aev = []
    for i, atom_i in enumerate(atom_list):
        G2 = []
        G3 = []
        for j, atom_j in enumerate(atom_list):
            if i != j:
                r_ij = np.linalg.norm(coordinates[i] - coordinates[j])
                G2.extend(radial_symmetry_func(r_ij, eta, Rs, Rc))
                for k, atom_k in enumerate(atom_list):
                    if i != k and j < k:
                        r_ik = np.linalg.norm(coordinates[i] - coordinates[k])
                        r_jk = np.linalg.norm(coordinates[j] - coordinates[k])
                        theta_ijk = np.arccos(np.dot(coordinates[i] - coordinates[j], coordinates[i] - coordinates[k]) / (r_ij * r_ik))
                        G3.append(angular_symmetry_func(r_ij, r_ik, r_jk, theta_ijk, eta, zeta, lambda_, Rc))
        aev.append(np.concatenate((G2, G3)))
    return aev

# LMBTR Descriptor
def lmbtr_descriptor(atoms_list, coordinates):
    # Create an ASE Atoms object from the atoms_list and coordinates
    molecule = Atoms(atoms_list, positions=coordinates)

    # Get unique species from atoms_list
    unique_species = list(set(atoms_list))

    #setup LMBTR
    lmbtr = LMBTR(
        species=unique_species,
        geometry={"function": "distance"},
        grid={"min": 0, "max": 5, "n": 100, "sigma": 0.1},
        weighting={"function": "exp", "scale": 0.5, "threshold": 1e-3},
        periodic=False,
        normalization="l2",
    )

    # Create LMBTR output for the molecule
    lmbtr_output = lmbtr.create(molecule)

    return lmbtr_output


def acsf_descriptor(atoms_list, coordinates):
    # Create an Atoms object from the atoms_list and coordinates
    molecule = Atoms(atoms_list, positions=coordinates)
    unique_species = list(set(atoms_list))
    #setup ACSF
    acsf = ACSF(
        species=unique_species,
        r_cut=6.0,
        g2_params=[[1, 1], [1, 2], [1, 3]],
        g4_params=[[1, 1, 1], [1, 2, 1], [1, 1, -1], [1, 2, -1]],
    )

    # Create ACSF output for the molecule
    acsf_output = acsf.create(molecule)

    return acsf_output



def sinematrix_descriptor(atoms_list, coordinates):
    # Create an Atoms object from the atoms_list and coordinates
    molecule = Atoms(atoms_list, positions=coordinates)
    

    # Setting up the sine matrix descriptor
    sm = SineMatrix(
        n_atoms_max=len(atoms_list),
        permutation="sorted_l2",
        sparse=False
    )

    # Create SineMatrix output for the molecule
    sm_output = sm.create(molecule)

    return sm_output



def valleoganov_descriptor(atoms_list, coordinates):    
    # Create an Atoms object from the atoms_list and coordinates
    molecule = Atoms(atoms_list, positions=coordinates)
    unique_species = list(set(atoms_list))

    # Setup
    vo = ValleOganov(
        species=unique_species,
        function="distance",
        sigma=10**(-0.5),
        n=100,
        r_cut=5
    )

    # Create ValleOganov output for the molecule
    vo_output = vo.create(molecule)

    return vo_output



def mbtr_descriptor(atoms_list, coordinates):
    # Create an Atoms object from the atoms_list and coordinates
    molecule = Atoms(atoms_list, positions=coordinates)

    # Get unique species from atoms_list
    unique_species = list(set(atoms_list))

    # Setup MBTR
    mbtr = MBTR(
        species=unique_species,
        geometry={"function": "inverse_distance"},
        grid={"min": 0, "max": 1, "n": 100, "sigma": 0.1},
        weighting={"function": "exp", "scale": 0.5, "threshold": 1e-3},
        periodic=False,
        normalization="l2",
)

    # Create LMBTR output for the molecule
    mbtr_output = mbtr.create(molecule)

    return mbtr_output






# def mbtr_descriptor(atoms_list, coordinates):
#     # Create an Atoms object from the atoms_list and coordinates
#     molecule = Atoms(atoms_list, positions=coordinates)

#     # Get unique species from atoms_list
#     unique_species = list(set(atoms_list))
    
#     k2min, k2max, k2n = 0.7, 2.0, 100
    
#     mbtr = MBTR(
#         species=unique_species,
#         geometry={"function": "distance"},
#         grid={"min": k2min, "max": k2max, "n": k2n, "sigma": 0.000000001},
#         weighting={"function": "exp", "scale": 0.5, "threshold": 3e-3},
#         periodic=False,
#         normalization="l2_each",
#     )

#     # Create MBTR output for the molecule
#     mbtr_output = mbtr.create(molecule)

#     return mbtr_output


def soap_descriptor(atoms_list, coordinates):
    # Create an Atoms object from the atoms_list and coordinates
    molecule = Atoms(atoms_list, positions=coordinates)
    # Get unique species from atoms_list
    unique_species = list(set(atoms_list))
    # Setup SOAP
    soap = SOAP(
        species=unique_species,
        periodic=False,
        r_cut=5,
        n_max=8,
        l_max=8,
        average="off",
        sparse=True
    )
    # Create SOAP output for the molecule
    soap_output = soap.create(molecule)
    return soap_output

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
