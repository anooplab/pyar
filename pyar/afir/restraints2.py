from itertools import product

import numpy as np

from pyar.data.units import kilojoules2atomic_units, angstrom2bohr


def isotropic(mol, force):
    assert len(mol.fragments) == 2

    parameter = 6.0  # inverse distance weighting parameter

    fragment1, fragment2 = mol.fragments_coordinates
    al1, al2 = mol.fragments_atoms_list
    radius_one, radius_two = [[mol.covalent_radius[r] for r in f] for f in mol.fragments]

    inter_atomic_distance_ij = np.array([np.linalg.norm(a - b) for a, b in product(fragment1, fragment2)])
    sum_of_covalent_radii__ij = np.array([(a + b) for a, b in product(radius_one, radius_two)])

    nom = np.sum((sum_of_covalent_radii__ij / inter_atomic_distance_ij) ** parameter * inter_atomic_distance_ij)
    den = np.sum((sum_of_covalent_radii__ij / inter_atomic_distance_ij) ** parameter)

    alpha = calculate_alpha(force)

    restraint_energy = alpha * nom / den

    gradients_one = np.zeros((len(fragment1), 3))
    for i, rc1 in enumerate(zip(radius_one, fragment1)):
        Ri, ci = rc1
        gradient_i = np.zeros(3)
        for j, rc2 in enumerate(zip(radius_two, fragment2)):
            Rj, cj = rc2
            sum_of_covalent_radii_ni = Ri + Rj
            inter_atomic_distance_ni = np.linalg.norm(cj - ci)
            if inter_atomic_distance_ni == 0.0:
                gradient_i += 0.0
            else:
                dq = cj - ci
                first_term = ((1 - parameter) / den * sum_of_covalent_radii_ni ** parameter *
                              inter_atomic_distance_ni ** (-1 - parameter))
                second_term = (parameter * nom / den ** 2 * sum_of_covalent_radii_ni ** parameter *
                               inter_atomic_distance_ni ** (-2 - parameter))
                this_gradient = first_term * dq + second_term * dq
                gradient_i += this_gradient
        gradients_one[i] = gradient_i

    gradients_two = np.zeros((len(fragment2), 3))
    for i, rc1 in enumerate(zip(radius_two, fragment2)):
        Ri, ci = rc1
        gradient_i = np.zeros(3)
        for j, rc2 in enumerate(zip(radius_one, fragment1)):
            Rj, cj = rc2
            sum_of_covalent_radii_ni = Ri + Rj
            inter_atomic_distance_ni = np.linalg.norm(cj - ci)
            if inter_atomic_distance_ni == 0.0:
                gradient_i += 0.0
            else:
                dq = cj - ci
                first_term = ((1 - parameter) / den * sum_of_covalent_radii_ni ** parameter *
                              inter_atomic_distance_ni ** (-1 - parameter))
                second_term = (parameter * nom / den ** 2 * sum_of_covalent_radii_ni ** parameter *
                               inter_atomic_distance_ni ** (-2 - parameter))
                this_gradient = first_term * dq + second_term * dq
                gradient_i += this_gradient
        gradients_two[i] = gradient_i
    gradients = np.concatenate((gradients_one, gradients_two), axis=0)

    from pprint import pprint
    pprint(gradients)
    return restraint_energy, gradients


def calculate_alpha(force):
    gamma = kilojoules2atomic_units(force)
    epsilon = kilojoules2atomic_units(1.0061)
    r_zero = angstrom2bohr(3.8164)
    #    eqn. 3 JCTC 2011,7,2335
    alpha = gamma / ((2 ** (-1.0 / 6.0) - (1 + np.sqrt(1 + gamma / epsilon)) ** (-1.0 / 6.0)) * r_zero)
    return alpha


def main():
    import sys
    from pyar.Molecule import Molecule
    from pyar.tabu import merge_two_molecules as merge

    a, b = sys.argv[1:3]
    force = float(sys.argv[3])
    mol_a = Molecule.from_xyz(a)
    mol_b = Molecule.from_xyz(b)
    mol = merge([1, 0, 0, 0, 0, 0], mol_a, mol_b)

    e1, g1 = isotropic(mol, force)
    print(e1, g1)


if __name__ == '__main__':
    main()
