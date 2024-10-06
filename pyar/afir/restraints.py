from itertools import product

import autograd.numpy as np
from autograd import grad

import pyar.data.units
from pyar.Molecule import Molecule

# covalent radii (taken from Pyykko and Atsumi, Chem. Eur. J. 15, 2009, 188-197)
# values for metals decreased by 10 %

covalent_radii = {'x ': 0.00,
                  'h': 0.000001, 'he': 0.46, 'li': 1.20, 'be': 0.94, 'b': 0.77,
                  'c': 0.75, 'n': 0.71, 'o': 0.63, 'f': 0.64, 'ne': 0.67,
                  'na': 1.40, 'mg': 1.25, 'al': 1.13, 'si': 1.04, 'p': 1.10,
                  's': 1.02, 'cl': 0.99, 'ar': 0.96, 'k': 1.76, 'ca': 1.54,
                  'sc': 1.33, 'ti': 1.22, 'v': 1.21, 'cr': 1.10, 'mn': 1.07,
                  'fe': 1.04, 'co': 1.00, 'ni': 0.99, 'cu': 1.01, 'zn': 1.09,
                  'ga': 1.12, 'ge': 1.09, 'as': 1.15, 'se': 1.10, 'br': 1.14,
                  'kr': 1.17, 'rb': 1.89, 'sr': 1.67, 'y': 1.47, 'zr': 1.39,
                  'nb': 1.32, 'mo': 1.24, 'tc': 1.15, 'ru': 1.13, 'rh': 1.13,
                  'pd': 1.08, 'ag': 1.15, 'cd': 1.23, 'in': 1.28, 'sn': 1.26,
                  'sb': 1.26, 'te': 1.23, 'i': 1.32, 'xe': 1.31, 'cs': 2.09,
                  'ba': 1.76, 'la': 1.62, 'ce': 1.47, 'pr': 1.58, 'nd': 1.57,
                  'pm': 1.56, 'sm': 1.55, 'eu': 1.51, 'gd': 1.52, 'tb': 1.51,
                  'dy': 1.50, 'ho': 1.49, 'er': 1.49, 'tm': 1.48, 'yb': 1.53,
                  'lu': 1.46, 'hf': 1.37, 'ta': 1.31, 'w': 1.23, 're': 1.18,
                  'os': 1.16, 'ir': 1.11, 'pt': 1.12, 'au': 1.13, 'hg': 1.32,
                  'tl': 1.30, 'pb': 1.30, 'bi': 1.36, 'po': 1.31, 'at': 1.38,
                  'rn': 1.42, 'fr': 2.01, 'ra': 1.81, 'ac': 1.67, 'th': 1.58,
                  'pa': 1.52, 'u': 1.53, 'np': 1.54, 'pu': 1.55}


def get_covalent_radius(z):
    """
    :param z: Atomic Symbol
    :type z: str
    :return: covalent radii
    """
    return pyar.data.units.angstrom2bohr(covalent_radii[z.lower()])

def get_data_structure(atoms, max_cycles):
    coord_size = 3 * len(atoms)
    _1d = (max_cycles,)
    _2d = (max_cycles, coord_size)
    _3d = (max_cycles, coord_size, coord_size)

    get_data_structure = {
        "cart_coords": _2d,
        "energy": _1d,
        "forces": _2d,
        "hessian": _3d,
        "true_energy": _1d,
        "true_forces": _2d,
        "true_hessian": _3d,
    }

    return get_data_structure



def isotropic(fragment_indices, atoms_list, coordinates, force):
    parameter = 6.0  # inverse distance weighting parameter
    epsilon = pyar.data.units.kilojoules2atomic_units(1.0061)
    r_zero = pyar.data.units.angstrom2bohr(3.8164)
    gamma = pyar.data.units.kilojoules2atomic_units(force)
    #    eqn. 3 JCTC 2011,7,2335
    if gamma == 0.0:
        alpha = 0.0
    else:
        alpha = gamma / ((2 ** (-1.0 / 6.0) - (1 + np.sqrt(1 + gamma / epsilon)) ** (-1.0 / 6.0)) * r_zero)

    fragment_one, fragment_two = [coordinates[fragment_list, :] for fragment_list in fragment_indices]
    al1, al2 = [np.array(atoms_list)[fragment_list] for fragment_list in fragment_indices]

    if isinstance(al1, str):
        radius_one = [get_covalent_radius(al1)]
    else:
        radius_one = [get_covalent_radius(z) for z in al1]
    if isinstance(al2, str):
        radius_two = [get_covalent_radius(al2)]
    else:
        radius_two = [get_covalent_radius(z) for z in al2]

    def calculate_restraint_energy(f_one, f_two):
        inter_atomic_distance = np.array([np.linalg.norm(a - b) for a, b in product(f_one, f_two)])
        sum_of_covalent_radii = np.array([(a + b) for a, b in product(radius_one, radius_two)])
        v = np.sum((sum_of_covalent_radii / inter_atomic_distance) ** parameter * inter_atomic_distance)
        w = np.sum((sum_of_covalent_radii / inter_atomic_distance) ** parameter)
        return alpha * v / w

    restraint_energy = calculate_restraint_energy(fragment_one, fragment_two)
    calculate_restraint_gradient = grad(calculate_restraint_energy)
    g_one = calculate_restraint_gradient(fragment_one, fragment_two)
    g_two = calculate_restraint_gradient(fragment_two, fragment_one)
    restraint_gradient = -np.concatenate((g_one, g_two))
    return restraint_energy, restraint_gradient


def main():
    import sys
    from pyar.tabu import merge_two_molecules as merge
    a, b = sys.argv[1:3]
    force = float(sys.argv[3])
    mol_a = Molecule.from_xyz(a)
    mol_b = Molecule.from_xyz(b)
    mol = merge([1, 1, 1, 45, 45, 45], mol_a, mol_b)
    mol.mol_to_xyz('t.xyz')

    e1, g1 = isotropic(mol.fragments, mol.atoms_list, mol.coordinates, force)
    print(e1)
    print(g1)


if __name__ == '__main__':
    main()
