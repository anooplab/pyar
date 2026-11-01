import logging
import math

import numpy as np

import pyar.data.units
from pyar.core.molecule import Molecule

afir_logger = logging.getLogger("pyar.biases.afir")

# covalent radii (taken from Pyykko and Atsumi, Chem. Eur. J. 15, 2009, 188-197)
# values for metals decreased by 10 %

from pyar.constants import AFIR_COVALENT_RADII_ANGSTROM as covalent_radii

__all__ = [
    "get_covalent_radius",
    "get_data_structure",
    "isotropic",
    "alpha_from_gamma",
    "resolve_gamma",
]


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


def resolve_gamma(value, fallback=100.0):
    """Return a finite numeric AFIR gamma value."""
    if value is None:
        return float(fallback)
    try:
        gamma = float(value)
    except (TypeError, ValueError):
        raise ValueError(f"Invalid AFIR gamma value: {value!r}") from None
    if not math.isfinite(gamma):
        raise ValueError(f"Invalid AFIR gamma value: {value!r}")
    return gamma


def alpha_from_gamma(gamma_kjmol):
    """Map the AFIR collision parameter in kJ mol-1 to alpha in Ha Bohr-1."""
    gamma = pyar.data.units.kilojoules2atomic_units(gamma_kjmol)
    if gamma == 0.0:
        return 0.0
    epsilon = pyar.data.units.kilojoules2atomic_units(1.0061)
    r_zero = pyar.data.units.angstrom2bohr(3.8164)
    return gamma / ((2 ** (-1.0 / 6.0) - (1 + np.sqrt(1 + gamma / epsilon)) ** (-1.0 / 6.0)) * r_zero)


def isotropic(fragment_indices, atoms_list, coordinates, force):
    parameter = 6.0  # inverse distance weighting parameter
    alpha = alpha_from_gamma(force)

    from pyar.biases.collective_coordinates import evaluate_contact_coordinate

    q, dq_dR, _ = evaluate_contact_coordinate(
        fragment_indices, atoms_list, coordinates, kind="afir", distance_power=parameter
    )
    return alpha * q, -alpha * np.asarray(dq_dR, dtype=float)


def main():
    import sys
    from pyar.sampling.trial_generator import merge_two_molecules as merge
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
