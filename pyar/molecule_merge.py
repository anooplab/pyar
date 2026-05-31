"""Merge helpers for :mod:`pyar.core.molecule`."""

from __future__ import annotations

import numpy as np

__all__ = ["combine_multiplicity", "merged_with"]


def combine_multiplicity(first, second) -> int:
    """Preserve the existing heuristic for combined multiplicity."""
    total_multiplicity = first + second
    if total_multiplicity == 3:
        return 2
    if total_multiplicity == 4:
        return 1 if first == 2 else 3
    return 1


def merged_with(first, second):
    """Return a molecule containing ``first`` and ``second`` as fragments."""
    from pyar.core.molecule import Molecule

    atoms_list = first.atoms_list + second.atoms_list
    coordinates = np.concatenate((first.coordinates, second.coordinates), axis=0)
    merged = Molecule(atoms_list, coordinates)
    atoms_in_self = list(range(first.number_of_atoms))
    atoms_in_other = list(range(first.number_of_atoms, merged.number_of_atoms))
    merged.fragments = [atoms_in_self, atoms_in_other]
    merged.fragments_coordinates = [first.coordinates, second.coordinates]
    merged.fragments_atoms_list = [first.atoms_list, second.atoms_list]
    merged.name = f"{first.name} + {second.name}"
    merged.title = f"{first.title} + {second.title}"
    merged.charge = first.charge + second.charge
    merged.multiplicity = combine_multiplicity(first.multiplicity, second.multiplicity)
    merged.scftype = "rhf" if first.scftype == "rhf" and second.scftype == "rhf" else "uhf"
    return merged
