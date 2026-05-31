"""Compatibility shim for the legacy ``pyar.afir.restraints`` path.

Use :mod:`pyar.biases.afir` instead.
"""

from pyar.biases.afir import (
    covalent_radii,
    get_covalent_radius,
    get_data_structure,
    isotropic,
    main,
    resolve_gamma,
)

__all__ = [
    "covalent_radii",
    "get_covalent_radius",
    "get_data_structure",
    "isotropic",
    "main",
    "resolve_gamma",
]
