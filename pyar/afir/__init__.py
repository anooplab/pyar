"""Compatibility package for legacy ``pyar.afir`` imports."""

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
