"""Compatibility shim for the legacy ``pyar.molecule_geometry`` import path.

Use :mod:`pyar.core.geometry` instead.
"""

from pyar.core.geometry import (
    align,
    moments_of_inertia_tensor,
    move_to_centre_of_mass,
    move_to_origin,
    rotate_3d,
    rotation_matrix,
    translate,
)

__all__ = [
    "align",
    "moments_of_inertia_tensor",
    "move_to_centre_of_mass",
    "move_to_origin",
    "rotate_3d",
    "rotation_matrix",
    "translate",
]
