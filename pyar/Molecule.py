"""Compatibility shim for the legacy ``pyar.Molecule`` import path.

Use :mod:`pyar.core.molecule` instead.
"""

from pyar.core.molecule import Molecule, XYZParseError, parse_xyz, read_xyz

__all__ = [
    "Molecule",
    "XYZParseError",
    "parse_xyz",
    "read_xyz",
]
