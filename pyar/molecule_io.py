"""Compatibility shim for the legacy ``pyar.molecule_io`` import path.

Use :mod:`pyar.io.xyz` instead.
"""

from pyar.io.xyz import (
    XYZParseError,
    as_coordinates_array,
    parse_xyz,
    read_xyz,
    validate_fragments,
    write_turbomole_coord,
    write_xyz,
)

__all__ = [
    "XYZParseError",
    "as_coordinates_array",
    "parse_xyz",
    "read_xyz",
    "validate_fragments",
    "write_turbomole_coord",
    "write_xyz",
]
