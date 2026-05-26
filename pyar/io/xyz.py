"""I/O helpers for :mod:`pyar.core.molecule`."""

from __future__ import annotations

import logging
import re
import sys
from pathlib import Path
from typing import Optional, Tuple

import numpy as np

molecule_logger = logging.getLogger("pyar.molecule")

__all__ = [
    "as_coordinates_array",
    "validate_fragments",
    "XYZParseError",
    "parse_xyz",
    "read_xyz",
    "write_xyz",
    "write_turbomole_coord",
]


class XYZParseError(ValueError):
    """Raised when an XYZ file cannot be parsed."""


def as_coordinates_array(value) -> np.ndarray:
    """Validate and return Cartesian coordinates as an ``(N, 3)`` array."""
    arr = np.asarray(value, dtype=float)
    if arr.ndim != 2 or arr.shape[1] != 3:
        raise ValueError(f"coordinates must have shape (N, 3); got {arr.shape!r}")
    return arr


def validate_fragments(fragments, number_of_atoms: int) -> list:
    """Validate atom-index fragment definitions for a molecule."""
    checked = []
    for fragment in fragments:
        indices = list(fragment)
        if any(index < 0 or index >= number_of_atoms for index in indices):
            raise ValueError("fragment atom index is outside the molecule")
        checked.append(indices)
    return checked


def parse_xyz(filename: str) -> Tuple[list, np.ndarray, str, str, Optional[float]]:
    """Parse an XYZ file and raise :class:`XYZParseError` on malformed input."""
    path = Path(filename)
    try:
        text = path.read_text().splitlines()
    except OSError as exc:
        raise XYZParseError(f"Could not read XYZ file {filename!r}: {exc}") from exc

    if not text:
        raise XYZParseError(f"XYZ file {filename!r} is empty")

    try:
        number_of_atoms = int(text[0].strip())
    except Exception as exc:
        raise XYZParseError(
            f"{filename!r} should have number of atoms in the first line; found {text[0]!r}"
        ) from exc

    if number_of_atoms < 1:
        raise XYZParseError(f"{filename!r} has invalid atom count {number_of_atoms}")

    if len(text) < 2:
        raise XYZParseError(f"{filename!r} is missing the title/comment line")

    mol_title = text[1].rstrip("\n")
    energy_match = re.search(
        r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?\s*$",
        mol_title,
    )
    energy = float(energy_match.group()) if energy_match else None

    geometry_lines = []
    for line in text[2:]:
        fields = line.split()
        if len(fields) >= 4:
            geometry_lines.append(fields)
    if len(geometry_lines) != number_of_atoms:
        raise XYZParseError(
            f"{filename!r} has {len(geometry_lines)} coordinate lines but declares {number_of_atoms} atoms"
        )

    atoms_list = []
    coords = []
    for i, fields in enumerate(geometry_lines, start=1):
        try:
            symbol = fields[0].capitalize()
            x_coord = float(fields[1])
            y_coord = float(fields[2])
            z_coord = float(fields[3])
        except Exception as exc:
            raise XYZParseError(f"Bad XYZ geometry line {i} in {filename!r}: {fields!r}") from exc
        atoms_list.append(symbol)
        coords.append([x_coord, y_coord, z_coord])

    mol_coordinates = as_coordinates_array(coords)
    mol_name = str(path)[:-4] if str(path).lower().endswith(".xyz") else str(path)
    return atoms_list, mol_coordinates, mol_name, mol_title, energy


def read_xyz(filename: str):
    """Read an XYZ file using PyAR's historical exit-on-error interface."""
    try:
        return parse_xyz(filename)
    except XYZParseError as exc:
        molecule_logger.error(str(exc))
        sys.exit(f"Error in reading {filename}")


def write_xyz(molecule, file_name):
    """Write a molecule and optional energy label in XYZ format."""
    title = getattr(molecule, "title", getattr(molecule, "name", "Molecule"))
    energy = getattr(molecule, "energy", 0.0)
    with open(file_name, "w") as fp:
        fp.write("{:3d}\n".format(molecule.number_of_atoms))
        fp.write(f"{title}: {energy}\n")
        for element_symbol, atom_coordinate in zip(molecule.atoms_list, molecule.coordinates):
            fp.write(
                "%-2s%12.5f%12.5f%12.5f\n"
                % (element_symbol, atom_coordinate[0], atom_coordinate[1], atom_coordinate[2])
            )


def write_turbomole_coord(molecule, file_name="coord"):
    """Write a molecule in Turbomole ``coord`` format."""
    with open(file_name, "w") as fp:
        fp.write("$coord\n")
        coords = molecule.coordinates
        atoms_list = molecule.atoms_list
        for i in range(molecule.number_of_atoms):
            fp.write(
                "%20.14f  %20.14f  %20.14f  %6s\n"
                % (coords[i][0], coords[i][1], coords[i][2], atoms_list[i].lower())
            )
        fp.write("$end\n")
