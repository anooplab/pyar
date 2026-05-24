"""PyAR Molecule object.

This module defines :class:`~pyar.Molecule.Molecule`, a lightweight domain
object used throughout PyAR. It is intentionally small and dependency-light.

Compatibility notes
-------------------

PyAR has historically treated XYZ parsing errors as fatal and exited the
process. That behavior is preserved by :func:`read_xyz`. New code should prefer
the exception-based :func:`parse_xyz` helper for improved testability.
"""

import itertools
import sys
from copy import deepcopy
from typing import Optional, Tuple

import numpy as np

import pyar.data.new_atomic_data as atomic_data
import pyar.property
import pyar.molecule_geometry as molecule_geometry
import pyar.molecule_io as molecule_io
import pyar.molecule_merge as molecule_merge
from pyar.molecule_io import XYZParseError


def parse_xyz(filename: str) -> Tuple[list, np.ndarray, str, str, Optional[float]]:
    return molecule_io.parse_xyz(filename)


class Molecule(object):
    """Domain object for a molecular geometry.

    Instances are created from XYZ files or directly from atom lists and
    coordinates. The object keeps PyAR-specific metadata such as charge,
    multiplicity, fragments, and optional electronic energy.
    """

    def __init__(
        self,
        atoms_list,
        coordinates,
        name=None,
        title=None,
        fragments=None,
        charge=0,
        multiplicity=1,
        scftype='rhf',
        energy=None,
    ):
        """Create a molecule from atom symbols and Cartesian coordinates."""
        self.atoms_list = [c.capitalize() for c in atoms_list]
        self._coordinates = None
        self.coordinates = coordinates
        self.charge = charge
        self.multiplicity = multiplicity
        self.scftype = scftype

        # self.elements = [element(z) for z in self.atoms_list]
        # self.atomic_number = [data.new_atomic_data. for n in self.atoms_list]

        # self.covalent_radius = [n.covalent_radius / 100.0 for n in
        #                         self.elements]
        # self.vdw_radius = [n.vdw_radius / 100.0 for n in self.elements]

        self.atomic_number = [atomic_data.atomic_number[z] for z in self.atoms_list]
        self.atomic_mass = [atomic_data.mass[n] for n in self.atoms_list]
        self.covalent_radius = [atomic_data.covalent_radius[n] for n in
                                self.atoms_list]
        self.vdw_radius = [atomic_data.vdw_radius[n] for n in
                           self.atoms_list]

        self.energy = energy

        # self.centre_of_mass = pyar.property.get_centre_of_mass(
        #     self.coordinates, self.atomic_mass)
        # self.average_radius = pyar.property.get_average_radius(
        #     self.coordinates, self.centroid)
        # self.std_of_radius = pyar.property.get_std_of_radius(
        #     self.coordinates, self.centroid)
        # self.distance_list = pyar.property.get_distance_list(
        #     self.coordinates)

        self.name = 'Molecule' if name is None else name
        self.title = 'Title' if title is None else title
        if fragments is None:
            self.fragments = []
        else:
            self.fragments = molecule_io.validate_fragments(fragments, self.number_of_atoms)
            self.fragments_coordinates = self.split_coordinates()
            self.fragments_atoms_list = self.split_atoms_lists()

    @property
    def coordinates(self) -> Optional[np.ndarray]:
        """Cartesian coordinates in Angstrom, or ``None`` after failure."""
        return self._coordinates

    @coordinates.setter
    def coordinates(self, value):
        if value is None:
            self._coordinates = None
            self.number_of_atoms = len(self.atoms_list)
            return
        coordinates = molecule_io.as_coordinates_array(value)
        if coordinates.shape[0] != len(self.atoms_list):
            raise ValueError(
                "atoms_list and coordinates must contain the same number of atoms"
            )
        self._coordinates = coordinates
        self.number_of_atoms = int(self._coordinates.shape[0])

    @property
    def centroid(self) -> Optional[np.ndarray]:
        """Return the current geometry centroid, if coordinates exist."""
        if self._coordinates is None:
            return None
        return pyar.property.get_centroid(self._coordinates)

    def _require_coordinates(self) -> np.ndarray:
        """Return coordinates or raise a clear error when they are missing."""
        if self.coordinates is None:
            raise ValueError("molecule has no coordinates")
        return self.coordinates

    def __str__(self):
        return f"Molecule(name={self.name!r}, atoms={self.number_of_atoms})"

    def __repr__(self):
        return (
            f"Molecule(name={self.name!r}, atoms={self.number_of_atoms}, "
            f"energy={self.energy!r})"
        )

    def __iter__(self):
        return iter(zip(self.atoms_list, self.coordinates))

    def __len__(self):
        return self.number_of_atoms

    def __add__(self, other):
        """Return ``self.merged_with(other)``."""
        return self.merged_with(other)

    def merged_with(self, other):
        """Return a molecule containing this molecule and ``other`` as fragments."""
        return molecule_merge.merged_with(self, other)

    def copy(self):
        """Return a deep copy of this molecule."""
        coordinates = (
            None
            if self.coordinates is None
            else np.array(self.coordinates, dtype=float, copy=True)
        )
        new = Molecule(
            list(self.atoms_list),
            coordinates,
            name=self.name,
            title=self.title,
            fragments=None if coordinates is None else deepcopy(getattr(self, "fragments", [])),
            charge=self.charge,
            multiplicity=self.multiplicity,
            scftype=self.scftype,
            energy=self.energy,
        )
        if coordinates is None:
            new.fragments = deepcopy(getattr(self, "fragments", []))
            for attr in ("fragments_coordinates", "fragments_atoms_list"):
                if hasattr(self, attr):
                    setattr(new, attr, deepcopy(getattr(self, attr)))
        # Preserve optional workflow attributes when present.
        for attr in ("optimized_coordinates",):
            if hasattr(self, attr):
                setattr(new, attr, deepcopy(getattr(self, attr)))
        return new

    @classmethod
    def from_xyz(cls, filename):
        """Create a molecule from an XYZ file."""
        atoms_list, mol_coordinates, mol_name, mol_title, energy = read_xyz(filename)
        return cls(atoms_list, mol_coordinates, name=mol_name, title=mol_title, energy=energy)

    def split_coordinates(self, coordinates=None):
        """Return coordinate blocks for each fragment."""
        if coordinates is None:
            coordinates = self.coordinates
        return [coordinates[fragment_atoms, :] for fragment_atoms in self.fragments]

    def split_atoms_lists(self):
        """Return atom-symbol blocks for each fragment."""
        atoms = np.array(self.atoms_list)
        return [atoms[fragment_atoms].tolist() for fragment_atoms in self.fragments]

    def split_covalent_radii_list(self):
        radii = np.array(self.covalent_radius)
        return [radii[fragment_identifiers] for fragment_identifiers in self.fragments]

    def mol_to_xyz(self, file_name):
        """Write the molecule to an XYZ file."""
        molecule_io.write_xyz(self, file_name)

    def mol_to_turbomole_coord(self):
        """Write the geometry in Turbomole ``coord`` format."""
        molecule_io.write_turbomole_coord(self)

    def is_bonded(self, bond_scale=1.3):
        """Return whether any inter-fragment pair forms a covalent bond.

        ``bond_scale`` follows the connectivity convention used elsewhere in
        PyAR; an unscaled covalent-radius sum is shorter than normal X-H and
        H-H bond lengths and therefore misses real products.
        """
        fragment_one, fragment_two = self.split_coordinates()
        radius_one, radius_two = self.split_covalent_radii_list()
        radii_sum = [
            r1 + r2 for r1, r2 in itertools.product(radius_one, radius_two)
        ]
        distances = [
            np.linalg.norm(a - b)
            for a, b in itertools.product(fragment_one, fragment_two)
        ]
        return any(
            distance < bond_scale * radius
            for distance, radius in zip(distances, radii_sum)
        )

    @property
    def moments_of_inertia_tensor(self):
        """Return the inertia tensor without mutating the molecule."""
        return molecule_geometry.moments_of_inertia_tensor(self)

    def rotate_3d(self, angles):
        """Rotate the molecule in place using Z-X-Z' Euler angles."""
        return molecule_geometry.rotate_3d(self, angles)

    def move_to_origin(self):
        """Shift the molecule so its centroid lies at the origin."""
        return molecule_geometry.move_to_origin(self)

    def move_to_centre_of_mass(self):
        """Shift the molecule so its center of mass lies at the origin."""
        return molecule_geometry.move_to_centre_of_mass(self)

    def translate(self, magnitude):
        """Translate the molecule in place by ``magnitude``."""
        return molecule_geometry.translate(self, magnitude)

    def align(self):
        """Align the molecule with its principal axes in place."""
        return molecule_geometry.align(self)

    # Non-mutating transform helpers (additive API).

    def translated(self, magnitude):
        """Return a translated copy of the molecule."""
        return self.copy().translate(magnitude)

    def rotated_3d(self, angles):
        """Return a rotated copy of the molecule."""
        return self.copy().rotate_3d(angles)

    def aligned(self):
        """Return an aligned copy of the molecule."""
        return self.copy().align()

    def moved_to_origin(self):
        """Return a copy shifted so its centroid lies at the origin."""
        return self.copy().move_to_origin()

    def moved_to_centre_of_mass(self):
        """Return a copy shifted so its center of mass lies at the origin."""
        return self.copy().move_to_centre_of_mass()

    def to_ase_atoms(self):
        """Convert the molecule to an ASE ``Atoms`` object."""
        from ase import Atoms  # optional dependency; imported lazily
        return Atoms(self.atoms_list, positions=self.coordinates)

    @classmethod
    def from_ase_atoms(cls, atoms, name=None, title=None, charge=0, multiplicity=1, scftype="rhf", energy=None):
        """Create a molecule from an ASE ``Atoms`` object."""
        atoms_list = atoms.get_chemical_symbols()
        coords = np.asarray(atoms.get_positions(), dtype=float)
        return cls(
            atoms_list,
            coords,
            name=name,
            title=title,
            fragments=None,
            charge=charge,
            multiplicity=multiplicity,
            scftype=scftype,
            energy=energy,
        )


def read_xyz(filename):
    return molecule_io.read_xyz(filename)


def main():
    xyz_filename = sys.argv[1]
    mol = Molecule.from_xyz(xyz_filename)


if __name__ == '__main__':
    main()
