# -*- coding: utf-8 -*-
"""
    Molecule Class
    --------------

    While running from command-line interface,
    the molecule object is created by reading the
    .xyz file and stores the atom lists and coordinates.
    During the job, new molecule objects are created by
    the aggregator or reactor modules.
"""
import itertools
import logging
import re
import sys
from math import cos, sin

import numpy as np

import pyar.data.new_atomic_data as atomic_data
import pyar.property

molecule_logger = logging.getLogger('pyar.molecule')


class Molecule(object):
    """
    Class used for representing molecule.

    Created either by reading from the xyz file
    or by other modules.

    The following are the main attributes:

    :type number_of_atoms: int
    :param number_of_atoms: The total number of atoms in the molecule
    :type atoms_list: list
    :pa ram atoms_list: The list of Atomic symbols of the atoms.
    :type coordinates: ndarray
    :param coordinates: 2D array of the cartesian coordinates of atoms
        [[x1, y1, z1] [x2, y2, z2]]
        in angstroms units.
    :type name: str
    :param name: The name of the molecule
    :type title: str
    :param title: Usually the second line in the xyz file.
    :type fragments: list
    :pa ram fragments: The list of atoms in each fragment required
        for the Reaction module.

    The following attributes are read from the
    data(atomic_data.py) and stored as a list in the
    order as the atoms list.

    param atomic_number: list
        The list of atomic numbers
    atomic_mass: list
        Atomic masses. Required for calculating the
        centre of mass.
    covalent_radius: list
        The covalent radii of atoms.  Required for
        calculating the close contact of banded-status.
    vdw_radius: list
        The van der Waals radii of atoms. Required for
        the placing the second fragment near the first.

    The following attributes are calculated using the
    above data.

    centroid: ndarray
        The centroid of the molecule (x,y,z)
    centre_of_mass: ndarray
        The centre of mass (x, y, z)
    average_radius: float
        The average radius of all the distances of
        atoms= from the centroid.
    std_of_radius: float
        The standard deviation of all the distances
        of atoms from the centroid.
    distance_list: ndarray
        The 2D array of interatomic distances.
        Useful for calculating fingerprints, coulomb
        matrix etc.

    Following attributes are electronic properties.  This will be
    available after a calculation

    energy: float
        The energy of the molecule. This is added after
        any quantum chemical calculations.

    """

    def __init__(self, atoms_list, coordinates, name=None,
                 title=None, fragments=None, charge=0,
                 multiplicity=1, scftype='rhf', energy=None):
        """
        Init function for Molecule

        :type multiplicity: int
        :type atoms_list: list
        :param atoms_list: The list of atomic symbols
        :type coordinates: ndarray
        :param coordinates: atomic coordinates
        :type name: str
        :param name: The name of the molecule
        :type title: str
        :param title: Usually the second line in the xyz file.
        :type fragments: list
        :param fragments: The list of atoms in each fragment required for the Reaction module.

        """

        self.number_of_atoms = len(coordinates)
        self.atoms_list = [c.capitalize() for c in atoms_list]
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

        self.centroid = pyar.property.get_centroid(self.coordinates)
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
            self.fragments = fragments
            self.fragments_coordinates = self.split_coordinates()
            self.fragments_atoms_list = self.split_atoms_lists()

    def __str__(self):
        return f"Name: {self.name}\n Coordinates:{self.coordinates}"

    def __repr__(self):
        return f"Molecule.from_xyz('{self.name}.xyz')"

    def __iter__(self):
        pass

    def __len__(self):
        return self.number_of_atoms

    def __add__(self, other):
        """
        Merge the 'other' molecule with 'self'

        Merges two molecules objects.

        :type other: Molecule
        :param other: Molecule object to be merged with self.
        :return: Merged Molecule object
        :rtype: object

        """
        atoms_list = self.atoms_list + other.atoms_list
        coordinates = np.concatenate((self.coordinates, other.coordinates),
                                     axis=0)
        merged = Molecule(atoms_list, coordinates)
        atoms_in_self = list(range(self.number_of_atoms))
        atoms_in_other = list(
            range(self.number_of_atoms, merged.number_of_atoms))
        merged.fragments = [atoms_in_self, atoms_in_other]
        merged.fragments_coordinates = [self.coordinates, other.coordinates]
        merged.fragments_atoms_list = [self.atoms_list, other.atoms_list]
        if hasattr(self, 'fragments_history'):
            merged.fragments_history = self.fragments_history + atoms_in_other
        else:
            merged.fragments_history = [atoms_in_self, atoms_in_other]
        merged.name = f"{self.name} + {other.name}"
        merged.title = f"{self.title} + {other.title}"

        merged.charge = self.charge + other.charge
        total_multiplicity = self.multiplicity + other.multiplicity
        if total_multiplicity == 3:
            new_multiplicity = 2
        elif total_multiplicity == 4:
            new_multiplicity = 1 if self.multiplicity == 2 else 3
        else:
            new_multiplicity = 1  # Complicated case.
        merged.multiplicity = new_multiplicity

        if self.scftype == 'rhf' and other.scftype == 'rhf':
            combined_scftype = 'rhf'
        else:
            combined_scftype = 'uhf'
        merged.scftype = combined_scftype
        return merged

    @classmethod
    def from_xyz(cls, filename):
        """
        Instantiates Molecule object from .xyz file

        Reads .xyz files, extracts list of atoms
        and coordinates. The name is set as the name
        of the molecule. The comment line of .xyz file
        is stored as title.

        :param filename: name of .xyz file
        :type filename: str
        :return: Molecule object
        :rtype: object

        """

        atoms_list, mol_coordinates, mol_name, mol_title, energy = read_xyz(
            filename)
        return cls(atoms_list, mol_coordinates, name=mol_name,
                   title=mol_title, energy=energy)

    def split_coordinates(self, coordinates=None):
        """
        This function splits coordinates into a list of coordinates
        for each fragment in self.fragments. 
        
        Inputs:
            coordinates (optional): Coordinates to be split. 
                                    Defaults to self.coordinates.
        
        Outputs:
            split_coordinates: List of coordinates for each fragment.
        
        Steps:
            1. If no coordinates are specified, use self.coordinates 
               as the coordinates to split.
            2. Iterate through each fragment in self.fragments.
            3. For each fragment, extract the coordinates for the 
               corresponding atoms.
            4. Append the extracted coordinates to the list of 
               split_coordinates. 
            5. Return the list of split_coordinates.

        :type coordinates: ndarray
        :param coordinates: coordinates

        """

        if coordinates is None:
            coordinates = self.coordinates
        return [coordinates[fragment_atoms, :] for fragment_atoms in
                self.fragments]

    def split_atoms_lists(self):
        """
        Split the list of atoms in to different fragment.

        :return: A list of list aof atoms in fragments.
        :rtype: list

        """

        return [self.atoms_list[fragment_atoms, :] for fragment_atoms in
                self.fragments]

    def split_covalent_radii_list(self):
        return [np.array(self.covalent_radius)[fragment_identifiers] for
                fragment_identifiers in self.fragments]

    def mol_to_xyz(self, file_name):
        """
        Write an xyz file of the Molecule.

        :param file_name: Output .xyz file
        :type file_name: str

        """
        with open(file_name, 'w') as fp:
            fp.write("{:3d}\n".format(self.number_of_atoms))
            fp.write(f"{self.title}: {self.energy}\n")
            for element_symbol, atom_coordinate in zip(self.atoms_list, self.coordinates):
                fp.write(("%-2s%12.5f%12.5f%12.5f\n" % (element_symbol, atom_coordinate[0], atom_coordinate[1], atom_coordinate[2])))

    def mol_to_turbomole_coord(self):
        """
        Write 'coord' file in Turbomole format

        """
        with open('coord', 'w') as fp:
            fp.write("$coord\n")
            coords = self.coordinates
            atoms_list = self.atoms_list
            for i in range(self.number_of_atoms):
                fp.write("%20.14f  %20.14f  %20.14f  %6s\n" % (
                    coords[i][0], coords[i][1], coords[i][2],
                    atoms_list[i].lower()))
            fp.write("$end\n")

    def is_bonded(self):
        """
        Checks if there is any bond between two fragments in the molecule.

        This is a simple distance check.  If any of the distance between
        any two atoms of different fragments
        distance is smaller than the sum of covalent radii, it is
        considered as a bond

        :return: boolean
        """
        fragment_one, fragment_two = self.split_coordinates()
        radius_one, radius_two = self.split_covalent_radii_list()
        if isinstance(radius_one, np.float) and isinstance(radius_two,
                                                           np.float64):
            # noinspection PyPep8Naming
            R = radius_one + radius_two
            r = np.linalg.norm(fragment_one - fragment_two)
            return r < R
        else:
            # noinspection PyPep8Naming
            R = [x + y for x, y in itertools.product(radius_one, radius_two)]
            r = [np.linalg.norm(a - b) for a, b in
                 itertools.product(fragment_one, fragment_two)]
            return any(a < b for a, b in zip(r, R))

    @property
    def moments_of_inertia_tensor(self):
        """
        Calculate moments of inertia

        :return: ndarray
        """
        mass = self.atomic_mass
        self.move_to_centre_of_mass()
        x = self.coordinates[:, 0]
        y = self.coordinates[:, 1]
        z = self.coordinates[:, 2]
        i_xx = np.sum(mass * (y * y + z * z))
        i_yy = np.sum(mass * (x * x + z * z))
        i_zz = np.sum(mass * (x * x + y * y))
        i_xy = -np.sum(mass * (x * y))
        i_yz = -np.sum(mass * (y * z))
        i_xz = -np.sum(mass * (x * z))
        return np.array([[i_xx, i_xy, i_xz],
                         [i_xy, i_yy, i_yz],
                         [i_xz, i_yz, i_zz]]) / self.number_of_atoms

    def rotate_3d(self, angles):
        """
        This function will rotate a molecule by Euler rotation theorem.
        The Z-X-Z' rotation convention is followed. The range of phi, theta
        and psi should be (0,360),(0,180) and (0,360) respectively.
        This function will first translate the molecule to its center
        of mass(centroid). Then, it rotate the molecule and translate
        to its original position. So, to rotate a molecule around the origin,
        (0.,0.,0.), set_origin usage is necessary
        """
        phi, theta, psi = angles
        matrix_d = np.array(((cos(phi), sin(phi), 0.),
                             (-sin(phi), cos(phi), 0.),
                             (0., 0., 1.)))
        matrix_c = np.array(((1., 0., 0.),
                             (0., cos(theta), sin(theta)),
                             (0., -sin(theta), cos(theta))))
        matrix_b = np.array(((cos(psi), sin(psi), 0.),
                             (-sin(psi), cos(psi), 0.),
                             (0., 0., 1.)))
        matrix_a = np.dot(matrix_b, np.dot(matrix_c, matrix_d))
        new_coordinates = np.dot(matrix_a, np.transpose(self.coordinates))
        self.coordinates = np.transpose(new_coordinates) + self.centroid
        return self

    def move_to_origin(self):
        self.translate(pyar.property.get_centroid(self.coordinates))
        return self

    def move_to_centre_of_mass(self):
        self.translate(pyar.property.get_centre_of_mass(self.coordinates,
                                                        self.atomic_mass))
        return self

    def translate(self, magnitude):
        self.coordinates -= magnitude
        return self

    def align(self):
        """
        Align the molecules to the principal axis

        :return: aligned coordinates
        """
        moi = self.moments_of_inertia_tensor
        eigenvalues, eigen_vectors = np.linalg.eig(moi)
        transformed_coordinates = np.array(
            np.dot(self.coordinates, eigen_vectors))
        order = [0, 1, 2]
        for p in range(3):
            for q in range(p + 1, 3):
                if moi.item(p, p) < moi.item(q, q):
                    temp = order[p]
                    order[p] = order[q]
                    order[q] = temp
        move_axes = np.zeros((3, 3))
        for p in range(3):
            move_axes[p][order[p]] = 1.0
        self.coordinates = np.dot(transformed_coordinates, move_axes)
        return self


def read_xyz(filename):
    with open(filename) as fp:
        f = fp.readlines()
    try:
        number_of_atoms = int(f[0])
    except Exception as e:
        molecule_logger.error(e)
        molecule_logger.error(f"{filename} should have number of atoms in the first line")

        molecule_logger.error("but we found\n %s" % f[0])
        molecule_logger.error("Is it an xyz file?")
        sys.exit(f'Error in reading {filename}')
    mol_title = f[1].rstrip()
    try:
        energy = float(re.split(':|=|\s+', mol_title)[1])
    except Exception as e:
        molecule_logger.debug(f"No energy found\n{e}")
        energy = None
    try:
        geometry_section = [each_line.split() for each_line in f[2:] if len(each_line) >= 4]

    except Exception as e:
        molecule_logger.error(e)
        molecule_logger.error("Something wrong with reading the geometry section")
        sys.exit(f'Error in reading {filename}')
    if len(geometry_section) != number_of_atoms:
        molecule_logger.error("Number of geometric coordinates is not equal to number of atoms")

        molecule_logger.error("Is something wrong?")
        sys.exit(f'Error in reading {filename}')
    atoms_list = []
    coordinates = []
    for i, c in enumerate(geometry_section):
        try:
            symbol = c[0].capitalize()
            x_coord = float(c[1])
            y_coord = float(c[2])
            z_coord = float(c[3])
        except Exception as e:
            molecule_logger.error(e)
            molecule_logger.error("Something wrong in line: %d" % (i + 1))
            molecule_logger.error(c)
            sys.exit(f'Error in reading {filename}')
        atoms_list.append(symbol)
        coordinates.append([x_coord, y_coord, z_coord])
    mol_coordinates = np.array(coordinates)
    mol_name = filename[:-4]
    return atoms_list, mol_coordinates, mol_name, mol_title, energy


def main():
    xyz_filename = sys.argv[1]
    mol = Molecule.from_xyz(xyz_filename)


if __name__ == '__main__':
    main()
