from math import radians, cos, sin, pi
from scipy.spatial import distance as scipy_distance

import numpy as np
from itertools import combinations, product

from atomic_data import atomic_masses, atomic_numbers, covalent_radii, vdw_radii
import itertools
import logging
molecule_logger = logging.getLogger('pyar.molecule')


class Molecule(object):
    def __init__(self, atoms_list, coordinates, name=None, title=None, fragments=None):
        self.number_of_atoms = len(coordinates)
        self.atoms_list = [c.capitalize() for c in atoms_list]
        self.atomic_number = self.get_atomic_number()
        self.atomic_mass = self.get_atomic_mass()
        self.covalent_radius = self.get_covalent_radius()
        self.vdw_radius = self.get_vdw_radius()

        self.coordinates = coordinates
        self.centroid = self.get_centroid()
        self.centre_of_mass = self.get_centre_of_mass()
        self.average_radius = self.get_average_radius()
        self.std_of_radius = self.get_std_of_radius()
        self.distance_list = self.get_distance_list()

        if name is None:
            self.name = ''
        else:
            self.name = name
        if title is None:
            self.title = ''
        else:
            self.title = title
        if fragments is None:
            self.fragments = ()
        else:
            self.fragments = fragments
            self.fragments_coordinates = self.split_coordinates()
            self.fragments_atoms_list = self.split_atoms_lists()

    def __str__(self):
        return "Name: {}\n Coordinates:{}".format(self.name, self.coordinates)

    def __repr__(self):
        return "Molecule.from_xyz('{}')".format(self.name + '.xyz')

    def __iter__(self):
        pass

    def __len__(self):
        return self.number_of_atoms

    def __add__(self, other):
        atoms_list = self.atoms_list + other.atoms_list
        coordinates = np.concatenate((self.coordinates, other.coordinates), axis=0)
        merged = Molecule(atoms_list, coordinates)
        atoms_in_self = [i for i in range(self.number_of_atoms)]
        atoms_in_other = [i for i in range(self.number_of_atoms, merged.number_of_atoms)]
        merged.fragments = [atoms_in_self, atoms_in_other]
        merged.fragments_coordinates = [self.coordinates, other.coordinates]
        merged.fragments_atoms_list = [self.atoms_list, other.atoms_list]
        if hasattr(self, 'fragments_history'):
            merged.fragments_history = self.fragments_history + [atoms_in_other]
        else:
            merged.fragments_history = [atoms_in_self, atoms_in_other]
        return merged

    def split_coordinates(self, coordinates=None):
        if coordinates is None:
            coordinates = self.coordinates
        fragments_coordinates = [coordinates[flist, :] for flist in self.fragments]
        return fragments_coordinates

    def split_atoms_lists(self):
        fragments_atoms_list = [self.atoms_list[flist, :] for flist in self.fragments]
        return fragments_atoms_list

    @classmethod
    def from_xyz(cls, filename):
        import sys
        with open(filename) as fp:
            f = fp.readlines()
        try:
            number_of_atoms = int(f[0])
        except:
            molecule_logger.error("%s should have number of atoms in the first line" % filename)
            molecule_logger.error("but we found\n %s" % f[0])
            molecule_logger.error("Is it an xyz file?")
            sys.exit('Error in reading %s' % filename)
        mol_title = f[1].rstrip()
        try:
            geometry_section = [each_line.split() for each_line in f[2:] if len(each_line) >= 4]
        except:
            molecule_logger.error("Something wrong with reading the geometry section")
            sys.exit('Error in reading %s' % filename)
        if len(geometry_section) != number_of_atoms:
            molecule_logger.error("Number of geometric coordinates is not equal to number of atoms")
            molecule_logger.error("Is something wrong?")
            sys.exit('Error in reading %s' % filename)
        atoms_list = []
        coordinates = []
        for i, c in enumerate(geometry_section):
            try:
                symbol = c[0].capitalize()
                x_coord = float(c[1])
                y_coord = float(c[2])
                z_coord = float(c[3])
            except:
                molecule_logger.error("Something wrong in line: %d" % (i + 1))
                molecule_logger.error(c)
                sys.exit('Error in reading %s' % filename)
            atoms_list.append(symbol)
            coordinates.append([x_coord, y_coord, z_coord])

        mol_coordinates = np.array(coordinates)
        mol_name = filename[:-4]
        return cls(atoms_list, mol_coordinates, name=mol_name, title=mol_title)

    def mol_to_xyz(self, file_name):
        """

        :param file_name: Output .xyz file
        """
        if not hasattr(self, 'energy'):
           self.energy = 0.0
        with open(file_name, 'w') as fp:
            fp.write("{:3d}\n".format(self.number_of_atoms))
            fp.write("{}: {}\n".format(self.title, self.energy))
            for l, c in zip(self.atoms_list, self.coordinates):
                fp.write("%-2s%12.5f%12.5f%12.5f\n" % (l, c[0], c[1], c[2]))

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
                    coords[i][0], coords[i][1], coords[i][2], atoms_list[i].lower()))
            fp.write("$end\n")

    def get_atomic_number(self):
        return [atomic_numbers[c.capitalize()] for c in self.atoms_list]

    def get_atomic_mass(self):
        return [atomic_masses[i] for i in self.atomic_number]

    def get_covalent_radius(self):
        return [covalent_radii[i] for i in self.atomic_number]

    def get_vdw_radius(self):
        result = []
        for i in self.atomic_number:
            if np.isnan(vdw_radii[i]):
                result.append(covalent_radii[i] * 1.5)
            else:
                result.append(vdw_radii[i])
        return result

    def get_centroid(self):
        return self.coordinates.mean(0)

    def get_centre_of_mass(self):
        return np.average(self.coordinates, axis=0, weights=self.atomic_mass)

    def get_average_radius(self):
        return np.mean(np.array([distance(self.centroid, coord_i)
                                 for coord_i in self.coordinates]))

    def get_std_of_radius(self):
        return np.std([distance(self.centroid, coord_i)
                       for coord_i in self.coordinates])

    @property
    def moments_of_inertia_tensor(self):
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

    def get_principal_axes(self):
        return np.linalg.eigvals(self.moments_of_inertia_tensor)

    def get_distance_list(self):
        return scipy_distance.pdist(self.coordinates)

    def get_distance_matrix(self):
        dm = np.zeros((len(self.coordinates), len(self.coordinates)))
        for i, c in enumerate(self.coordinates):
            for j, d in enumerate(self.coordinates):
                dm[i,j] = np.linalg.norm(c-d)
        return dm

    def get_bond_matrix(self):
        """return bond matrix"""
        dm = self.get_distance_matrix()
        bm = np.zeros((len(self.coordinates), len(self.coordinates)), dtype=int)
        for i, c in enumerate(self.coordinates):
            for j, d in enumerate(self.coordinates):
                socr = (self.covalent_radius[i] + self.covalent_radius[j])*1.3
                if i != j and dm[i,j] < socr:
                    bm[i,j] = 1
        return bm

    def calculate_angle(self, i, j, k):

        a1 = self.coordinates[i]
        b1 = self.coordinates[k]
        c1 = self.coordinates[j]

        v1 = c1 - b1
        v2 = c1 - a1
        a = np.arccos(np.dot(v1, v2) / (
                    np.linalg.norm(v1) * np.linalg.norm(v2))) * 180 / pi

        return a

    def calculate_dihedral(self, i, j, k, l):
        p0 = self.coordinates[i]
        p1 = self.coordinates[j]
        p2 = self.coordinates[k]
        p3 = self.coordinates[l]

        b0 = -1.0 * (p1 - p0)
        b1 = p2 - p1
        b2 = p3 - p2
        b1 /= np.linalg.norm(b1)
        # vector rejections
        # v = projection of b0 onto plane perpendicular to b1
        #   = b0 minus component that aligns with b1
        # w = projection of b2 onto plane perpendicular to b1
        #   = b2 minus component that aligns with b1
        v = b0 - np.dot(b0, b1) * b1
        w = b2 - np.dot(b2, b1) * b1
        # angle between v and w in a plane is the torsion angle
        # v and w may not be normalized but that's fine since tan is y/x
        x = np.dot(v, w)
        y = np.dot(np.cross(v, b1), w)
        d = np.degrees(np.arctan2(y, x))
        return d

    def make_internal_coordinates(self):
        """

        :rtype: list of internal coordinates
        """
        dm = self.get_distance_matrix()
        bm = self.get_bond_matrix()
        bl = []
        for i in range(len(self.coordinates)):
            for j in range(len(self.coordinates)):
                if bm[i,j]:
                    seq = [i, j]
                    if seq[::-1] not in [i[0] for i in bl]:
                        bl.append([seq, dm[i,j]])
        al = []
        for i, j, k in product(range(len(self.coordinates)), repeat=3):
           if i != j and j != k and i != k:
                if bm[i,j] and bm[j,k]:
                    seq = [i, j, k]
                    if seq[::-1] not in [i[0] for i in al]:
                        al.append([seq,self.calculate_angle(i, j, k)])

        dl = []
        for i, j, k, l in product(range(len(self.coordinates)), repeat=4):
            if i !=j and i != k and i != l and j != k and j != l and k != l:
                if bm[i,j] and bm[j,k] and bm[k,l]:
                    seq = [i, j, k, l]
                    if seq[::-1] not in [i[0] for i in dl]:
                        dl.append([seq, self.calculate_dihedral(i, j, k, l)])

        return [bl, al, dl]

    def split_covalent_radii_list(self):
        return [np.array(self.covalent_radius)[flist] for flist in self.fragments]

    def is_bonded(self):
        fragment_one, fragment_two = self.split_coordinates()
        radius_one, radius_two = self.split_covalent_radii_list()
        if isinstance(radius_one, np.float) and isinstance(radius_two, np.float64):
            R = radius_one + radius_two
            r = np.linalg.norm(fragment_one-fragment_two)
            if r < R:
                return True
            else:
                return False
        else:
            R = [x + y for x, y in itertools.product(radius_one, radius_two)]
            r = [np.linalg.norm(a - b) for a, b in itertools.product(fragment_one, fragment_two)]
            for l, m in zip(r, R):
                if l < m:
                    return True
            else:
                return False

    @property
    def coulomb_matrix(self):
        """ 
        Author: Debankur Bhattacharyya
        """
        charges = [atomic_numbers[c.capitalize()] for c in self.atoms_list]
        N = self.number_of_atoms
        coords = self.coordinates
        coulomb_matrix = np.zeros((N, N))
        for i in range(N):
            for j in range(N):
                if i == j:
                    coulomb_matrix[i, j] = 0.5 * (charges[i] ** 2.4)
                else:
                    R_ij = np.linalg.norm(coords[i, :] - coords[j, :])
                    coulomb_matrix[i, j] = (charges[i] * charges[j]) / R_ij
        return coulomb_matrix

    @property
    def fingerprint(self):
        eigenvalues = np.linalg.eigvals(self.coulomb_matrix)
        eigenvalues[::-1].sort()
        return eigenvalues

    @property
    def sorted_coulomb_matrix(self):
        """
        From: https://github.com/pythonpanda/coulomb_matrix/
        Takes in a Coulomb matrix of (mxn) dimension and performs a
        row-wise sorting such that ||C(j,:)|| > ||C(j+1,:)||, J=
        0,1,.......,(m-1) Finally returns a vectorized (m*n,1) column 
        matrix .
        """
        summation = np.array([sum(x ** 2) for x in self.coulomb_matrix])
        sorted_matrix = self.coulomb_matrix[np.argsort(summation)[::-1, ], :]
        return sorted_matrix.ravel()

    def rotate_3d(self, vector):
        """This function will rotate a molecule by Euler rotation theorem.
        The Z-X-Z' rotation convention is followed. The range of phi, theta
        and psi should be (0,360),(0,180) and (0,360) respectively.
        This function will first translate the molecule to its center
        of mass(centroid). Then, it rotate the molecule and translate
        to its original position. So, to rotate a molecule around the origin,
        (0.,0.,0.), set_origin usage is necessary"""
        phi, theta, psi = vector
        D = np.array(((cos(phi), sin(phi), 0.),
                      (-sin(phi), cos(phi), 0.),
                      (0., 0., 1.)))
        C = np.array(((1., 0., 0.),
                      (0., cos(theta), sin(theta)),
                      (0., -sin(theta), cos(theta))))
        B = np.array(((cos(psi), sin(psi), 0.),
                      (-sin(psi), cos(psi), 0.),
                      (0., 0., 1.)))
        A = np.dot(B, np.dot(C, D))
        new_coordinates = np.dot(A, np.transpose(self.coordinates))
        self.coordinates = np.transpose(new_coordinates) + self.centroid
        return self

    def move_to_origin(self):
        self.translate(self.get_centroid())
        return self

    def move_to_centre_of_mass(self):
        self.translate(self.get_centre_of_mass())
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
        transformed_coordinates = np.array(np.dot(self.coordinates, eigen_vectors))
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


def distance(coords_a, coords_b):
    return np.linalg.norm(coords_a - coords_b)


def main():
    import sys
    molfile1 = sys.argv[1]
    molobj = Molecule.from_xyz(molfile1)
    for i in molobj.make_internal_coordinates():
        for j in i:
            print(j)

if __name__ == '__main__':
    main()


