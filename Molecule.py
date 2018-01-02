from math import radians, cos, sin
from scipy.spatial import distance as scipy_distance

import numpy as np

from atomic_data import atomic_masses, atomic_numbers, covalent_radii


class Molecule(object):
    def __init__(self, atoms_list, coordinates, name=None, title=None, fragments=None):
        self.number_of_atoms = len(coordinates)
        self.atoms_list = atoms_list
        self.coordinates = coordinates
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

    def __str__(self):
        return "Name: {}\n Coordinates:{}".format(self.name, self.coordinates)

    def __repr__(self):
        return "Molecule.from_xyz('{}')".format(self.name + '.xyz')

    def __iter__(self):
        pass

    def __len__(self):
        return self.number_of_atoms

    @classmethod
    def from_xyz(cls, filename):
        import sys
        with open(filename) as fp:
            f = fp.readlines()
        try:
            number_of_atoms = int(f[0])
        except:
            print(filename, "should have number of atoms in the first line")
            print("but we found\n", f[0])
            print("Is it an xyz file?")
            sys.exit()
        mol_title = f[1].rstrip()
        try:
            geometry_section = [each_line.split() for each_line in f[2:] if len(each_line) >= 4]
        except:
            print("Something wrong with reading the geometry section")
            sys.exit()
        if len(geometry_section) != number_of_atoms:
            print("Number of geometric coordinates is not equal to number of atoms")
            print("Is something wrong?")
            sys.exit()
        atoms_list = []
        coordinates = []
        for i, c in enumerate(geometry_section):
            try:
                symbol = c[0]
                x_coord = float(c[1])
                y_coord = float(c[2])
                z_coord = float(c[3])
            except:
                print("Something wrong in line: ", i + 1)
                print(c)
                sys.exit()
            atoms_list.append(symbol)
            coordinates.append([x_coord, y_coord, z_coord])

        mol_coordinates = np.array(coordinates)
        mol_name = filename[:-4]
        return cls(atoms_list, mol_coordinates, name=mol_name, title=mol_title)

    def mol_to_xyz(self, file_name):
        """

        :param file_name: Output .xyz file
        """
        with open(file_name, 'w') as fp:
            fp.write("{:3d}\n".format(self.number_of_atoms))
            fp.write("{}\n".format(self.title))
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

    def __add__(self, molecule):
        atoms_list = self.atoms_list + molecule.atoms_list
        coordinates = np.concatenate((self.coordinates, molecule.coordinates), axis=0)
        merged = Molecule(atoms_list, coordinates)
        atoms_in_self = [i for i in range(self.number_of_atoms)]
        atoms_in_molecule = [i for i in range(self.number_of_atoms, merged.number_of_atoms)]
        merged.fragments = [atoms_in_self, atoms_in_molecule]
        return merged

    @property
    def atomic_number(self):
        return [atomic_numbers[c.capitalize()] for c in self.atoms_list]

    @property
    def atomic_mass(self):
        return [atomic_masses[i] for i in self.atomic_number]

    @property
    def covalent_radius(self):
        return [covalent_radii[i] for i in self.atomic_number]

    @property
    def centroid(self):
        return self.coordinates.mean(0)

    @property
    def centre_of_mass(self):
        return np.average(self.coordinates, axis=0, weights=self.atomic_mass)

    @staticmethod
    def distance(coords_a, coords_b):
        return np.linalg.norm(coords_a - coords_b)

    @property
    def average_radius(self):
        return np.mean(np.array([self.distance(self.centroid, coord_i)
                                 for coord_i in self.coordinates]))

    @property
    def std_of_radius(self):
        return np.std([self.distance(self.centroid, coord_i)
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

    @property
    def principal_axes(self):
        return np.linalg.eigvals(self.moments_of_inertia_tensor)

    @property
    def distance_matrix(self):
        return scipy_distance.pdist(self.coordinates)

    def is_bonded(self):
        fragment_one = np.array([self.coordinates[i, :] for i in self.fragments[0]])
        fragment_two = np.array([self.coordinates[i, :] for i in self.fragments[1]])
        radius_one = [covalent_radii[i] for i in self.fragments[0]]
        radius_two = [covalent_radii[i] for i in self.fragments[1]]

        for i, a in enumerate(fragment_one):
            for j, b in enumerate(fragment_two):
                dist = np.linalg.norm(a - b)
                sum_of_radii = radius_one[i] + radius_two[j]
                if dist <= sum_of_radii:
                    return True
        else:
            return False

    @property
    def coulomb_matrix(self):
        """ 
        Author: Debankur Bhattacharyya
        """
        charges = [atomic_numbers[c] for c in self.atoms_list]
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
        Takes in a Coloumb matrix of (mxn) dimension and performs a 
        rowwise sorting such that ||C(j,:)|| > ||C(j+1,:)||, J= 
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
        phi, theta, psi = map(radians, vector)
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

    def move_to_origin(self):
        self.translate(self.centroid)

    def move_to_centre_of_mass(self):
        self.translate(self.centre_of_mass)

    def translate(self, magnitude):
        self.coordinates = self.coordinates - magnitude

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
        return np.dot(transformed_coordinates, move_axes)


def main():
    import sys
    mol = Molecule.from_xyz(sys.argv[1])
    print(mol)


if __name__ == '__main__':
    main()
