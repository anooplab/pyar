"""
tabu.py

Functions to merge two molecules.
"""
import collections
import copy
import itertools
import logging

import numpy as np
from numpy import pi, cos, sin

from pyar.Molecule import Molecule
# from pyar.property import get_connectivity
import networkx as nx
from typing import List, Tuple

tabu_logger = logging.getLogger('pyar.tabu')


def polar_to_cartesian(r, theta, phi):
    """

    :param r: distance r
    :param theta: angle theta
    :param phi: angle phi
    :return: (x, y, z)
    """
    x = r * sin(theta) * cos(phi)
    y = r * sin(theta) * sin(phi)
    z = r * cos(theta)
    return x, y, z


def check_close_contact(mol_1, mol_2, factor):
    """
    Checks for close contacts between two molecules.  If the distance between
    any pair of atoms is less than the sum of van der Waals radii (scaled by
    factor), then it is assumed to have a close contact.

    :return: boolean
    :type mol_1: Callable[..., Molecule]
    :type mol_2: () -> Molecule
    :type factor: float
    """
    fragment_one, fragment_two = mol_1.coordinates, mol_2.coordinates
    radius_one, radius_two = mol_1.covalent_radius, mol_2.covalent_radius
    for f1, r1 in zip(fragment_one, radius_one):
        for f2, r2 in zip(fragment_two, radius_two):
            interatomic_distance = np.linalg.norm(f1 - f2)
            sum_of_radii = (r1 + r2) * factor
            if interatomic_distance < sum_of_radii:
                return True
    return False


def minimum_separation(mol_1, mol_2):
    """
    Find the minimum separation between two fragments.

    :return: boolean
    :type mol_1: Molecule.Molecule
    :type mol_2: Molecule.Molecule
    """
    return np.min([np.linalg.norm(a - b) for a, b in
                   itertools.product(mol_1.coordinates, mol_2.coordinates)])


def merge_two_molecules(vector, seed_input, monomer_input, freeze_fragments=False, site=None,  distance_scaling=1.5):
    """

    :param ndarray vector: the direction for placing the monomer
    :param seed_input: Seed () -> Molecule
    :type seed_input: Seed molecule
    :param molecule monomer_input: Monomer molecule
    :param bool freeze_fragments: Whether to rotate the monomer
    :param site: weighted atom centres for placing
    :type site: Union[list, None]
    :param float distance_scaling: minimum separation between the
        two fragments in merged molecule is sum_of_covalent_radii
        * distance_scaling
    :return: merged molecule
    """

    x, y, z, theta, phi, psi = vector
    direction = np.array([x, y, z])
    tiny_steps = direction / 10
    seed = copy.deepcopy(seed_input)
    monomer = copy.deepcopy(monomer_input)
    tabu_logger.debug('Merging two molecules')
    if freeze_fragments is False:
        seed.move_to_origin()
    monomer.move_to_origin()

    if monomer.number_of_atoms > 1:
        monomer.rotate_3d((theta, phi, psi))
    tabu_logger.debug('checking close contact')

    tabu_logger.debug('Taking a large first step')
    r_max_of_seed = np.max([np.linalg.norm(c) for c in seed.coordinates])
    r_max_of_monomer = np.max([np.linalg.norm(c) for c in monomer.coordinates])
    maximum_distance_to_move = r_max_of_seed + r_max_of_monomer + 1.0
    move_to = direction * maximum_distance_to_move
    monomer.translate(move_to)

    if contact := check_close_contact(seed, monomer, distance_scaling):
        tabu_logger.debug("Large step was not enough")
        while contact:
            tabu_logger.debug("moving by small steps")
            monomer.translate(tiny_steps)
            contact = check_close_contact(seed, monomer, distance_scaling)
    else:
        while contact:
            monomer.translate(-1 * tiny_steps)
            contact = check_close_contact(seed, monomer, distance_scaling)

    #     lower_distance = np.zeros(3)
    #     upper_distance = move_to
    #     move_to = (upper_distance + lower_distance) / 2
    #     tabu_logger.debug("Binary steps")
    #     for _ in range(100):
    #         monomer.move_to_origin()
    #         monomer.translate(move_to)
    #         step_counter += 1
    #         contact = check_close_contact(seed, monomer, distance_scaling)
    #         if contact:
    #             lower_distance = move_to
    #         else:
    #             upper_distance = move_to
    #         move_to = (upper_distance + lower_distance) / 2
    #         if np.linalg.norm(upper_distance - lower_distance) < 0.01:
    #             break
    # tabu_logger.debug(f"Total steps taken {step_counter}")

    # is_in_cage = True
    # while check_close_contact(seed, monomer, distance_scaling) or is_in_cage:
    #     print("in loop")
    #     minimum_sep_1 = minimum_separation(seed, monomer)
    #     monomer.translate(translate_by)
    #     minimum_sep_2 = minimum_separation(seed, monomer)
    #     if minimum_sep_2 > minimum_sep_1:
    #         if is_in_cage is True:
    #             is_in_cage = 'BorW'
    #         if is_in_cage == 'BorW':
    #             is_in_cage = False

    orientation = seed + monomer
    if site is not None:
        atoms_in_self = site[0]
        atoms_in_other = site[1]
    else:
        atoms_in_self = list(range(seed.number_of_atoms))
        atoms_in_other = list(range(seed.number_of_atoms,
                                    orientation.number_of_atoms))
    orientation.fragments = [atoms_in_self, atoms_in_other]
    tabu_logger.debug('Merged.')
    return orientation


def check_tabu_status(point_n_angle, d_threshold, a_threshold, tabu_list,
                      angle_tabu):
    """
    Check if the given point (point_n_angle) is within the proximity (set by
    d_threshold and a_threshold

    :type point_n_angle: Numpy.ndarray
    :param point_n_angle: six numbers (x, y, z, theta, phi, psi)
    :type d_threshold: float
    :param d_threshold: minimum distance for Tabu
    :type a_threshold: float
    :param a_threshold: minimum angle (in radian) for Tabu
    :type tabu_list: ndarray
    :param tabu_list: the saved points and angles
    :param bool angle_tabu: boolean if True consider also angles for checking the proximity
    :type angle_tabu: bool
    :param angle_tabu: boolean if True consider also angles for checking the proximity
    :rtype: bool
    :return: True if the new points is within the proximity of saved points
    """
    if tabu_list is None:
        return False
    tabu = False
    for each_saved_entry in tabu_list:
        distance = np.linalg.norm(each_saved_entry[:3] - point_n_angle[:3])
        if distance < d_threshold:
            tabu = True
        if tabu and angle_tabu is True:
            delta_theta = abs(each_saved_entry[3] - point_n_angle[3])
            delta_phi = abs(each_saved_entry[4] - point_n_angle[4])
            delta_psi = abs(each_saved_entry[5] - point_n_angle[5])
            tabu = (
                    delta_theta < a_threshold
                    and delta_phi < a_threshold
                    and delta_psi < a_threshold
            )

    return tabu


def create_composite_molecule(seed, monomer, points_and_angles, d_scale):
    """

    :param seed: seed molecule
    :type seed: Molecule.Molecule
    :type monomer: Molecule.Molecule
    :param monomer: monomer molecule
    :param points_and_angles: a set of 3 points and 3 angles
    :type points_and_angles: ndarray
    :param d_scale: the minimum distance between two fragment is scaled by
        sumo of covalent radii * d_scale
    :type d_scale: float
    :returns a composite molecule
    :rtype Molecule.Molecule
    """
    composite = copy.deepcopy(seed)
    for vector in points_and_angles:
        composite = merge_two_molecules(vector, composite, monomer,
                                        freeze_fragments=False,
                                        distance_scaling=d_scale)
    return composite


def distribute_points_uniformly(points_and_angles):
    """
    Code from Simon Tatham
    https://www.chiark.greenend.org.uk/~sgtatham/polyhedra/
    modified

    :param points_and_angles:
    :return: ndarray containing three points and angles
    :rtype ndarray
    """
    angles = points_and_angles[:, 3:]
    points = points_and_angles[:, :3]
    for _ in range(100000):
        forces = []
        number_of_points = len(points)
        for ith in range(number_of_points):
            p = points[ith]
            forces_on_p = np.zeros(3)
            for jth in range(number_of_points):
                if ith == jth:
                    continue
                q = points[jth]
                v = p - q
                r = np.linalg.norm(v)
                f = v / r ** 3
                forces_on_p += f
            forces.append(forces_on_p)
        forces = np.array(forces)
        total_forces = np.sqrt(np.sum(forces ** 2))
        scale_force = 0.25 / total_forces if total_forces > 0.25 else 1
        dist = 0
        for ith in range(len(points)):
            p = points[ith]
            f = forces[ith]
            moved_point = (p + f * scale_force)
            moved_point /= np.linalg.norm(moved_point)
            dist += np.linalg.norm(p - moved_point)
            points[ith] = moved_point

        if dist < 1e-6:
            return np.concatenate((points, angles), axis=1)


def make_5d_coords(grid_location):
    """

    :type grid_location: ndarray
    :param grid_location: In grid based generation of points, this will choose the grid.
    :return: a set of three points and three angles
    :rtype: ndarray
    """
    xyz = np.random.uniform(-0.5, 0.5, size=3)
    xyz += grid_location
    xyz /= np.linalg.norm(xyz, axis=0)
    theta_phi = np.random.uniform(0, 2 * pi, size=3)
    return np.concatenate((xyz, theta_phi), axis=0)


def generate_points(number_of_orientations, tabu_on, grid_on, check_angle,
                    d_threshold=0.3, a_threshold=15.0):
    """
    Generate points

    :type a_threshold: float
    :param a_threshold:
    :type d_threshold: float
    :param d_threshold:
    :param bool check_angle: Should it check tabu list for angles
    :param int number_of_orientations: Number of orientations
    :param bool grid_on: Use Grid or not
    :param bool tabu_on: Use tabu or not
    :return:
    """

    if grid_on:
        choose_grid = itertools.cycle([(0.5, 0.5, 0.5), (-0.5, -0.5, -0.5),
                                       (0.5, 0.5, -0.5), (-0.5, -0.5, 0.5),
                                       (0.5, -0.5, 0.5), (-0.5, 0.5, -0.5),
                                       (0.5, -0.5, -0.5), (-0.5, 0.5, 0.5)])
    else:
        choose_grid = itertools.cycle([0, 0, 0])

    tabu_list = [make_5d_coords(next(choose_grid))]
    for _ in range(number_of_orientations - 1):
        current_grid = next(choose_grid)
        point_n_angle = make_5d_coords(current_grid)
        if tabu_on:
            tries = 1
            while check_tabu_status(point_n_angle, d_threshold, a_threshold,
                                    tabu_list, check_angle) is True:
                point_n_angle = make_5d_coords(current_grid)
                tries += 1
                if tries > 10000:
                    d_threshold *= 0.95

        tabu_list.append(point_n_angle)
    return np.array(tabu_list)


def create_trial_geometries(molecule_id, seed, monomer,
                            number_of_orientations,
                            tabu_on, grid_on, site):
    """

    :param grid_on: Use grid for making points
    :type grid_on: bool
    :param tabu_on: Use tabu
    :type tabu_on:
    :type molecule_id: str
    :param molecule_id: An ID for molecule
    :param seed: seed Molecule object
    :type seed: Any
    :param monomer: monomer Molecule object
    :type monomer: Any
    :param number_of_orientations: Number of trial geometries
    :type number_of_orientations: int
    :param site:
    :type site: list[int, int] or None
    :return: A list of trial geometries
    :rtype: list
    """
    if monomer.number_of_atoms == 1:
        tabu_check_for_angles = False
        proximity_factor = 1.2
    else:
        tabu_check_for_angles = True
        proximity_factor = 1.5

    tabu_logger.debug('Generating points')
    points_and_angles = generate_points(number_of_orientations, tabu_on,
                                        grid_on, tabu_check_for_angles)
    tabu_logger.debug('Generated points')
    write_tabu_list(points_and_angles, 'tabu.dat')
    # plot_points(pts)

    tabu_logger.debug('generate orientations from points and angles')

    orientations = []
    filename_prefix = 'trial_'
    for counter, vector in enumerate(points_and_angles):
        new_orientation = merge_two_molecules(vector, seed, monomer, site=site,
                                              distance_scaling=proximity_factor)
        new_orientation_id = f"{counter:03d}_{molecule_id}"
        new_orientation.title = f'trial orientation {new_orientation_id}'
        new_orientation.name = new_orientation_id
        new_orientation_xyz_file = f'{filename_prefix}{new_orientation_id}.xyz'
        new_orientation.mol_to_xyz(new_orientation_xyz_file)
        orientations.append(new_orientation)

    return orientations


def plot_points(points, location):
    """have to run with python -i """

    import matplotlib.pyplot as plt
    from matplotlib import style
    style.use('dark_background')

    phi = np.linspace(0, np.pi, 60)
    theta = np.linspace(0, 2 * np.pi, 90)
    x = np.outer(np.sin(theta), np.cos(phi))
    y = np.outer(np.sin(theta), np.sin(phi))
    z = np.outer(np.cos(theta), np.ones_like(phi))
    fig, ax = plt.subplots(1, 1, subplot_kw={'projection': '3d'})
    ax.plot_wireframe(x, y, z, color='blue', rstride=1, cstride=1,
                      linewidth=0.1)
    ax.scatter(points[:, 0], points[:, 1], points[:, 2], s=100, c='r',
               zorder=10)
    fig.show()
    fig.savefig(f'{location}/points.png')


def write_tabu_list(tabu_list, tabu_file):
    """
    Save the tabu_list in tabu_file
    :param tabu_list: The tabu list
    :type tabu_list: list or ndarray
    :param tabu_file: The file name to write tabu file
    :type tabu_file: str
    :rtype: None
    :returns None
    """
    with open(tabu_file, 'a') as tf:
        for line in tabu_list:
            tf.write(f'{str(line)} ')
        tf.write('\n')


def main():
    """Nothing now"""
    pass


def get_connectivity(coordinates: np.ndarray, covalent_radii: List[float],
                     tolerance: float = 0.4) -> List[Tuple[int, int]]:
    """
    Calculate connectivity based on interatomic distances and covalent radii.

    :param coordinates: numpy array of shape (n_atoms, 3) with atomic coordinates
    :param covalent_radii: list of covalent radii for each atom
    :param tolerance: tolerance factor for bond detection
    :return: list of tuples representing bonds (atom_index1, atom_index2)
    """
    n_atoms = len(coordinates)
    bonds = []
    for i in range(n_atoms):
        for j in range(i+1, n_atoms):
            distance = np.linalg.norm(coordinates[i] - coordinates[j])
            if distance < (covalent_radii[i] + covalent_radii[j]) * (1 + tolerance):
                bonds.append((i, j))
    return bonds

def broken(molobj) -> bool:
    """
    Check if the molecule is fragmented.

    :param molobj: object(Molecule)
    :return: Is the molecule fragmented?
    :rtype: bool
    """
    # Create a graph from the molecular structure
    G = nx.Graph()
    G.add_nodes_from(range(len(molobj.atoms_list)))
    G.add_edges_from(get_connectivity(molobj.coordinates, molobj.covalent_radius))

    # Check if the graph is connected
    return not nx.is_connected(G)

if __name__ == "__main__":
    import sys


    input_xyz = sys.argv[1]
    mol = Molecule.from_xyz(input_xyz)
    if broken(mol):
        print("is a fragment")
    else:
        print("Fine")
