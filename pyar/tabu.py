#!/usr/bin/env python3
import copy
import itertools
import logging

import numpy as np
from numpy import pi, cos, sin

tabu_logger = logging.getLogger('pyar.tabu')


def polar_to_cartesian(r, theta, phi):
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
    :type mol_1: Molecule.Molecule
    :type mol_2: Molecule.Molecule
    :type factor: float
    """
    fragment_one, fragment_two = mol_1.coordinates, mol_2.coordinates
    radius_one, radius_two = mol_1.covalent_radius, mol_2.covalent_radius
    for i in range(len(fragment_one)):
        for j in range(len(fragment_two)):
            interatomic_distance = np.linalg.norm(fragment_one[i] - fragment_two[j])
            sum_of_radii = (radius_one[i] + radius_two[j]) * factor
            if interatomic_distance < sum_of_radii:
                return True
    else:
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


def merge_two_molecules(vector, seed_input, monomer_input,
                        freeze_fragments=False, site=None):
    x, y, z, theta, phi, psi = vector
    translate_by = np.array([x, y, z]) / 10
    seed = copy.deepcopy(seed_input)
    monomer = copy.deepcopy(monomer_input)
    tabu_logger.debug('Merging two molecules')
    if freeze_fragments is False:
        seed.move_to_origin()
    monomer.move_to_origin()

    if monomer.number_of_atoms > 1:
        monomer.rotate_3d((theta, phi, psi))
    tabu_logger.debug('checking close contact')

    is_in_cage = True
    while check_close_contact(seed, monomer, 2.3) or is_in_cage:
        minimum_sep_1 = minimum_separation(seed, monomer)
        monomer.translate(translate_by)
        minimum_sep_2 = minimum_separation(seed, monomer)
        if minimum_sep_2 > minimum_sep_1:
            if is_in_cage is True:
                is_in_cage = 'BorW'
            if is_in_cage == 'BorW':
                is_in_cage = False

    orientation = seed + monomer
    if site is not None:
        atoms_in_self = site[0]
        atoms_in_other = site[1]
    else:
        atoms_in_self = [i for i in range(seed.number_of_atoms)]
        atoms_in_other = [i for i in range(seed.number_of_atoms, orientation.number_of_atoms)]
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
        if tabu is True and angle_tabu is True:
            delta_theta = abs(each_saved_entry[3] - point_n_angle[3])
            delta_phi = abs(each_saved_entry[4] - point_n_angle[4])
            delta_psi = abs(each_saved_entry[5] - point_n_angle[5])
            if delta_theta < a_threshold and delta_phi < a_threshold and delta_psi < a_threshold:
                tabu = True
            else:
                tabu = False
    return tabu


def generate_composite_molecule(seed, monomer, points_and_angles):
    composite = copy.deepcopy(seed)
    for vector in points_and_angles:
        composite = merge_two_molecules(vector, composite, monomer,
                                        freeze_fragments=False)
    return composite


def distribute_points_uniformly(points_and_angles):
    """
        Code from Simon Tatham
    https://www.chiark.greenend.org.uk/~sgtatham/polyhedra/
    modified

    :param points_and_angles:
    :return:
    """
    angles = points_and_angles[:, 3:]
    points = points_and_angles[:, :3]
    for _ in range(100000):
        forces = []
        number_of_points = len(points)
        for i in range(number_of_points):
            p = points[i]
            forces_on_p = np.zeros(3)
            for j in range(number_of_points):
                if i == j:
                    continue
                q = points[j]
                v = p - q
                r = np.linalg.norm(v)
                f = v / r ** 3
                forces_on_p += f
            forces.append(forces_on_p)
        forces = np.array(forces)
        total_forces = np.sqrt(np.sum(forces ** 2))
        if total_forces > 0.25:
            scale_force = 0.25 / total_forces
        else:
            scale_force = 1
        dist = 0
        for i in range(len(points)):
            p = points[i]
            f = forces[i]
            moved_point = (p + f * scale_force)
            moved_point /= np.linalg.norm(moved_point)
            dist += np.linalg.norm(p - moved_point)
            points[i] = moved_point

        if dist < 1e-6:
            return np.concatenate((points, angles), axis=1)


def make_5d_coords(grid_location):
    xyz = np.random.uniform(-0.5, 0.5, size=3)
    xyz += grid_location
    xyz /= np.linalg.norm(xyz, axis=0)
    theta_phi = np.random.uniform(0, 2 * pi, size=3)
    return np.concatenate((xyz, theta_phi), axis=0)


def new_gen(number_of_orientations, tabu_on, grid_on, check_angle, d_threshold=0.3, a_threshold=15.0):
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
    :type seed: object
    :param monomer: monomer Molecule object
    :type monomer: object
    :param number_of_orientations: Number of trial geometries
    :type number_of_orientations: int
    :param site:
    :type site: list[int, int] or None
    :return: A list of trial geometries
    :rtype: list
    """
    if monomer.number_of_atoms == 1:
        tabu_check_for_angles = False
    else:
        tabu_check_for_angles = True

    tabu_logger.debug('Generating points')
    points_and_angles = new_gen(number_of_orientations, tabu_on,
                                grid_on, tabu_check_for_angles)
    tabu_logger.debug('Generated points')
    write_tabu_list(points_and_angles, 'tabu.dat')
    # plot_points(pts)

    tabu_logger.debug('generate orientations from points and angles')

    orientations = []
    filename_prefix = 'trial_'
    for i, vector in enumerate(points_and_angles):
        new_orientation = merge_two_molecules(vector, seed, monomer, site=site)
        new_orientation_id = f"{i:03d}_{molecule_id}"
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
    ax.plot_wireframe(x, y, z, color='blue', rstride=1, cstride=1, linewidth=0.1)
    ax.scatter(points[:, 0], points[:, 1], points[:, 2], s=100, c='r', zorder=10)
    fig.show()
    fig.savefig(f'{location}/points.png')


def main():
    pass


if __name__ == "__main__":
    pts = new_gen(4, True, True, False, 0.3, 15.0)
    for i in pts:
        print(i)
    plot_points(pts, '/home/anoop')


def write_tabu_list(tabu_list, tabu_file):
    with open(tabu_file, 'a') as tf:
        for line in tabu_list:
            tf.write(str(line) + ' ')
        tf.write('\n')
