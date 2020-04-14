#!/usr/bin/env python3
import copy
import itertools
import logging
import os.path
import random
import time

import numpy as np
from numpy import pi, cos, sin, degrees

tabu_logger = logging.getLogger('pyar.tabu')


def load_tabu_list(tabu_file='tabu.dat'):
    if not os.path.exists(tabu_file):
        return []
    with open(tabu_file, 'r') as tf:
        tabu_list = [map(float, i.split()) for i in tf.readlines()]
        return tabu_list


def write_tabu_list(tabu_list, tabu_file):
    with open(tabu_file, 'a') as tf:
        for i in tabu_list:
            tf.write(str(i) + ' ')
        tf.write('\n')


def gen_a_set_of_angles(which_octant):
    sa = 0.0
    ri = pi / 2.0
    ga = pi
    ma = 3 * pi / 2.0
    pa = 2 * pi

    octant = {
        1: [sa, ri, sa, ri], 2: [ri, ga, sa, ri], 3: [sa, ri, ri, ga], 4: [ri, ga, ri, ga],
        5: [sa, ri, ga, ma], 6: [ri, ga, ga, ma], 7: [sa, ri, ma, pa], 8: [ri, ga, ma, pa]
    }

    do, rea, me, fa = octant[which_octant]

    theta = random.uniform(do, rea)
    phi = random.uniform(me, fa)
    alpha = random.uniform(sa, pa)
    beta = random.uniform(sa, ga)
    gamma = random.uniform(sa, pa)

    return list(map(degrees, [theta, phi, alpha, beta, gamma]))


def polar_to_cartesian(r, theta, phi):
    x = r * sin(theta) * cos(phi)
    y = r * sin(theta) * sin(phi)
    z = r * cos(theta)
    return x, y, z


def check_similarity(angles, previous_angles, d_threshold, a_threshold):
    r = 1.0
    theta, phi, alpha, beta, gamma = angles
    current_point = np.array(polar_to_cartesian(r, theta, phi))
    for i in previous_angles:
        previous_theta, previous_phi = i[:2]
        previous_alpha, previous_beta, previous_gamma = i[2:]
        previous_point = np.array(polar_to_cartesian(r, previous_theta, previous_phi))
        distance = np.linalg.norm(current_point - previous_point)
        diff_alpha = abs(previous_alpha - alpha)
        diff_beta = abs(previous_beta - beta)
        diff_gamma = abs(previous_gamma - gamma)
        if abs(distance) < d_threshold and \
                diff_alpha < a_threshold and \
                diff_beta < a_threshold and \
                diff_gamma < a_threshold:
            return False, angles
    return True, angles


def gen_vectors(number_of_orientations):
    vectors = []
    a_threshold = pi / 2.0
    d_threshold = 1.414
    for i in range(int(int(number_of_orientations) / 8)):
        for j in range(8):
            accepted = False
            one_set = None
            while not accepted:
                accepted, one_set = check_similarity(gen_a_set_of_angles(j + 1), vectors, d_threshold, a_threshold)
            vectors.append(one_set)
        d_threshold *= 0.95
    return vectors


def generate_orientations_old(molecule_id, seed, monomer, number_of_orientations):
    number_of_atoms_in_seed = seed.number_of_atoms
    number_of_atoms_in_monomer = monomer.number_of_atoms
    filename_prefix = 'trial_'
    if number_of_atoms_in_seed > 1:
        seed.align()
    if number_of_atoms_in_monomer > 1:
        monomer.align()
    seed.move_to_origin()

    orientations = []

    if number_of_atoms_in_seed == 1 and number_of_atoms_in_monomer == 1:
        vector = [0.0, 0.0, 0.0, 0.0, 0.0]
        orientation = merge_monomer_and_seed(vector, monomer, seed)
        orientation_id = "%03d_" % 0 + molecule_id
        orientation.title = 'trial orientation ' + orientation_id
        orientation.name = orientation_id

        orientation_xyz_file = filename_prefix + orientation_id + '.xyz'
        orientation.mol_to_xyz(orientation_xyz_file)
        orientations.append(orientation)
    else:
        vectors = []
        a_threshold = pi / 2.0
        d_threshold = 1.414
        for i in range(int(int(number_of_orientations) / 8)):
            for j in range(8):
                accepted = False
                one_set = None
                while not accepted:
                    accepted, one_set = check_similarity(gen_a_set_of_angles(j + 1), vectors, d_threshold, a_threshold)
                orientation = merge_monomer_and_seed(one_set, monomer, seed)
                vectors.append(one_set)
                orientation_id = "%03d_" % (i * 8 + j) + molecule_id
                orientation.title = 'trial orientation ' + orientation_id
                orientation.name = orientation_id

                orientation_xyz_file = filename_prefix + orientation_id + '.xyz'
                orientation.mol_to_xyz(orientation_xyz_file)
                orientations.append(orientation)

                d_threshold *= 0.95
                a_threshold *= 0.95
    return orientations


def merge_monomer_and_seed(each_vector, monomer, seed, site=None):
    theta, phi, alpha, beta, gamma = each_vector
    monomer.move_to_origin()
    if monomer.number_of_atoms > 1:
        monomer.rotate_3d((alpha, beta, gamma))

    r = 0.1

    while close_contact(seed, monomer, 1.0):
        r += 0.1
        if seed.number_of_atoms == 2:
            phi = 0.0
        monomer.translate(polar_to_cartesian(r, theta, phi))

    orientation = seed + monomer
    if site is not None:
        atoms_in_self = site
    else:
        atoms_in_self = [i for i in range(seed.number_of_atoms)]
    atoms_in_other = [i for i in range(seed.number_of_atoms, orientation.number_of_atoms)]
    orientation.fragments = [atoms_in_self, atoms_in_other]
    write_tabu_list(each_vector, 'tabu.dat')
    return orientation


def close_contact(mol_1, mol_2, factor):
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


def check_tabu(point_n_angle, d_threshold, a_threshold, saved_points_and_angles,
               angle_tabu):
    """
    Check if the given point (point_n_angle) is within the proximity (set by
    d_threshold and a_threshold

    :param point_n_angle: numpy array of six numbers (x, y, z, theta, phi, psi)
    :param d_threshold: float minimum distance for Tabu
    :param a_threshold: float minimum angle (in radian) for Tabu
    :param saved_points_and_angles: numpy array containing the saved points and
     angles
    :param angle_tabu: boolean if True consider also angles for checking the
    proximity
    :return: boolean True if the new points is within the proximity of saved
    points
    """
    if saved_points_and_angles is None:
        return False
    tabu = False
    for each_saved_entry in saved_points_and_angles:
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


def make_an_untabooed_point(a_threshold, angle_tabu, d_threshold,
                            octant_chooser, saved_points_and_angles):
    """

    Make a new point (x, y, z, theta, phi, psi) which is not in the Tabu list

    :param a_threshold:
    :param angle_tabu:
    :param d_threshold:
    :param octant_chooser:
    :param saved_points_and_angles:
    :return:
    """
    tries = 1
    point_n_angle = make_point_n_angles(octant_chooser)
    while check_tabu(point_n_angle, d_threshold, a_threshold, saved_points_and_angles, angle_tabu) is True:
        point_n_angle = make_point_n_angles(octant_chooser)
        tries += 1
        if tries > 10000:
            d_threshold *= 0.95
    return point_n_angle


def rotating_octant(number_of_points, distance_tabu=True, angle_tabu=True,
                    remember_points_and_angles=None):
    """
    Make N points ((x, y, z, theta, phi, psi) using 'rotating octant' method.

    :type remember_points_and_angles: numpy.array
    :param number_of_points: int
    :param distance_tabu: float
    :param angle_tabu: float
    :return: numpy.array
    """
    if remember_points_and_angles is None:
        saved_points_and_angles = []
    else:
        saved_points_and_angles = remember_points_and_angles[:]

    d_threshold = np.sqrt(2) / 2.0
    a_threshold = pi / 2.0
    for i in range(number_of_points // 8):
        for j in itertools.product([1, -1], repeat=3):
            octant_chooser = np.array(j)
            if distance_tabu is False:
                saved_points_and_angles.append(make_point_n_angles(octant_chooser))
            if distance_tabu is True:
                if len(saved_points_and_angles) == 0:
                    saved_points_and_angles.append(make_point_n_angles(octant_chooser))
                else:
                    point_n_angle = make_an_untabooed_point(a_threshold, angle_tabu, d_threshold, octant_chooser, remember_points_and_angles)
                    saved_points_and_angles.append(point_n_angle)

    return np.array(saved_points_and_angles)


def make_point_n_angles(octant):
    a_point = make_a_random_point(octant)
    two_angles = make_three_random_angles()
    return np.concatenate((a_point, two_angles), axis=0)


def make_three_random_angles():
    return np.random.uniform(0, 2 * pi, size=3)


def make_a_random_point(octant):
    p = np.random.uniform(0, 1.0, size=3)
    p /= np.linalg.norm(p, axis=0)
    p *= octant
    return p


def make_random_points_and_angles(number_of_points):
    """https://stackoverflow.com/questions/33976911/generate-a-random-sample-of-points-distributed-on-the-surface-of-a-unit-sphere"""
    vec_1 = np.random.randn(3, number_of_points)
    vec_1 /= np.linalg.norm(vec_1, axis=0)
    vec_2 = np.random.uniform(2 * pi, size=(number_of_points, 3))
    vec = np.concatenate((vec_1.transpose(), vec_2), axis=1)
    return vec


def uniformly_distributed_points(N):
    """
    Code from Simon Tatham
    https://www.chiark.greenend.org.uk/~sgtatham/polyhedra/
    modified
    """

    points_and_angles = make_random_points_and_angles(N)

    angles = points_and_angles[:, 3:]
    points = points_and_angles[:, :3]

    for _ in range(100000):
        forces = []
        for i in range(N):
            p = points[i]
            forces_on_p = np.zeros(3)
            for j in range(N):
                if i == j: continue
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


def plot_points(pts):
    '''have to run with python -i '''

    import matplotlib.pyplot as plt
    from matplotlib import style
    style.use('dark_background')

    phi = np.linspace(0, np.pi, 60)
    theta = np.linspace(0, 2 * np.pi, 90)
    x = np.outer(np.sin(theta), np.cos(phi))
    y = np.outer(np.sin(theta), np.sin(phi))
    z = np.outer(np.cos(theta), np.ones_like(phi))
    fig, ax = plt.subplots(1, 1, subplot_kw={'projection': '3d', 'aspect': 'equal'})
    ax.plot_wireframe(x, y, z, color='blue', rstride=1, cstride=1, linewidth=0.1)
    ax.scatter(pts[:, 0], pts[:, 1], pts[:, 2], s=100, c='r', zorder=10)
    fig.show()
    fig.savefig('points.png')


def merge_two_molecules(vector, seed_input, monomer_input, freeze_fragments=False, site=None):
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
    while close_contact(seed, monomer, 1.0) or is_in_cage:
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


def generate_orientations_from_points_and_angles(seed, monomer,
                                                 points_and_angles, site=None):
    orientations = []
    tabu_logger.debug('generate orientations from points and angles')
    for vector in points_and_angles:
        orientations.append(merge_two_molecules(vector, seed, monomer, site=site))
    return orientations


def generate_composite_molecule(seed, monomer, points_and_angles):
    composite = copy.deepcopy(seed)
    for vector in points_and_angles:
        composite = merge_two_molecules(vector, composite, monomer,
                                        freeze_fragments=False)
    return composite


def generate_orientations(molecule_id, seed, monomer, number_of_orientations,
                          method='rotating'):
    if monomer.number_of_atoms == 1:
        tabu_check_for_angles = False
    else:
        tabu_check_for_angles = True

    tabu_logger.debug('Generating points')
    pts = wrapper_of_make_points_and_angles(method, number_of_orientations,
                                            tabu_check_for_angles)
    tabu_logger.debug('Generated points')
    write_tabu_list(pts, 'tabu.dat')
    # plot_points(pts)

    filename_prefix = 'trial_'
    orientations = generate_orientations_from_points_and_angles(seed, monomer, pts)
    for i, each_orientation in enumerate(orientations):
        each_orientation_id = "%03d_" % (i) + molecule_id
        each_orientation.title = 'trial orientation ' + each_orientation_id
        each_orientation.name = each_orientation_id
        each_orientation_xyz_file = filename_prefix + each_orientation_id + '.xyz'
        each_orientation.mol_to_xyz(each_orientation_xyz_file)

    return orientations


def wrapper_of_make_points_and_angles(method, number_of_orientations,
                                      tabu_check_for_angles):
    if method == 'rotating':
        pts = rotating_octant(number_of_orientations,
                              angle_tabu=tabu_check_for_angles)
    elif method == 'random':
        pts = make_random_points_and_angles(number_of_orientations)
    elif method == 'uniform':
        pts = uniformly_distributed_points(number_of_orientations)
    else:
        tabu_logger.info('using default method: random')
        pts = make_random_points_and_angles(number_of_orientations)
    return pts


def generate_guess_for_bonding(molecule_id, seed, monomer, a, b,
                               number_of_orientations):
    if monomer.number_of_atoms == 1:
        tabu_check_for_angles = False
    else:
        tabu_check_for_angles = True

    saved_pts = []

    orientations = []
    for i in range(number_of_orientations):
        t1 = time.clock()
        pts = rotating_octant(32, angle_tabu=tabu_check_for_angles,
                              remember_points_and_angles=saved_pts)
        t2 = time.clock()
        tabu_logger.debug('Created points: in {} seconds'.format(t2 - t1))
        t1 = time.clock()
        current_orientations = generate_orientations_from_points_and_angles(seed, monomer, pts, site=[a, b])
        t2 = time.clock()
        tabu_logger.debug('Created orientations {} seconds'.format(t2 - t1))
        t1 = time.clock()
        stored_orientations = {}
        for j, each_orientation in enumerate(current_orientations):
            coords = each_orientation.coordinates
            dist = np.linalg.norm(coords[a] - coords[b])
            stored_orientations[j] = dist
        best_orientation = min(stored_orientations, key=stored_orientations.get)
        best_point = pts[best_orientation]
        tabu_logger.debug("{} {}".format(best_orientation, stored_orientations[best_orientation]))
        saved_pts.append(best_point)
        orientations.append(current_orientations[best_orientation])
        t2 = time.clock()
        tabu_logger.debug('Found best orientation in {} seconds'.format(t2 - t1))

    t1 = time.clock()
    filename_prefix = 'trial_'
    for i, each_orientation in enumerate(orientations):
        each_orientation_id = "%03d_" % (i) + molecule_id
        each_orientation.title = 'trial orientation ' + each_orientation_id
        each_orientation.name = each_orientation_id
        each_orientation_xyz_file = filename_prefix + each_orientation_id + '.xyz'
        each_orientation.mol_to_xyz(each_orientation_xyz_file)
    t2 = time.clock()
    tabu_logger.debug('Wrote files in {} seconds'.format(t2 - t1))
    write_tabu_list(saved_pts, 'tabu.dat')
    plot_points(np.array(saved_pts))

    return orientations


def main():
    pass


if __name__ == "__main__":
    main()
