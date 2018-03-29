from __future__ import print_function

import math
import os.path
import random
import itertools

import numpy as np
import logging
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
    pi = math.pi
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

    return list(map(math.degrees, [theta, phi, alpha, beta, gamma]))


def polar_to_cartesian(r, theta, phi):
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
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


def close_contact(mol_1, mol_2, factor):
    fragment_one, fragment_two = mol_1.coordinates, mol_2.coordinates
    radius_one, radius_two = mol_1.vdw_radius, mol_2.vdw_radius
    sum_of_radii = np.array([(x + y) for x, y in itertools.product(radius_one, radius_two)])
    inter_atomic_distance = np.array([np.linalg.norm(a - b) for a, b in itertools.product(fragment_one, fragment_two)])
    for l, m in zip(inter_atomic_distance, sum_of_radii):
        if l < m:
            return True
        else:
            print(l)
            return False


def gen_vectors(number_of_orientations):
    vectors = []
    a_threshold = math.pi / 2.0
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


def generate_orientations(molecule_id, seed, monomer, number_of_orientations):
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
        a_threshold = math.pi / 2.0
        d_threshold = 1.414
        for i in range(int(int(number_of_orientations) / 8)):
            for j in range(8):
                accepted = False
                tries = 1
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


def generate_orientations_from_given_vectors(molecule_id, seed, monomer, hm_orientations):
    noa = seed.number_of_atoms
    nob = monomer.number_of_atoms
    if noa == 1 and nob == 1:
        vectors = [[0.0, 0.0, 0.0, 0.0, 0.0]]
    else:
        vectors = gen_vectors(hm_orientations)
    orientations = {}
    filename_prefix = 'trial_'
    if seed.number_of_atoms > 1:
        seed.align()
    if monomer.number_of_atoms > 1:
        monomer.align()
    seed.move_to_origin()

    for i, each_vector in enumerate(vectors):
        orientation_id = "%03d_" % i + molecule_id
        orientation = merge_monomer_and_seed(each_vector, monomer, seed)

        orientation.title = 'trial orientation ' + orientation_id
        orientation.name = orientation_id

        orientation_xyz_file = filename_prefix + orientation_id + '.xyz'
        orientation.mol_to_xyz(orientation_xyz_file)
        orientations[orientation_id] = orientation

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


def test_gen_vec():
    x = gen_vectors(8 * 32)
    with open('tmp.xyz', 'w') as fp:
        fp.write(str(len(x)) + '\n\n')
        for i in x:
            x, y, z = polar_to_cartesian(5.0, i[0], i[1])
            fp.write('H ' + str(x) + ' ' + str(y) + ' ' + str(z) + '\n')
    return True


def main():
    pass


if __name__ == "__main__":
    main()
