from __future__ import print_function

import math
import os.path
import random
import sys

import numpy as np

from atomic_data import atomic_numbers, vdw_radii, covalent_radii


def load_tabu_list(tabu_file='tabu.dat'):
    if not os.path.exists(tabu_file):
        return []
    with open(tabu_file, 'r') as tf:
        tabu_list = [map(float, i.split()) for i in tf.readlines()]
        return tabu_list


def write_tabu_list(tabu_list, tabu_file):
    with open(tabu_file, 'w') as tf:
        for i in tabu_list:
            for j in i:
                tf.write(str(j) + ' ')
            tf.write('\n')


def gen_a_set_of_angles(which_octant):
    pi = math.pi
    a = 0.0
    b = pi / 2.0
    c = pi
    d = 3 * pi / 2.0
    e = 2 * pi

    octant = {
        1: [a, b, a, b], 2: [b, c, a, b], 3: [a, b, b, c], 4: [b, c, b, c], \
        5: [a, b, c, d], 6: [b, c, c, d], 7: [a, b, d, e], 8: [b, c, d, e] \
        }

    i, j, k, l = octant[which_octant]

    theta = random.uniform(i, j)
    phi = random.uniform(k, l)
    alpha = random.uniform(a, e)
    beta = random.uniform(a, c)
    gamma = random.uniform(a, e)

    return theta, phi, alpha, beta, gamma


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
        previous_apha, previous_beta, previous_gamma = i[2:]
        previous_point = np.array(polar_to_cartesian(r, previous_theta, previous_phi))
        distance = np.linalg.norm(current_point - previous_point)
        diff_alpha = abs(previous_apha - alpha)
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


def close_contact(mol_1, mol_2, contact_type):
    if contact_type == 'vdw':
        mol_1.radii = [vdw_radii[atomic_numbers[c]] for c in mol_1.atoms_list]
        mol_2.radii = [vdw_radii[atomic_numbers[c]] for c in mol_2.atoms_list]
    elif contact_type == 'covalent':
        mol_1.radii = [covalent_radii[atomic_numbers[c]] for c in mol_1.atoms_list]
        mol_2.radii = [covalent_radii[atomic_numbers[c]] for c in mol_2.atoms_list]

    for i in range(mol_1.number_of_atoms):
        for j in range(mol_2.number_of_atoms):
            R_AB = mol_1.radii[i] + mol_2.radii[j]
            distance = np.linalg.norm(mol_1.coordinates[i] - mol_2.coordinates[j])
            if distance <= R_AB:
                return True
    else:
        return False


def generate_orientations(molecule_id, seed, monomer, hm_orientations):
    noa = seed.number_of_atoms
    nob = monomer.number_of_atoms
    if noa == 1 and nob == 1:
        vectors = [[0.0, 0.0, 0.0, 0.0, 0.0]]
    else:
        vectors = gen_vectors(hm_orientations)
    orientations = {}
    tabu_list = []
    filename_prefix = 'trial_'
    if seed.number_of_atoms > 1:
        seed.align()
    if monomer.number_of_atoms > 1:
        monomer.align()
    seed.move_to_origin()

    for i, each_vector in enumerate(vectors):
        orientation_id = "%03d_" % i + molecule_id
        theta, phi, alpha, beta, gamma = each_vector
        monomer.move_to_origin()
        if monomer.number_of_atoms > 1:
            monomer.rotate_3d((alpha, beta, gamma))
        r = 0.1
        while close_contact(seed, monomer, 'vdw') is True:
            r += 0.1
            if seed.number_of_atoms == 2: phi = 0.0
            monomer.translate(polar_to_cartesian(r, theta, phi))

        orientation = seed + monomer
        orientation.title = 'trial orientation ' + orientation_id
        orientation.name = orientation_id
        tabu_list.append([r, theta, phi, alpha, beta, gamma])
        orientation_xyz_file = filename_prefix + orientation_id + '.xyz'
        orientation.mol_to_xyz(orientation_xyz_file)
        orientations[orientation_id] = orientation

    write_tabu_list(tabu_list, 'tabu.dat')
    return orientations


def test_gen_vec():
    x = gen_vectors(8 * 32)
    with open('tmp.xyz', 'w') as fp:
        fp.write(str(len(x)) + '\n\n')
        for i in x:
            x, y, z = polar_to_cartesian(5.0, i[0], i[1])
            fp.write('H ' + str(x) + ' ' + str(y) + ' ' + str(z) + '\n')
    return True


if __name__ == "__main__":
    number_of_atoms, title, atoms_list, coords_1 = xyz_file.read(sys.argv[1])
    coords_2 = coords_1[:]
    x = gen_vectors(256)
    generate_orientations(coords_1, coords_2, x)
