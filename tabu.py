import os.path
import random
import itertools

import numpy as np
from numpy import pi, arcsin, cos, sin, degrees
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
        import copy
        new_seed = copy.copy(seed)
        vectors = []
        a_threshold = pi / 2.0
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
    radius_one, radius_two = mol_1.vdw_radius, mol_2.vdw_radius
    status = False
    for i in range(len(fragment_one)):
        for j in range(len(fragment_two)):
            interatomic_distance = np.linalg.norm(fragment_one[i]-fragment_two[j])
            sum_of_radii = (radius_one[i] + radius_two[j])*factor
            if interatomic_distance < sum_of_radii:
                status = True
    return status


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


def rotating_octants(number_of_points, distance_tabu=True, angle_tabu=True):
    """
    Make N points ((x, y, z, theta, phi, psi) using 'rotating octants' method.

    :param number_of_points: int
    :param distance_tabu: float
    :param angle_tabu: float
    :return: numpy.array
    """
    saved_points_and_angles = []
    d_threshold = np.sqrt(2)/2.0
    a_threshold = pi/2.0
    for i in range(number_of_points//8):
        for j in itertools.product([1, -1], repeat=3):
            octant_chooser = np.array(j)
            if distance_tabu is False:
                saved_points_and_angles.append(make_point_n_angles(octant_chooser))
            if distance_tabu is True:
                if len(saved_points_and_angles) == 0:
                    saved_points_and_angles.append(make_point_n_angles(octant_chooser))
                else:
                    point_n_angle = make_an_untabooed_point(a_threshold, angle_tabu, d_threshold, octant_chooser, saved_points_and_angles)
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
    vec_2 = np.random.uniform(2*pi, size=(number_of_points, 3))
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
                v = p-q
                r = np.linalg.norm(v)
                f = v/r**3
                forces_on_p += f
            forces.append(forces_on_p)
        forces = np.array(forces)
        total_forces = np.sqrt(np.sum(forces**2))
        if total_forces > 0.25:
            fscale = 0.25 / total_forces
        else:
            fscale = 1
        dist = 0
        for i in range(len(points)):
            p = points[i]
            f = forces[i]
            moved_point = (p + f*fscale)
            moved_point /=  np.linalg.norm(moved_point)
            dist += np.linalg.norm(p-moved_point)
            points[i] = moved_point

        if dist < 1e-6:
            return np.concatenate((points, angles), axis=1)


def plot_points(pts):
    '''have to run with python -i '''
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt

    phi = np.linspace(0, np.pi, 20)
    theta = np.linspace(0, 2 * np.pi, 40)
    x = np.outer(np.sin(theta), np.cos(phi))
    y = np.outer(np.sin(theta), np.sin(phi))
    z = np.outer(np.cos(theta), np.ones_like(phi))

    fig, ax = plt.subplots(1, 1, subplot_kw={'projection':'3d', 'aspect':'equal'})
    ax.plot_wireframe(x, y, z, color='gray', rstride=1, cstride=1)
    ax.scatter(pts[:,0], pts[:,1], pts[:,2], s=100, c='r', zorder=10)
    # fig.show()
    fig.savefig('points.png')


def merge_two_molecules(vector, seed, monomer, freeze_fragments=False, site=None):
    x, y, z, theta, phi, psi = vector

    if freeze_fragments is False:
        seed.move_to_origin()
    monomer.move_to_origin()

    if monomer.number_of_atoms > 1:
        monomer.rotate_3d((theta, phi, psi))
    while close_contact(seed, monomer, 1.0):
        monomer.translate(np.array([x, y, z])/10)
    orientation = seed + monomer
    if site is not None:
        atoms_in_self = site
    else:
        atoms_in_self = [i for i in range(seed.number_of_atoms)]
    atoms_in_other = [i for i in range(seed.number_of_atoms, orientation.number_of_atoms)]
    orientation.fragments = [atoms_in_self, atoms_in_other]
    return orientation


def generate_orientations_from_points_and_angles(seed, monomer, points_and_angles):
    orientations = []
    for vector in points_and_angles:
        orientations.append(merge_two_molecules(vector, seed, monomer))
    return orientations


def generate_composite_molecule(seed, monomer, points_and_angles):
    import copy
    composite = copy.copy(seed)
    for vector in points_and_angles:
        composite = merge_two_molecules(vector, composite, monomer,
                                        freeze_fragments=True)
    return composite


def generate_orientations(molecule_id, seed, monomer, number_of_orientations,
                          method='rotating'):

    if method == 'rotating':
        if monomer.number_of_atoms ==1:
            pts = rotating_octants(number_of_orientations, angle_tabu=False)
        else:
            pts = rotating_octants(number_of_orientations, angle_tabu=True)
    elif method == 'random':
        pts = make_random_points_and_angles(number_of_orientations)
    elif method == 'uniform':
        pts = uniformly_distributed_points(number_of_orientations)
    else:
        tabu_logger.info('using default method: random')
        pts = make_random_points_and_angles(number_of_orientations)

    write_tabu_list(pts, 'tabu.dat')
    plot_points(pts)

    filename_prefix = 'trial_'
    orientations = generate_orientations_from_points_and_angles(seed, monomer, pts)
    for i, each_orientation in enumerate(orientations):
        each_orientation_id = "%03d_" % (i) + molecule_id
        each_orientation.title = 'trial orientation ' + each_orientation_id
        each_orientation.name = each_orientation_id
        each_orientation_xyz_file = filename_prefix + each_orientation_id + '.xyz'
        each_orientation.mol_to_xyz(each_orientation_xyz_file)

    return orientations


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-N', type=int, required=True,
                        help='number of points/configurations')
    parser.add_argument('-mpt', action='store_true',
                        help='generate random points on the surface of a unit '
                             'sphere')
    parser.add_argument('-mmc', action='store_true',
                        help='generate N configurations of two molecules')
    parser.add_argument('-mcm',  action='store_true',
                        help='generate a composite molecule with seed and N '
                             'monomers')
    parser.add_argument('-i', type=str, nargs=2,
                        dest='file', help="two molecules, seed and monomer, "
                                          "in xyz file format")
    parser.add_argument('-method',
                        choices=['random','rotating', 'uniform'],
                        help="method for generating points")
    parser.add_argument('-spr', type=str,
                        help="create a molecules from N atoms of given elememts")
    parser.add_argument('-best', type=int, nargs=2,
                        help="create the best orientation with lowest a b "
                             "distance")

    args = parser.parse_args()
    tabu_logger.debug(args)

    N = args.N

    if args.file:
        from Molecule import Molecule
        seed = Molecule.from_xyz(args.file[0])
        monomer = Molecule.from_xyz(args.file[1])
    else:
        seed = None
        monomer = None

    if args.method == 'rotating':
        if args.file and monomer.number_of_atoms ==1:
            pts = rotating_octants(N, angle_tabu=False)
        else:
            pts = rotating_octants(N, angle_tabu=True)
    elif args.method == 'random':
        pts = make_random_points_and_angles(N)
    elif args.method == 'uniform':
        pts = uniformly_distributed_points(N)
    else:
        print('using default method: random')
        pts = make_random_points_and_angles(N)

    if args.mpt:
        plot_points(pts)

    if args.mmc:
        orientations = generate_orientations_from_points_and_angles(seed, monomer, pts)
        for i, each_orientation in enumerate(orientations):
            each_orientation.mol_to_xyz('mol'+str(i)+'.xyz')

    if args.mcm:
        result = generate_composite_molecule(seed, monomer, pts)
        result.mol_to_xyz('mol.xyz')

    if args.spr:
        import scipy.spatial.distance as sdist

        from Molecule import Molecule
        from atomic_data import atomic_numbers, covalent_radii
        radii = covalent_radii[atomic_numbers[args.spr.capitalize()]]
        pts = pts[:, :3]
        while np.min(sdist.pdist(pts)) < 2*radii:
            pts += pts / 10

        atoms = [args.spr for _ in range(N)]
        result = Molecule(atoms, pts)
        result.mol_to_xyz('mol.xyz')
    if args.best:
        a = args.best[0]
        b = seed.number_of_atoms+args.best[1]
        orientations = generate_orientations_from_points_and_angles(seed, monomer, pts)
        dictor = {}
        for i, each_orientation in enumerate(orientations):
            coords = each_orientation.coordinates
            dist = np.linalg.norm(coords[a]-coords[b])
            dictor[i] = dist

        best = min(dictor, key=dictor.get)
        print(best, dictor[best])
        orientations[best].mol_to_xyz('mol.xyz')
        import optimiser
        method_here = {'software':'xtb_turbo',
                       method={'software':'xtb_turbo', 'charge':0,
                               'scftype':'rhf', 'multiplicity':1}
        optimiser.optimise(orientations[best], method=method_here , gamma=1000)

if __name__ == "__main__":
    main()
