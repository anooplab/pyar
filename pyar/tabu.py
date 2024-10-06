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
from scipy.spatial.distance import cosine
from scipy.stats import qmc

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
    Checks for close contacts between two molecules.

    :param mol_1: First molecule
    :param mol_2: Second molecule
    :param factor: Scaling factor for van der Waals radii
    :return: boolean indicating close contact
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


def merge_two_molecules(vector, seed_input, monomer_input, freeze_fragments=False, site=None, distance_scaling=1.5):
    """
    Merge two molecules with improved proximity control.

    :param vector: the direction for placing the monomer
    :param seed_input: Seed molecule
    :param monomer_input: Monomer molecule
    :param freeze_fragments: Whether to rotate the monomer
    :param site: weighted atom centres for placing
    :param distance_scaling: minimum separation between fragments
    :return: merged molecule
    """
    x, y, z, theta, phi, psi = vector
    direction = np.array([x, y, z])
    tiny_steps = direction / 20  # Reduced step size for finer control
    seed = copy.deepcopy(seed_input)
    monomer = copy.deepcopy(monomer_input)
    tabu_logger.debug('Merging two molecules')
    
    if not freeze_fragments:
        seed.move_to_origin()
    monomer.move_to_origin()

    if monomer.number_of_atoms > 1:
        monomer.rotate_3d((theta, phi, psi))
    
    tabu_logger.debug('Checking close contact')

    # Calculate initial placement distance
    r_max_of_seed = np.max([np.linalg.norm(c) for c in seed.coordinates])
    r_max_of_monomer = np.max([np.linalg.norm(c) for c in monomer.coordinates])
    initial_distance = r_max_of_seed + r_max_of_monomer + 0.5  # Reduced initial distance
    
    move_to = direction * initial_distance
    monomer.translate(move_to)

    # Gradually move monomer closer to seed
    while not check_close_contact(seed, monomer, distance_scaling):
        monomer.translate(-1 * tiny_steps)
    
    # Fine-tune the position
    monomer.translate(tiny_steps)

    orientation = seed + monomer
    if site is not None:
        atoms_in_self = site[0]
        atoms_in_other = site[1]
    else:
        atoms_in_self = list(range(seed.number_of_atoms))
        atoms_in_other = list(range(seed.number_of_atoms, orientation.number_of_atoms))
    orientation.fragments = [atoms_in_self, atoms_in_other]
    tabu_logger.debug('Merged.')
    return orientation


def ellipsoid_wall_potential(coordinates, a, b, c, k=100.0):
    """
    Apply an ellipsoid wall potential to keep molecules within bounds.

    :param coordinates: numpy array of shape (n_atoms, 3) with atomic coordinates
    :param a: semi-major axis in x direction
    :param b: semi-major axis in y direction
    :param c: semi-major axis in z direction
    :param k: force constant for the wall potential
    :return: potential energy from the wall
    """
    x, y, z = coordinates.T
    r = np.sqrt((x/a)**2 + (y/b)**2 + (z/c)**2)
    return np.sum(k * np.maximum(r - 1, 0)**2)


def create_composite_molecule(seed, monomer, points_and_angles, d_scale):
    """
    Create a composite molecule with ellipsoid wall potential.

    :param seed: seed molecule
    :param monomer: monomer molecule
    :param points_and_angles: a set of 3 points and 3 angles
    :param d_scale: distance scaling factor
    :return: a composite molecule
    """
    composite = copy.deepcopy(seed)
    a, b, c = 1.0, 1.0, 1.0  # Initial ellipsoid parameters
    
    for vector in points_and_angles:
        composite = merge_two_molecules(vector, composite, monomer,
                                        freeze_fragments=False,
                                        distance_scaling=d_scale)
        
        # Update ellipsoid parameters based on current molecule size
        max_coords = np.max(composite.coordinates, axis=0)
        min_coords = np.min(composite.coordinates, axis=0)
        a, b, c = (max_coords - min_coords) / 2 + 1.0  # Add buffer
        
        # Apply wall potential
        ellipsoid_wall_potential(composite.coordinates, a, b, c)
    
    return composite


def spherical_to_cartesian(r, theta, phi):
    """
    Convert spherical coordinates to Cartesian coordinates.

    :param r: Radius
    :param theta: Polar angle
    :param phi: Azimuthal angle
    :return: Cartesian coordinates (x, y, z)
    """
    x = r * sin(theta) * cos(phi)
    y = r * sin(theta) * sin(phi)
    z = r * cos(theta)
    return np.array([x, y, z])

def check_tabu_status(point_n_angle, d_threshold, tabu_list):
    """
    Check if the given point (point_n_angle) is within the proximity (set by
    d_threshold) of any point in the tabu list using cosine similarity.

    :type point_n_angle: ndarray
    :param point_n_angle: Five numbers (theta, phi, angle1, angle2, angle3)
    :type d_threshold: float
    :param d_threshold: Minimum cosine similarity for Tabu check
    :type tabu_list: list
    :param tabu_list: The saved points and angles
    :rtype: bool
    :return: True if the new point is within the proximity of saved points
    """
    for each_saved_entry in tabu_list:
        if 1 - cosine(each_saved_entry, point_n_angle) > d_threshold:
            return True
    return False

def generate_points(number_of_orientations, tabu_on, d_threshold=0.95):
    """
    Generate points using LHS for both spherical coordinates and angles, and convert to Cartesian coordinates.

    :type number_of_orientations: int
    :param number_of_orientations: Number of orientations
    :type tabu_on: bool
    :param tabu_on: Use tabu or not
    :type d_threshold: float
    :param d_threshold: Cosine similarity threshold for Tabu check
    :return: Numpy array of generated points in Cartesian coordinates and three angles
    :rtype: ndarray
    """
    # Create a Latin Hypercube sampler for 5 dimensions (theta, phi, angle1, angle2, angle3)
    sampler = qmc.LatinHypercube(d=5)
    samples = sampler.random(n=number_of_orientations * 100)  # Generate more samples initially

    # Scale the first dimension to range [0, pi] for theta
    samples[:, 0] = np.arccos(1 - 2 * samples[:, 0])  # theta
    # Scale the second dimension to range [0, 2*pi] for phi
    samples[:, 1] = samples[:, 1] * 2 * pi  # phi

    # Scale the last three dimensions to range [0, 2*pi] for angles
    samples[:, 2] = samples[:, 2] * 2 * pi  # angle1
    samples[:, 3] = samples[:, 3] * 2 * pi  # angle2
    samples[:, 4] = samples[:, 4] * 2 * pi  # angle3

    tabu_list = []
    valid_samples = 0
    i = 0

    while valid_samples < number_of_orientations and i < len(samples):
        point_n_angle = samples[i]
        if tabu_on:
            if not check_tabu_status(point_n_angle, d_threshold, tabu_list):
                tabu_list.append(point_n_angle)
                valid_samples += 1
        else:
            tabu_list.append(point_n_angle)
            valid_samples += 1
        i += 1

    if valid_samples < number_of_orientations:
        raise ValueError("Unable to generate enough valid points with the given constraints.")

    # Convert to Cartesian coordinates with three angles
    result_points = []
    for point in tabu_list:
        theta, phi, angle1, angle2, angle3 = point
        cartesian_coords = spherical_to_cartesian(1, theta, phi)
        result_points.append(np.concatenate((cartesian_coords, [angle1, angle2, angle3])))

    return np.array(result_points)

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
    points_and_angles = generate_points(number_of_orientations, tabu_on)
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




