"""
trial_generation.py

Functions to merge two molecules.
"""
import itertools
import logging

import numpy as np
from numpy import pi, cos, sin

from pyar.Molecule import Molecule
from pyar import orientation_sampling
# from pyar.property import get_connectivity
import networkx as nx
from typing import List, Tuple

trial_generation_logger = logging.getLogger('pyar.trial_generation')


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


def _minimum_contact_margin(mol_1, mol_2, factor):
    """Return the smallest clearance and vector from fragment two to one."""
    minimum_margin = float('inf')
    closest_vector = None
    for coordinate_one, radius_one in zip(mol_1.coordinates, mol_1.covalent_radius):
        for coordinate_two, radius_two in zip(mol_2.coordinates, mol_2.covalent_radius):
            vector = coordinate_one - coordinate_two
            separation = np.linalg.norm(vector)
            margin = separation - (radius_one + radius_two) * factor
            if margin < minimum_margin:
                minimum_margin = margin
                closest_vector = vector
    return minimum_margin, closest_vector


def _place_monomer_at_contact(seed, monomer, direction, distance_scaling):
    """Move a monomer to a non-overlapping contact position around a seed."""
    direction_norm = np.linalg.norm(direction)
    if direction_norm == 0.0:
        raise ValueError("trial placement direction cannot be a zero vector")
    direction = direction / direction_norm
    step_length = 0.05
    radial_step = direction * step_length

    maximum_contact_distance = max(
        (radius_one + radius_two) * distance_scaling
        for radius_one in seed.covalent_radius
        for radius_two in monomer.covalent_radius
    )
    seed_radius = max(np.linalg.norm(coordinate) for coordinate in seed.coordinates)
    monomer_radius = max(np.linalg.norm(coordinate) for coordinate in monomer.coordinates)
    outside_distance = seed_radius + monomer_radius + maximum_contact_distance + step_length

    # Retain the existing direction convention: placement begins along -direction.
    monomer.translate(direction * outside_distance)
    closest_margin = float('inf')
    closest_coordinates = monomer.coordinates.copy()
    maximum_radial_steps = int(np.ceil(2.0 * outside_distance / step_length)) + 2

    for _ in range(maximum_radial_steps):
        margin, _ = _minimum_contact_margin(seed, monomer, distance_scaling)
        if margin <= 0.0:
            monomer.translate(radial_step)
            return
        if margin < closest_margin:
            closest_margin = margin
            closest_coordinates = monomer.coordinates.copy()
        monomer.translate(-radial_step)

    # A radial line can pass through a gap in an irregular cluster. From its
    # closest point, approach the nearest atom pair to guarantee termination.
    monomer.coordinates = closest_coordinates
    margin, vector_to_seed = _minimum_contact_margin(seed, monomer, distance_scaling)
    if vector_to_seed is None or np.linalg.norm(vector_to_seed) == 0.0:
        raise RuntimeError("could not identify a contact direction for trial geometry")
    contact_step = vector_to_seed / np.linalg.norm(vector_to_seed) * step_length
    maximum_contact_steps = int(np.ceil(max(margin, 0.0) / step_length)) + 2
    trial_generation_logger.debug(
        "Radial placement missed the seed surface; using nearest-contact fallback."
    )
    for _ in range(maximum_contact_steps):
        if check_close_contact(seed, monomer, distance_scaling):
            monomer.translate(contact_step)
            return
        monomer.translate(-contact_step)
    raise RuntimeError("failed to place monomer at a non-overlapping contact position")


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
    direction = np.array([x, y, z], dtype=float)
    seed = seed_input.copy()
    monomer = monomer_input.copy()
    trial_generation_logger.debug('Merging two molecules')
    
    if not freeze_fragments:
        seed.move_to_origin()
    monomer.move_to_origin()

    if monomer.number_of_atoms > 1:
        monomer.rotate_3d((theta, phi, psi))
    
    trial_generation_logger.debug('Placing monomer at contact')
    _place_monomer_at_contact(seed, monomer, direction, distance_scaling)

    orientation = seed.merged_with(monomer)
    if site is not None:
        atoms_in_self = site[0]
        atoms_in_other = site[1]
    else:
        atoms_in_self = list(range(seed.number_of_atoms))
        atoms_in_other = list(range(seed.number_of_atoms, orientation.number_of_atoms))
    orientation.fragments = [atoms_in_self, atoms_in_other]
    trial_generation_logger.debug('Merged.')
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
    composite = seed.copy()
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

def generate_points(number_of_orientations, sequence_offset=0):
    """
    Generate trial directions and monomer orientations.

    :type number_of_orientations: int
    :param number_of_orientations: Number of orientations
    :param sequence_offset: Deterministic population index for distinct repeated samples
    :return: Numpy array of generated points in Cartesian coordinates and three angles
    :rtype: ndarray
    """
    return orientation_sampling.generate_trial_vectors(
        number_of_orientations,
        direction_method="fibonacci",
        rotation_method="halton",
        sequence_offset=sequence_offset,
        use_angles=True,
    )

def create_trial_geometries(molecule_id, seed, monomer,
                            number_of_orientations,
                            site):
    """
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
        proximity_factor = 1.2
    else:
        proximity_factor = 1.5

    trial_generation_logger.debug('Generating points')
    points_and_angles = orientation_sampling.generate_trial_vectors(
        number_of_orientations,
        direction_method="fibonacci",
        rotation_method="halton",
        use_angles=monomer.number_of_atoms > 1,
    )
    trial_generation_logger.debug('Generated points')
    write_trial_vectors(points_and_angles, 'trial_vectors.dat')
    # plot_points(pts)

    trial_generation_logger.debug('generate orientations from points and angles')

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


def write_trial_vectors(vectors, output_file):
    """
    Save generated trial direction and rotation vectors.

    :param vectors: Trial vectors to write.
    :type vectors: list or ndarray
    :param output_file: Output file name.
    :type output_file: str
    :rtype: None
    :returns: None
    """
    with open(output_file, 'a') as tf:
        for line in vectors:
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
