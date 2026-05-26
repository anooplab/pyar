"""Trial geometry generation and merge helpers."""

from __future__ import annotations

import itertools
import logging
from typing import List, Tuple

import networkx as nx
import numpy as np
from numpy import pi, cos, sin

from pyar.core.molecule import Molecule
from pyar.sampling.rotation import (
    canonicalize_quaternion,
    halton_quaternions,
    quaternions_to_euler_zxz,
)
from pyar.sampling.sphere import generate_directions

trial_generation_logger = logging.getLogger('pyar.trial_generation')


def polar_to_cartesian(r, theta, phi):
    """Convert spherical polar coordinates to Cartesian coordinates."""
    x = r * sin(theta) * cos(phi)
    y = r * sin(theta) * sin(phi)
    z = r * cos(theta)
    return x, y, z


def check_close_contact(mol_1, mol_2, factor):
    """Return whether two fragments contain a scaled close contact."""
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
    """Return the minimum inter-fragment Cartesian separation."""
    return np.min([np.linalg.norm(a - b) for a, b in itertools.product(mol_1.coordinates, mol_2.coordinates)])


def merge_two_molecules(vector, seed_input, monomer_input, freeze_fragments=False, site=None, distance_scaling=1.5):
    """Place and merge a monomer around a seed using one trial vector."""
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
    """Return an ellipsoidal wall penalty for coordinates outside the bounds."""
    x, y, z = coordinates.T
    r = np.sqrt((x/a)**2 + (y/b)**2 + (z/c)**2)
    return np.sum(k * np.maximum(r - 1, 0)**2)


def create_composite_molecule(seed, monomer, points_and_angles, d_scale):
    """Create a composite molecule by applying consecutive placements."""
    composite = seed.copy()
    a, b, c = 1.0, 1.0, 1.0
    for vector in points_and_angles:
        composite = merge_two_molecules(vector, composite, monomer, freeze_fragments=False, distance_scaling=d_scale)
        max_coords = np.max(composite.coordinates, axis=0)
        min_coords = np.min(composite.coordinates, axis=0)
        a, b, c = (max_coords - min_coords) / 2 + 1.0
        ellipsoid_wall_potential(composite.coordinates, a, b, c)
    return composite


def spherical_to_cartesian(r, theta, phi):
    """Convert spherical coordinates to a Cartesian vector."""
    x = r * sin(theta) * cos(phi)
    y = r * sin(theta) * sin(phi)
    z = r * cos(theta)
    return np.array([x, y, z])


def generate_trial_vectors(
    number_of_points: int,
    *,
    direction_method: str = "fibonacci",
    rotation_method: str = "halton",
    seed: int | None = None,
    sequence_offset: int = 0,
    oversample_factor: int = 8,
    use_angles: bool = True,
) -> np.ndarray:
    """Return direction and orientation vectors for trial geometries."""
    count = int(number_of_points)
    if count < 1:
        raise ValueError("number_of_points must be positive")
    directions = generate_directions(
        direction_method,
        count,
        seed=seed,
        sequence_offset=sequence_offset,
        oversample_factor=oversample_factor,
    )
    if use_angles:
        rotation_offset = None
        if sequence_offset or seed is not None:
            rotation_offset = int(sequence_offset) * count
            if seed is not None:
                rotation_offset += int(seed)
        if rotation_method == "halton":
            quaternions = halton_quaternions(count, seed=rotation_offset)
        elif rotation_method == "random":
            from pyar.sampling.rotation import random_quaternions
            quaternions = random_quaternions(count, seed=rotation_offset)
        else:
            raise ValueError(f"Unknown rotation method: {rotation_method}")
    else:
        quaternions = None
    if use_angles:
        canonical_quaternions = np.array([canonicalize_quaternion(q) for q in quaternions], dtype=float)
        angles = quaternions_to_euler_zxz(canonical_quaternions)
    else:
        angles = np.zeros((len(directions), 3), dtype=float)
    return np.column_stack((directions, angles))


def generate_points(number_of_orientations, sequence_offset=0):
    """Generate default trial direction and monomer-orientation vectors."""
    return generate_trial_vectors(
        number_of_orientations,
        direction_method="fibonacci",
        rotation_method="halton",
        sequence_offset=sequence_offset,
        use_angles=True,
    )


def create_trial_geometries(molecule_id, seed, monomer, number_of_orientations, site):
    """Generate and write merged trial geometries for one growth operation."""
    if monomer.number_of_atoms == 1:
        proximity_factor = 1.2
    else:
        proximity_factor = 1.5

    trial_generation_logger.debug('Generating points')
    points_and_angles = generate_trial_vectors(
        number_of_orientations,
        direction_method="fibonacci",
        rotation_method="halton",
        use_angles=monomer.number_of_atoms > 1,
    )
    trial_generation_logger.debug('Generated points')
    write_trial_vectors(points_and_angles, 'trial_vectors.dat')
    trial_generation_logger.debug('generate orientations from points and angles')

    orientations = []
    filename_prefix = 'trial_'
    for counter, vector in enumerate(points_and_angles):
        new_orientation = merge_two_molecules(vector, seed, monomer, site=site, distance_scaling=proximity_factor)
        new_orientation_id = f"{counter:03d}_{molecule_id}"
        new_orientation.title = f'trial orientation {new_orientation_id}'
        new_orientation.name = new_orientation_id
        new_orientation_xyz_file = f'{filename_prefix}{new_orientation_id}.xyz'
        new_orientation.mol_to_xyz(new_orientation_xyz_file)
        orientations.append(new_orientation)

    return orientations


def plot_points(points, location):
    """Plot sampled direction points on a unit sphere."""
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


def write_trial_vectors(vectors, output_file):
    """Append generated direction and rotation vectors to an output file."""
    with open(output_file, 'a') as tf:
        for line in vectors:
            tf.write(f'{str(line)} ')
        tf.write('\n')


def main():
    pass


def get_connectivity(coordinates: np.ndarray, covalent_radii: List[float], tolerance: float = 0.4) -> List[Tuple[int, int]]:
    """Return covalent-radius-based connectivity edges for a geometry."""
    n_atoms = len(coordinates)
    bonds = []
    for i in range(n_atoms):
        for j in range(i + 1, n_atoms):
            distance = np.linalg.norm(coordinates[i] - coordinates[j])
            if distance < (covalent_radii[i] + covalent_radii[j]) * (1 + tolerance):
                bonds.append((i, j))
    return bonds


def broken(molobj) -> bool:
    """Return whether a molecule is disconnected under its covalent-radius graph."""
    graph = nx.Graph()
    graph.add_nodes_from(range(len(molobj.atoms_list)))
    graph.add_edges_from(get_connectivity(molobj.coordinates, molobj.covalent_radius))
    return not nx.is_connected(graph)


if __name__ == "__main__":
    import sys

    molecule = Molecule.from_xyz(sys.argv[1])
    print("is a fragment" if broken(molecule) else "Fine")
