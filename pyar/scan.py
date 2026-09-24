import time

import numpy as np

import pyar.sampling.trial_generator as trial_generation
from pyar.selection import clustering


def generate_guess_for_bonding(molecule_id, seed, monomer, a, b,
                               number_of_orientations, d_scale):
    from scipy.optimize import differential_evolution as global_opt
    from functools import partial
    my_bounds = [(-0.5, 0.5), (-0.5, 0.5), (-0.5, 0.5), (0, 2 * np.pi),
                 (0, 2 * np.pi), (0, 2 * np.pi)]
    orientations = []
    fun = partial(ab_dist, a, b, monomer, seed)
    for i in range(number_of_orientations):
        x = global_opt(fun, my_bounds,
                       polish=True, disp=True, workers=-1)
        print(x.message)
        filename_prefix = "aai_"
        each_orientation = trial_generation.merge_two_molecules(x.x, seed, monomer,
                                                         site=[a, b],
                                                         distance_scaling=d_scale)
        each_orientation_id = f"{i:03d}_{molecule_id}"
        each_orientation.title = f'trial orientation {each_orientation_id}'
        each_orientation.name = each_orientation_id
        each_orientation.energy = 0.0
        each_orientation_xyz_file = filename_prefix + each_orientation_id + '.xyz'
        each_orientation.mol_to_xyz(each_orientation_xyz_file)
        orientations.append(each_orientation)
    try:
        return clustering.remove_similar(orientations)
    except Exception:
        return orientations


def ab_dist(a, b, monomer, seed, pts):
    orientation = trial_generation.merge_two_molecules(pts, seed, monomer, site=[a, b])
    coordinates = orientation.coordinates
    return np.linalg.norm(coordinates[a] - coordinates[b])


def generate_guess_for_bonding_brute_force(molecule_id, seed, monomer, a, b, number_of_orientations, d_scale):
    saved_pts = []
    orientations = []
    for population_index in range(number_of_orientations):
        t1 = time.time()
        pts = trial_generation.generate_points(32, sequence_offset=population_index)
        t2 = time.time()
        trial_generation.trial_generation_logger.debug(f'Created points: in {t2 - t1} seconds')
        t1 = time.time()
        current_orientations = [trial_generation.merge_two_molecules(vector, seed, monomer, site=[a, b]) for vector in pts]

        t2 = time.time()
        trial_generation.trial_generation_logger.debug(f'Created orientations {t2 - t1} seconds')
        t1 = time.time()
        stored_orientations = {}
        for j, each_orientation in enumerate(current_orientations):
            coords = each_orientation.coordinates
            dist = np.linalg.norm(coords[a] - coords[b])
            stored_orientations[j] = dist
        best_orientation = min(stored_orientations, key=stored_orientations.get)
        best_point = pts[best_orientation]
        trial_generation.trial_generation_logger.debug(f"{best_orientation} {stored_orientations[best_orientation]}")


        saved_pts.append(best_point)
        orientations.append(current_orientations[best_orientation])
        t2 = time.time()
        trial_generation.trial_generation_logger.debug(f'Found best orientation in {t2 - t1} seconds')

    t1 = time.time()
    filename_prefix = 'trial_'
    for i, each_orientation in enumerate(orientations):
        each_orientation_id = f"{i:03d}_{molecule_id}"
        each_orientation.title = f'trial orientation {each_orientation_id}'
        each_orientation.name = each_orientation_id
        each_orientation_xyz_file = filename_prefix + each_orientation_id + '.xyz'
        each_orientation.mol_to_xyz(each_orientation_xyz_file)
    t2 = time.time()
    trial_generation.trial_generation_logger.debug(f'Wrote files in {t2 - t1} seconds')
    trial_generation.write_trial_vectors(saved_pts, 'trial_vectors.dat')
    return orientations
