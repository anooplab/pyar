import os
import time

import numpy as np

import pyar.tabu
from pyar import old_optimiser
from pyar.data_analysis import clustering


def generate_guess_for_bonding(molecule_id, seed, monomer, a, b,
                               number_of_orientations, d_scale):
    # pts = pyar.tabu.generate_points(1, True, True, True, 0.3, 5.0)
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
        each_orientation = pyar.tabu.merge_two_molecules(x.x, seed, monomer,
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
    orientation = pyar.tabu.merge_two_molecules(pts, seed, monomer, site=[a, b])
    coordinates = orientation.coordinates
    return np.linalg.norm(coordinates[a] - coordinates[b])


def generate_guess_for_bonding_brute_force(molecule_id, seed, monomer, a, b, number_of_orientations, d_scale):
    tabu_check_for_angles = monomer.number_of_atoms != 1
    saved_pts = []
    orientations = []
    for _ in range(number_of_orientations):
        t1 = time.time()
        pts = pyar.tabu.generate_points(32, True, True, True, 0.3, 5.0)
        t2 = time.time()
        pyar.tabu.tabu_logger.debug(f'Created points: in {t2 - t1} seconds')
        t1 = time.time()
        current_orientations = [pyar.tabu.merge_two_molecules(vector, seed, monomer, site=[a, b]) for vector in pts]

        t2 = time.time()
        pyar.tabu.tabu_logger.debug(f'Created orientations {t2 - t1} seconds')
        t1 = time.time()
        stored_orientations = {}
        for j, each_orientation in enumerate(current_orientations):
            coords = each_orientation.coordinates
            dist = np.linalg.norm(coords[a] - coords[b])
            stored_orientations[j] = dist
        best_orientation = min(stored_orientations, key=stored_orientations.get)
        best_point = pts[best_orientation]
        pyar.tabu.tabu_logger.debug(f"{best_orientation} {stored_orientations[best_orientation]}")


        saved_pts.append(best_point)
        orientations.append(current_orientations[best_orientation])
        t2 = time.time()
        pyar.tabu.tabu_logger.debug(f'Found best orientation in {t2 - t1} seconds')

    t1 = time.time()
    filename_prefix = 'trial_'
    for i, each_orientation in enumerate(orientations):
        each_orientation_id = f"{i:03d}_{molecule_id}"
        each_orientation.title = f'trial orientation {each_orientation_id}'
        each_orientation.name = each_orientation_id
        each_orientation_xyz_file = filename_prefix + each_orientation_id + '.xyz'
        each_orientation.mol_to_xyz(each_orientation_xyz_file)
    t2 = time.time()
    pyar.tabu.tabu_logger.debug(f'Wrote files in {t2 - t1} seconds')
    pyar.tabu.write_tabu_list(saved_pts, 'tabu.dat')
    return orientations


def scan_distance(input_molecules, site_atoms, number_of_orientations,
                  quantum_chemistry_parameters):
    a_molecule = input_molecules[0]
    b_molecule = input_molecules[1]
    a_atom = site_atoms[0]
    b_atom = a_molecule.number_of_atoms + a_atom -  1
    proximity_factor = 1.5  #
    input_molecules = generate_guess_for_bonding('abc',
                                                 a_molecule, b_molecule,
                                                 a_atom, b_atom,
                                                 int(number_of_orientations),
                                                 d_scale=proximity_factor)

    for each_molecule in input_molecules:
        coordinates = each_molecule.coordinates
        start_dist = np.linalg.norm(coordinates[a_atom] - coordinates[b_atom])
        final_distance = each_molecule.covalent_radius[a_atom] + \
                         each_molecule.covalent_radius[b_atom]
        if quantum_chemistry_parameters['software'] == 'orca':
            step = int(abs(final_distance - start_dist) * 10)
            c_k = f'\n% geom\n    scan B {a_atom} {b_atom} = {start_dist}, ' \
                  f'{final_distance}, {step}\n        end\nend\n'
            quantum_chemistry_parameters['gamma'] = 0.0
            cwd = os.getcwd()
            job_dir = 'scans'
            from pyar import file_manager
            file_manager.make_directories(job_dir)
            os.chdir(job_dir)
            quantum_chemistry_parameters['custom_keyword'] = c_k
            old_optimiser.optimise(each_molecule, quantum_chemistry_parameters)
            os.chdir(cwd)
        else:
            print('Optimization with %s is not implemented '
                  'yet' % quantum_chemistry_parameters['software'])
