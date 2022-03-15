import os
import time

import numpy as np

import pyar.tabu
from pyar import optimiser


def generate_guess_for_bonding(molecule_id, seed, monomer, a, b,
                               number_of_orientations, d_scale):
    tabu_check_for_angles = monomer.number_of_atoms != 1
    saved_pts = []

    orientations = []
    for _ in range(number_of_orientations):
        t1 = time.time()
        pts = pyar.tabu.generate_points(32, True, True, True, 0.3, 5.0)
        t2 = time.time()
        pyar.tabu.tabu_logger.debug('Created points: in {} seconds'.format(t2 - t1))
        t1 = time.time()
        current_orientations = [
            pyar.tabu.merge_two_molecules(vector, seed, monomer, site=[a, b])
            for vector in pts
        ]

        t2 = time.time()
        pyar.tabu.tabu_logger.debug('Created orientations {} seconds'.format(t2 - t1))
        t1 = time.time()
        stored_orientations = {}
        for j, each_orientation in enumerate(current_orientations):
            coords = each_orientation.coordinates
            dist = np.linalg.norm(coords[a] - coords[b])
            stored_orientations[j] = dist
        best_orientation = min(stored_orientations, key=stored_orientations.get)
        best_point = pts[best_orientation]
        pyar.tabu.tabu_logger.debug("{} {}".format(best_orientation, stored_orientations[best_orientation]))
        saved_pts.append(best_point)
        orientations.append(current_orientations[best_orientation])
        t2 = time.time()
        pyar.tabu.tabu_logger.debug('Found best orientation in {} seconds'.format(t2 - t1))

    t1 = time.time()
    filename_prefix = 'trial_'
    for i, each_orientation in enumerate(orientations):
        each_orientation_id = f"{i:03d}_{molecule_id}_"
        each_orientation.title = f'trial orientation {each_orientation_id}'
        each_orientation.name = each_orientation_id
        each_orientation_xyz_file = filename_prefix + each_orientation_id + '.xyz'
        each_orientation.mol_to_xyz(each_orientation_xyz_file)
    t2 = time.time()
    pyar.tabu.tabu_logger.debug('Wrote files in {} seconds'.format(t2 - t1))
    pyar.tabu.write_tabu_list(saved_pts, 'tabu.dat')

    return orientations


def scan_distance(input_molecules, site_atoms, number_of_orientations, quantum_chemistry_parameters):
    a_molecule = input_molecules[0]
    b_molecule = input_molecules[1]
    a_atom = site_atoms[0]
    b_atom = a_molecule.number_of_atoms + a_atom
    proximity_factor = 2.3  # TODO: find the optimum value
    input_molecules = generate_guess_for_bonding('abc',
                                                 a_molecule, b_molecule,
                                                 a_atom, b_atom,
                                                 int(number_of_orientations),
                                                 d_scale=proximity_factor)
    for each_molecule in input_molecules:
        coordinates = each_molecule.coordinates
        start_dist = np.linalg.norm(coordinates[a_atom] - coordinates[b_atom])
        final_distance = each_molecule.covalent_radius[a_atom] + each_molecule.covalent_radius[b_atom]
        if quantum_chemistry_parameters['software'] == 'orca':
            step = int(abs(final_distance - start_dist) * 10)
            c_k = f'\n% geom\n    scan B {a_atom} {b_atom}= {start_dist}, ' \
                  f'{final_distance}, {step}\n        end\nend\n'
            quantum_chemistry_parameters['gamma'] = 0.0
            cwd = os.getcwd()
            job_dir = 'scans'
            from pyar import file_manager
            file_manager.make_directories(job_dir)
            os.chdir(job_dir)
            quantum_chemistry_parameters['custom_keywords'] = c_k
            optimiser.optimise(each_molecule, quantum_chemistry_parameters)
            os.chdir(cwd)
        else:
            print('Optimization with %s is not implemented '
                  'yet' % quantum_chemistry_parameters['software'])
