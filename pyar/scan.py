import os
import time

import numpy as np

from pyar import optimiser
from pyar.tabu import rotating_octant, tabu_logger, generate_orientations_from_points_and_angles, write_tabu_list


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

    return orientations


def scan_distance(input_molecules, site_atoms, number_of_orientations, quantum_chemistry_parameters):
    a_molecule = input_molecules[0]
    b_molecule = input_molecules[1]
    a_atom = site_atoms[0]
    b_atom = a_molecule.number_of_atoms + a_atom
    input_molecules = generate_guess_for_bonding('abc',
                                                 a_molecule, b_molecule,
                                                 a_atom, b_atom,
                                                 int(number_of_orientations))
    for each_molecule in input_molecules:
        coordinates = each_molecule.coordinates
        start_dist = np.linalg.norm(coordinates[a_atom] - coordinates[b_atom])
        final_distance = each_molecule.covalent_radius[a_atom] + each_molecule.covalent_radius[b_atom]
        step = int(abs(final_distance - start_dist) * 10)
        if quantum_chemistry_parameters['software'] == 'orca':
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
            logger.error('Optimization with %s is not implemented '
                         'yet' % quantum_chemistry_parameters['software'])
