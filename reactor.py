import copy
import os
import shutil
import sys

import numpy as np

import file_manager
import interface.babel
import tabu
from data_analysis import clustering
from optimiser import optimise
import logging
reactor_logger = logging.getLogger('pyar.reactor')

saved_products = {}
saved_inchi_strings = {}
saved_smile_strings = {}


def print_header(gamma_max, gamma_min, hm_orientations, software):
    reactor_logger.info('Starting PyAR 2.0\n')
    reactor_logger.info('{} orientations will be tried'.format(hm_orientations))
    reactor_logger.info(' Gamma (min): {}'.format(gamma_min))
    reactor_logger.info(' Gamma (max): {}'.format(gamma_max))
    reactor_logger.info(' Software   : {}'.format(software))


def react(reactant_a, reactant_b, gamma_min, gamma_max,
          hm_orientations, method,
          site, proximity_factor):
    """Run reactor module.  This is the outer loop
    generates all the orientations
    loop over all the gamma values
      optimize all orientations in each gamma
      after eliminating the products or failed geometries
    """
    cwd = os.getcwd()
    reactor_logger.debug('Current working directory: {}'.format(cwd))
    software = method['software']
    print_header(gamma_max, gamma_min, hm_orientations, software)
    # prepare job directories
    product_dir = cwd + '/products'
    reactor_logger.debug('Product directory: {}'.format(product_dir))
    file_manager.make_directories(product_dir)
    file_manager.make_directories('trial_geometries')
    os.chdir('trial_geometries')

    if site is None:
        all_orientations = tabu.generate_orientations('geom', reactant_a,
                                                      reactant_b,
                                                      hm_orientations)
    else:
        all_orientations = tabu.generate_guess_for_bonding('geom', reactant_a,
                                                           reactant_b,
                                                           site[0], site[1],
                                                           hm_orientations)

    os.chdir(cwd)

    gamma_list = np.linspace(gamma_min, gamma_max, num=10, dtype=float)
    orientations_to_optimize = all_orientations[:]

    for gamma in gamma_list:
        reactor_logger.info('  Current gamma : {}'.format(gamma))
        gamma_id = "%04d" % (int(gamma))
        gamma_home = cwd + '/gamma_' + gamma_id
        file_manager.make_directories(gamma_home)
        os.chdir(gamma_home)

        optimized_molecules = optimize_all(gamma_id, gamma,
                                           orientations_to_optimize,
                                           product_dir, method)

        reactor_logger.info("      {} geometries from this gamma cycle".format(len(optimized_molecules)))
        if len(optimized_molecules) == 0:
            reactor_logger.info("No orientations to be optimized for the next gamma cycle.")
            return
        if len(optimized_molecules) == 1:
            orientations_to_optimize = optimized_molecules[:]
        else:
            orientations_to_optimize = clustering.remove_similar(optimized_molecules)
        reactor_logger.info("Number of products found from gamma:{} = {}".format(gamma, len(saved_inchi_strings)))
        reactor_logger.info("{} geometries are considered for the next gamma cycle".format(len(orientations_to_optimize)))
        reactor_logger.debug("the keys of the molecules for next gamma cycle")
        for this_orientation in orientations_to_optimize:
            reactor_logger.debug("{}".format(this_orientation.name))
    os.chdir(cwd)
    return


def optimize_all(gamma_id, gamma, orientations_to_optimize,
                 product_dir, method):
    cwd = os.getcwd()
    table_of_optimized_molecules = []
    for this_molecule in orientations_to_optimize:
        job_key = this_molecule.name
        reactor_logger.info('   Orientation: {}'.format(job_key))
        o_key = "_{}".format(job_key[-8:])
        orientations_home = 'orientation' + o_key
        file_manager.make_directories(orientations_home)
        os.chdir(orientations_home)
        job_name = gamma_id + o_key
        this_molecule.name = job_name
        reactor_logger.info('Optimizing {}'.format(this_molecule.name))
        start_xyz_file_name = 'trial_' + this_molecule.name + '.xyz'
        this_molecule.mol_to_xyz(start_xyz_file_name)
        start_inchi = interface.babel.make_inchi_string_from_xyz(start_xyz_file_name)
        start_smile = interface.babel.make_smile_string_from_xyz(start_xyz_file_name)
        status = optimise(this_molecule, method, gamma=gamma)
        before_relax = copy.copy(this_molecule)
        this_molecule.name = job_name
        reactor_logger.info('... completed')
        if status is True or status == 'converged' or status == 'cycle_exceeded':
            reactor_logger.info("      E({}): {:12.7f}".format(job_name, this_molecule.energy))
            if this_molecule.is_bonded():
                reactor_logger.info("The fragments have close contracts. Going for relaxation")
                this_molecule.mol_to_xyz('trial_relax.xyz')
                this_molecule.name = 'relax'
                status = optimise(this_molecule, method)
                this_molecule.name = job_name
                if status is True or status == 'converged':
                    this_molecule.mol_to_xyz('result_relax.xyz')
                    current_inchi = interface.babel.make_inchi_string_from_xyz('result_relax.xyz')
                    current_smile = interface.babel.make_smile_string_from_xyz('result_relax.xyz')
                    reactor_logger.info('geometry relaxed')
                    reactor_logger.info("Checking for product formation with SMILE and InChi strings")
                    reactor_logger.info("Start SMILE: {} Current SMILE: {}".format(start_smile, current_smile))
                    reactor_logger.info("Start InChi: {} Current InChi: {}".format(start_inchi, current_inchi))

                    if start_inchi != current_inchi or start_smile != current_smile:
                        saved_products[job_name] = this_molecule
                        reactor_logger.info("       The geometry is different "
                                            "from the stating structure.")
                        reactor_logger.info("       Checking if this is a (new)"
                                            " products")
                        if current_inchi not in saved_inchi_strings.values() and \
                                current_smile not in saved_smile_strings.values():
                            reactor_logger.info("        New Product! Saving")
                            saved_inchi_strings[job_name] = current_inchi
                            saved_smile_strings[job_name] = current_smile
                            saved_products[job_name] = this_molecule
                            shutil.copy('result_relax.xyz',
                                        product_dir + '/' + job_name + '.xyz')
                            os.chdir(cwd)
                            continue
                        else:
                            reactor_logger.info("Both strings matches with "
                                                "those of already saved "
                                                "products. Discarded")
                            os.chdir(cwd)
                            continue
                    else:
                        table_of_optimized_molecules.append(before_relax)
                        reactor_logger.info('{} is added to the table to '
                                            'optimize with higher '
                                            'gamma'.format(job_name))
                elif status == 'cycle_exceeded':
                    table_of_optimized_molecules.append(before_relax)
                    reactor_logger.info('{} is added to the table to optimize '
                                        'with higher gamma'.format(job_name))
            else:
                table_of_optimized_molecules.append(this_molecule)
                reactor_logger.info('        no close contacts found')
                reactor_logger.info('       {} is added to the table to '
                                    'optimize with higher gamma'.format(job_name))
        os.chdir(cwd)
        sys.stdout.flush()
    return table_of_optimized_molecules


def main():
    pass


if __name__ == "__main__":
    main()
