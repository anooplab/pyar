import copy
import logging
import os
import shutil
import sys

import numpy as np
from pyar.checkpt import dumpchk, readchk, updtchk
import pyar.interface.babel
import pyar.scan
from pyar import tabu, file_manager
from pyar.data_analysis import clustering
from pyar.optimiser import optimise

reactor_logger = logging.getLogger('pyar.reactor')

saved_products = {}
saved_inchi_strings = {}
saved_smile_strings = {}


def print_header(gamma_max, gamma_min, hm_orientations, software):
    pass


def react(reactant_a, reactant_b, gamma_min, gamma_max, hm_orientations, qc_params,
          site, proximity_factor, tabu_on=None, grid_on=None):
    """
    The Reactor module

    This is the outer loop generates all the orientations
    loop over all the gamma values optimize all orientations
    in each gamma after eliminating the products or failed geometries.

    """
    global workdir
    workdir = os.getcwd()

    if readchk(workdir) is not None:
        chk = readchk(workdir)
        import shutil
        # shutil.move('pyar.log','pyar_old.log')
        reactor_logger.info('====================Reading from Checkpoint====================')
        gamma_list = list(chk.keys()).copy()
        orientations_to_optimize = chk[gamma_list[0]].copy()
        os.chdir('reaction')
        cwd = os.getcwd()
        product_dir = f'{cwd}/products'

    else:
        file_manager.make_directories('reaction')
        os.chdir('reaction')
        cwd = os.getcwd()

        reactor_logger.info('Starting Reactor')
        reactor_logger.info(f'{hm_orientations} orientations will be tried')
        reactor_logger.info(f' Gamma (min): {gamma_min}')
        reactor_logger.info(f' Gamma (max): {gamma_max}')

        reactor_logger.debug(f'Current working directory: {cwd}')

        software = qc_params['software']
        print_header(gamma_max, gamma_min, hm_orientations, software)
        # prepare job directories
        product_dir = f'{cwd}/products'
        reactor_logger.debug(f'Product directory: {product_dir}')
        file_manager.make_directories(product_dir)
        file_manager.make_directories('trial_geometries')
        os.chdir('trial_geometries')
        if site is None:
            all_orientations = tabu.create_trial_geometries('geom', reactant_a,
                                                            reactant_b,
                                                            hm_orientations,
                                                            tabu_on,
                                                            grid_on,
                                                            site)
        else:
            all_orientations = pyar.scan.generate_guess_for_bonding('geom', reactant_a,
                                                                    reactant_b,
                                                                    site[0], site[1],
                                                                    hm_orientations,
                                                                    d_scale=proximity_factor)

        os.chdir(cwd)

        gamma_list = np.linspace(gamma_min, gamma_max, num=10, dtype=float)
        gamma_list = [f"{int(gamma):04d}" for gamma in gamma_list]

        orientations_to_optimize = all_orientations[:]
        # print(k.name for k in orientations_to_optimize)
        chk = {gamma: orientations_to_optimize.copy() for gamma in gamma_list}
        dumpchk(chk, workdir, reactor_logger)

    for en, gamma in enumerate(gamma_list):
        qc_params['gamma'] = gamma
        reactor_logger.info(f'  Current gamma : {gamma}')
        gamma_id = f"{int(gamma):04d}"
        gamma_home = f'{cwd}/gamma_{gamma_id}'
        if not os.path.exists(gamma_home):
            file_manager.make_directories(gamma_home)
        os.chdir(gamma_home)

        optimized_molecules = optimize_all(gamma_id, orientations_to_optimize,chk,
                                           product_dir, qc_params)

        reactor_logger.info(
            f"      {len(optimized_molecules)} geometries from this gamma cycle")
        if len(optimized_molecules) == 0:
            reactor_logger.info(
                "No orientations to be optimized for the next gamma cycle.")
            chk.clear()
            return
        if len(optimized_molecules) == 1:
            orientations_to_optimize = optimized_molecules[:]
        else:
            orientations_to_optimize = clustering.remove_similar(
                optimized_molecules)
        if(en != len(gamma_list)-1):
            chk[gamma_list[en+1]] = orientations_to_optimize
        reactor_logger.info(f"Number of products found from gamma:{gamma} = {len(saved_inchi_strings)}")

        reactor_logger.info(f"{len(orientations_to_optimize)} geometries are considered for the next gamma cycle")

        reactor_logger.debug("the keys of the molecules for next gamma cycle")
        for this_orientation in orientations_to_optimize:
            reactor_logger.debug(f"{this_orientation.name}")
        updtchk(chk,'gamma',gamma,reactor_logger,workdir)

    os.chdir(workdir)
    os.remove('jobs.pkl')
    reactor_logger.info("Removed checkpoints!!")
    return


def optimize_all(gamma_id, orientations, chkdict, product_dir, qc_param):
    gamma = qc_param['gamma']
    cwd = os.getcwd()
    table_of_optimized_molecules = []
    for this_molecule in orientations:
        job_key = this_molecule.name
        reactor_logger.info(f'   Orientation: {job_key}')
        o_key = f"_{job_key[-8:]}"
        orientations_home = f'orientation{o_key}'
        file_manager.make_directories(orientations_home)
        os.chdir(orientations_home)
        job_name = gamma_id + o_key
        this_molecule.name = job_name
        reactor_logger.info(f'Optimizing {this_molecule.name}')
        start_xyz_file_name = f'trial_{this_molecule.name}.xyz'
        this_molecule.mol_to_xyz(start_xyz_file_name)
        start_inchi = pyar.interface.babel.make_inchi_string_from_xyz(start_xyz_file_name)

        start_smile = pyar.interface.babel.make_smile_string_from_xyz(start_xyz_file_name)
        # Update qc_param with the current gamma and index
        # qc_param['gamma'] = gamma
        # qc_param['index'] = len(this_molecule.atoms_list)-1

        status = optimise(this_molecule, qc_param)
        before_relax = copy.copy(this_molecule)
        this_molecule.name = job_name
        reactor_logger.info('... completed')
        if status is True or status == 'converged' or status == 'cycle_exceeded':
            reactor_logger.info("      E({}): {:12.7f}".format(job_name, this_molecule.energy))

            if this_molecule.is_bonded():
                reactor_logger.info("The fragments have close contracts. Going for relaxation")
                this_molecule.mol_to_xyz('trial_relax.xyz')
                this_molecule.name = 'relax'
                status = optimise(this_molecule, qc_param)
                this_molecule.name = job_name
                if status is True or status == 'converged':
                    this_molecule.mol_to_xyz('result_relax.xyz')
                    current_inchi = pyar.interface.babel.make_inchi_string_from_xyz('result_relax.xyz')

                    current_smile = pyar.interface.babel.make_smile_string_from_xyz('result_relax.xyz')

                    reactor_logger.info('geometry relaxed')
                    reactor_logger.info("Checking for product formation with SMILE and InChi strings")

                    reactor_logger.info(f"Start SMILE: {start_smile} Current SMILE: {current_smile}")

                    reactor_logger.info(f"Start InChi: {start_inchi} Current InChi: {current_inchi}")

                    if start_inchi == current_inchi and start_smile == current_smile:
                        table_of_optimized_molecules.append(before_relax)
                        reactor_logger.info(f'{job_name} is added to the table to optimize with higher gamma')

                    else:
                        saved_products[job_name] = this_molecule
                        reactor_logger.info("       The geometry is different from the stating structure.")

                        reactor_logger.info("       Checking if this is a (new) products")
                        if current_inchi in saved_inchi_strings.values() or current_smile in saved_smile_strings.values():
                            reactor_logger.info("Both strings matches with those of already saved products. Discarded")

                        else:
                            reactor_logger.info("        New Product! Saving")
                            saved_inchi_strings[job_name] = current_inchi
                            saved_smile_strings[job_name] = current_smile
                            saved_products[job_name] = this_molecule
                            shutil.copy('result_relax.xyz', f'{product_dir}/{job_name}.xyz')
                        os.chdir(cwd)
                        updtchk(chkdict, 'ori', job_name, reactor_logger, workdir)
                        continue
                elif status == 'cycle_exceeded':
                    table_of_optimized_molecules.append(before_relax)
                    reactor_logger.info(f'{job_name} is added to the table to optimize with higher gamma')

            else:
                table_of_optimized_molecules.append(this_molecule)
                reactor_logger.info('        no close contacts found')
                reactor_logger.info(f'        {job_name} is added to the table to optimize with higher gamma')

        updtchk(chkdict, 'ori', job_name, reactor_logger, workdir)
        os.chdir(cwd)
        sys.stdout.flush()
    return table_of_optimized_molecules


def main():
    pass


if __name__ == "__main__":
    main()
