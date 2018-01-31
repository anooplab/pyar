import copy
import operator
import os
import shutil
import sys

import numpy as np

import file_manager
import interface.babel
import tabu
from afir import fragment
from data_analysis import clustering
from optimiser import optimise

table_of_product_molecules = {}
table_of_product_inchi_strings = {}
table_of_product_smile_strings = {}


def print_header(gamma_max, gamma_min, hm_orientations, software):
    print('Starting PyAR 2.0\n')
    print(hm_orientations, 'orientations will be tried.')
    print(' Gamma (min): ', gamma_min)
    print(' Gamma (max): ', gamma_max)
    print(' Software   : ', software)


def react(reactant_a, reactant_b, gamma_min, gamma_max,
          hm_orientations, method,
          site, number_of_core_atom, proximity_factor):
    cwd = os.getcwd()
    software = method['software']
    print_header(gamma_max, gamma_min, hm_orientations, software)
    # prepare job directories
    product_dir = cwd + '/products'
    file_manager.make_directories(product_dir)
    file_manager.make_directories('trial_geometries')
    os.chdir('trial_geometries')

    all_orientations = tabu.new_func('geom',
                                     reactant_a,
                                     reactant_b,
                                     hm_orientations,
                                     site,
                                     number_of_core_atom,
                                     proximity_factor)

    os.chdir(cwd)

    gamma_list = np.linspace(gamma_min, gamma_max, num=10, dtype=float)
    orientations_to_optimize = all_orientations.copy()

    for gamma in gamma_list:
        print('  Current gamma :', gamma)
        gamma_id = "%04d" % (int(gamma))
        gamma_home = cwd + '/gamma_' + gamma_id
        file_manager.make_directories(gamma_home)
        os.chdir(gamma_home)

        optimized_molecules = optimize_all(gamma_id, gamma,
                                           orientations_to_optimize,
                                           product_dir, method)

        print("      ", len(optimized_molecules), "geometries from this gamma cycle")
        orientations_to_optimize = clustering.remove_similar(optimized_molecules)
        if len(orientations_to_optimize) == 0:
            print("No orientations to optimized for the next gamma cycle.")
            break
        print("      ", len(orientations_to_optimize), "geometries in the next gamma cycle")
        print("number of products found from gamma:", gamma
              , " = ", len(table_of_product_inchi_strings))
        for key, value in orientations_to_optimize.items():
            print("the key for next gamma cycle:", key)
    print("\n\n\n\n")
    os.chdir(cwd)
    return


def optimize_all(gamma_id, gamma, orientations_to_optimize,
                 product_dir, method):
    cwd = os.getcwd()
    table_of_optimized_molecules = {}
    for job_key, this_molecule in sorted(orientations_to_optimize.items(), key=operator.itemgetter(0)):
        print('   Orientation:', job_key)
        o_key = u"_{}".format(job_key[-8:])
        orientations_home = 'orientation' + o_key
        file_manager.make_directories(orientations_home)
        os.chdir(orientations_home)
        na, nb = [len(i) for i in this_molecule.fragments]
        fragment.make_fragment_file(na, nb)
        job_name = gamma_id + o_key
        this_molecule.name = job_name
        print('    Optimizing', this_molecule.name, ':')
        start_xyz_file_name = 'trial_' + this_molecule.name + '.xyz'
        this_molecule.mol_to_xyz(start_xyz_file_name)
        start_inchi = interface.babel.make_inchi_string_from_xyz(start_xyz_file_name)
        start_smile = interface.babel.make_smile_string_from_xyz(start_xyz_file_name)
        status = optimise(this_molecule, method, gamma=gamma)
        before_relax = copy.copy(this_molecule)
        this_molecule.name = job_name
        print('     job completed')
        if status is True or status == 'converged' or status == 'cycle_exceeded':
            print("      E(%s): %12.7f" % (job_name, this_molecule.energy))
            if this_molecule.is_bonded() is True:
                print("      The fragments have close contracts. Going for relaxation")
                this_molecule.mol_to_xyz('trial_relax.xyz')
                this_molecule.name = 'relax'
                status = optimise(this_molecule, method)
                this_molecule.name = job_name
                if status is True or status == 'converged' or status == 'cycle_exceeded':
                    current_inchi = interface.babel.make_inchi_string_from_xyz('result_relax.xyz')
                    current_smile = interface.babel.make_smile_string_from_xyz('result_relax.xyz')
                    print('      geometry relaxed')
                    print("Checking for product formation with SMILE and InChi strings")
                    print("Start SMILE:", start_smile, "Current SMILE:",
                          current_smile)
                    print("Start InChi:", start_inchi, "Current InChi:",
                          current_inchi)

                    if start_inchi != current_inchi or start_smile != current_smile:
                        table_of_product_molecules[job_name] = this_molecule
                        print("       The geometry is different from the stating structure.")
                        print("       Checking if this is a (new) products")
                        if current_inchi not in table_of_product_inchi_strings.values():
                            if current_smile not in table_of_product_smile_strings.values():
                                print("        New Product! Saving")
                                table_of_product_inchi_strings[job_name] = current_inchi
                                table_of_product_smile_strings[job_name] = current_smile
                                table_of_product_molecules[job_name] = this_molecule
                                shutil.copy('result_relax.xyz',
                                            product_dir + '/' + job_name + '.xyz')
                                os.chdir(cwd)
                                continue
                            else:
                                print("SMILE matches")
                                os.chdir(cwd)
                                continue
                        else:
                            print("InChi matches")
                            os.chdir(cwd)
                            continue
                    else:
                        table_of_optimized_molecules[job_name] = before_relax
                        print(job_name, 'is added to the table to optimize with higher gamma')
            else:
                table_of_optimized_molecules[job_name] = this_molecule
                print('        no close contacts found')
                print('       ', job_name, 'is added to the table to optimize with higher gamma')
        os.chdir(cwd)
        sys.stdout.flush()
    return table_of_optimized_molecules


def main():
    pass


if __name__ == "__main__":
    main()
