import os
import shutil

import file_manager
import tabu
from data_analysis import clustering
from optimiser import optimise
from tabu import proximity_check


def check_stop_signal():
    if os.path.exists('stop') or os.path.exists('STOP'):
        print("Found stop file, in ", os.getcwd())
        return 1


def aggregate(seeds, monomer, aggregate_size, hm_orientations,
              method, site, number_of_core_atoms, proximity_factor):
    """
    Input: a list of seed molecules, a monomer Molecule objects
    """
    if check_stop_signal():
        print("Function: aggregate")
        return StopIteration

    if hm_orientations == 'auto':
        number_of_orientations = 8
    else:
        number_of_orientations = hm_orientations

    starting_directory = os.getcwd()
    print("Starting Aggregation in\n {}".format(starting_directory))
    for aggregation_counter in range(2, aggregate_size + 2):
        aggregate_id = "{:03d}".format(aggregation_counter)
        aggregate_home = 'aggregate_' + aggregate_id
        file_manager.make_directories(aggregate_home)
        os.chdir(aggregate_home)

        print(" Starting aggregation cycle: {}".format(aggregation_counter))
        seeds = add_one(aggregate_id, seeds, monomer, number_of_orientations,
                        method, site, number_of_core_atoms, proximity_factor)
        print(" Aggregation cycle: {} completed\n".format(aggregation_counter))

        if hm_orientations == 'auto' and number_of_orientations <= 256:
            number_of_orientations *= 2
        os.chdir(starting_directory)
    return


def add_one(aggregate_id, seeds, monomer, hm_orientations, method,
            site, number_of_core_atoms, proximity_factor):
    """

    """
    if check_stop_signal():
        print("Function: add_one")
        return StopIteration
    print('  There are', len(seeds), 'seed molecules')
    cwd = os.getcwd()

    list_of_optimized_molecules = []
    for seed_count, each_seed in enumerate(seeds):
        if check_stop_signal():
            print("Function: add_one")
            return
        print('   Seed: {}'.format(seed_count))
        seed_id = "{:03d}".format(seed_count)
        seeds_home = 'seed_' + seed_id
        file_manager.make_directories(seeds_home)
        os.chdir(seeds_home)
        each_seed.mol_to_xyz('seed.xyz')
        monomer.mol_to_xyz('monomer.xyz')
        mol_id = '{0}_{1}'.format(seed_id, aggregate_id)

        all_orientations = tabu.generate_orientations(mol_id, seeds[seed_count], monomer,
                                                      hm_orientations, site,
                                                      number_of_core_atoms, proximity_factor)
        for molecule in all_orientations:
            o_status = optimise(molecule, method)
            if o_status is True:
                if proximity_check(molecule, site, number_of_core_atoms,
                                   proximity_factor) is False:
                    continue
                print("      E(%10s): %12.7f" % (molecule.name, molecule.energy))
                list_of_optimized_molecules.append(molecule)
            else:
                print('    Optimisation failed:', molecule.name, 'will be discarded')
        os.chdir(cwd)
    if len(list_of_optimized_molecules) < 2:
        return list_of_optimized_molecules
    print("  Clustering")
    selected_seeds = clustering.choose_geometries(list_of_optimized_molecules)
    file_manager.make_directories('selected')
    for each_file in selected_seeds:
        xyz_file = 'seed_' + each_file.name[4:7] + '/job_' + each_file.name + '/result_' + each_file.name + '.xyz'
        shutil.copy(xyz_file, 'selected/')
    return selected_seeds


def main():
    pass


if __name__ == "__main__":
    main()
