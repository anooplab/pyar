import operator
import os
import shutil

import file_manager
import tabu
from data_analysis import clustering
from optimiser import optimise


def check_stop_signal():
    if os.path.exists('stop') or os.path.exists('STOP'):
        print("Found stop file, in ", os.getcwd())
        return 1


def aggregate(seeds, monomer, aggregate_size=2, hm_orientations=8,
              method=None):
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
        seeds = add_one(aggregate_id, seeds, monomer, number_of_orientations, method)
        print(" Aggregation cycle: {} completed\n".format(aggregation_counter))

        if hm_orientations == 'auto' and number_of_orientations <= 256:
            number_of_orientations *= 2
        os.chdir(starting_directory)
    return


def add_one(aggregate_id, seeds, monomer, hm_orientations, method):
    """
    :type aggregate_id str
    :type seeds list of Molecules
    :type monomer Molecule.Molecule
    :type hm_orientations int how many orientations
    :type method dict containing charge, multiplicity, software
    """
    if check_stop_signal():
        print("Function: add_one")
        return StopIteration
    print('  There are', len(seeds), 'seed molecules')
    cwd = os.getcwd()

    dict_of_optimized_molecules = {}
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

        all_orientations = tabu.generate_orientations(mol_id, seeds[seed_count], monomer, hm_orientations)
        for name, molecule in sorted(all_orientations.items(), key=operator.itemgetter(0)):
            o_status = optimise(molecule, method)
            if o_status is True:
                print("      E(%10s): %12.7f" % (name, molecule.energy))
                dict_of_optimized_molecules[name] = molecule
            else:
                print('    Optimisation failed:', name, 'will be discarded')
        os.chdir(cwd)
    if len(dict_of_optimized_molecules) < 2:
        return list(dict_of_optimized_molecules.values())
    print("  Clustering")
    selected_seeds = clustering.choose_geometries(dict_of_optimized_molecules)
    file_manager.make_directories('selected')
    for each_file in selected_seeds.values():
        xyz_file = 'seed_' + each_file.name[4:7] + '/result_' + each_file.name + '.xyz'
        shutil.copy(xyz_file, 'selected/')
    return list(selected_seeds.values())

def main():
    pass

if __name__ == "__main__":
    main()
