import os
import shutil

from pyar import tabu, file_manager
from pyar.data_analysis import clustering
from pyar.optimiser import optimise
import logging

aggregator_logger = logging.getLogger('pyar.aggregator')


def check_stop_signal():
    if os.path.exists('stop') or os.path.exists('STOP'):
        aggregator_logger.info("Found stop file, in {}".format(os.getcwd()))
        return 1


def aggregate(seeds, monomer, aggregate_size, hm_orientations, method, maximum_number_of_seeds):
    """
    Input: a list of seed molecules, a monomer Molecule objects
    """
    if check_stop_signal():
        aggregator_logger.info("Function: aggregate")
        return StopIteration

    if hm_orientations == 'auto':
        number_of_orientations = 8
    else:
        number_of_orientations = int(hm_orientations)

    starting_directory = os.getcwd()
    aggregator_logger.info("Starting Aggregation in\n {}".format(starting_directory))
    for aggregation_counter in range(2, aggregate_size + 2):
        aggregate_id = "{:03d}".format(aggregation_counter)
        aggregate_home = 'aggregate_' + aggregate_id
        file_manager.make_directories(aggregate_home)
        os.chdir(aggregate_home)

        aggregator_logger.info(" Starting aggregation cycle: {}".format(aggregation_counter))
        seeds = add_one(aggregate_id, seeds, monomer, number_of_orientations, method, maximum_number_of_seeds)
        aggregator_logger.info(" Aggregation cycle: {} completed\n".format(aggregation_counter))

        if hm_orientations == 'auto' and number_of_orientations <= 256:
            number_of_orientations *= 2
        os.chdir(starting_directory)
    return


def add_one(aggregate_id, seeds, monomer, hm_orientations, method, maximum_number_of_seeds):
    """

    """
    if check_stop_signal():
        aggregator_logger.info("Function: add_one")
        return StopIteration

    aggregator_logger.info('  There are {} seed molecules'.format(len(seeds)))
    cwd = os.getcwd()

    list_of_optimized_molecules = []
    for seed_count, each_seed in enumerate(seeds):
        if check_stop_signal():
            aggregator_logger.info("Function: add_one")
            return
        aggregator_logger.info('   Seed: {}'.format(seed_count))
        seed_id = "{:03d}".format(seed_count)
        seeds_home = 'seed_' + seed_id
        file_manager.make_directories(seeds_home)
        os.chdir(seeds_home)
        each_seed.mol_to_xyz('seed.xyz')
        monomer.mol_to_xyz('monomer.xyz')
        mol_id = '{0}_{1}'.format(seed_id, aggregate_id)
        aggregator_logger.debug('Making orientations')
        all_orientations = tabu.generate_orientations(mol_id, seeds[seed_count], monomer, hm_orientations)
        aggregator_logger.debug('Orientations are made.')

        not_converged = all_orientations[:]
        for i in range(10):
            aggregator_logger.info("Round %d of block optimizations with %d molecules" % (i + 1, len(not_converged)))
            if len(not_converged) == 0:
                aggregator_logger.info("No more files")
                break
            status_list = [optimise(each_mol, method, max_cycles=100, convergence='loose') for each_mol in
                           not_converged]
            converged = [n for n, s in zip(not_converged, status_list) if s is True]
            list_of_optimized_molecules.extend(converged)
            not_converged = [n for n, s in zip(not_converged, status_list) if s == 'CycleExceeded']
            not_converged = clustering.remove_similar(not_converged)

        os.chdir(cwd)

    if len(list_of_optimized_molecules) < 2:
        return list_of_optimized_molecules
    aggregator_logger.info("  Clustering")
    selected_seeds = clustering.choose_geometries(list_of_optimized_molecules,
                                                  maximum_number_of_seeds=maximum_number_of_seeds)
    file_manager.make_directories('selected')
    for each_file in selected_seeds:
        status = optimise(each_file, method, max_cycles=100, convergence='normal')
        if status is True:
            xyz_file = 'seed_' + each_file.name[4:7] + '/job_' + each_file.name + '/result_' + each_file.name + '.xyz'
            shutil.copy(xyz_file, 'selected/')
        else:
            selected_seeds.remove(each_file)
    return selected_seeds


def binary_aggregate(seed_a_input, seed_b_input, a_n_max, b_n_max, hm_orientations,
                     method, maximum_number_of_seeds):
    """
    Input: a list of seed molecules, a monomer Molecule objects
    """
    if check_stop_signal():
        aggregator_logger.info("Function: aggregate")
        return StopIteration

    if hm_orientations == 'auto':
        number_of_orientations = 8
    else:
        number_of_orientations = int(hm_orientations)

    starting_directory = os.getcwd()

    parent_folder = 'binary_aggregates'

    if not os.path.exists(parent_folder):
        os.mkdir(parent_folder)

    aggregator_logger.info("Starting Aggregation in\n {}".format(starting_directory))

    current_seed_a = seed_a_input
    current_seed_b = seed_b_input
    tmp_holder = seed_a_input

    for a_counter in range(1, a_n_max + 1):
        for b_counter in range(1, b_n_max + 1):
            print(a_counter, b_counter, os.getcwd())
            a_n = "{:02d}".format(a_counter)
            b_n = "{:02d}".format(b_counter)

            aggregate_home = 'a_' + a_n + '_b_' + b_n
            aggregate_id = a_n + '_' + b_n

            if not (a_counter > 1 and b_counter == 1):

                os.chdir(parent_folder)
                file_manager.make_directories(aggregate_home)
                os.chdir(aggregate_home)

                aggregator_logger.info(" Starting aggregation cycle: {}".format(a_counter))

                seed_a_input = add_one(aggregate_id, seed_a_input,
                                       current_seed_b[0], number_of_orientations, method,
                                       maximum_number_of_seeds)

                aggregator_logger.info(" Aggregation cycle: {} completed\n".format(a_counter))

                os.chdir(starting_directory)
            if a_counter < a_n_max and b_counter == 1:
                os.chdir(parent_folder)
                file_manager.make_directories(aggregate_home)
                os.chdir(aggregate_home)
                tmp_holder = seed_a_input
                tmp_holder = add_one(aggregate_id, tmp_holder,
                                     current_seed_a[0], number_of_orientations, method,
                                     maximum_number_of_seeds)
                os.chdir(starting_directory)

            if b_counter == b_n_max:
                seed_a_input = tmp_holder
            if hm_orientations == 'auto' and number_of_orientations <= 256:
                number_of_orientations += 8


def main():
    pass


if __name__ == "__main__":
    main()
