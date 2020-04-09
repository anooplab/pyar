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


def aggregate_in_binary(seeds1, seeds2, aggregate_size1, aggregate_size2, hm_orientations,
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
    parent_folder = 'binaryAggregates'
    if not os.path.exists(parent_folder):
        os.mkdir(parent_folder)
    os.chdir(parent_folder)
    starting_directory = os.getcwd()

    aggregator_logger.info("Starting Aggregation in\n {}".format(starting_directory))
    seed_init1 = seeds1
    seed_init2 = seeds2
    seed_in1 = seeds1
    seed_in2 = seeds2
    for aggregation_counter1 in range(1, aggregate_size1 + 1):
        for aggregation_counter2 in range(1, aggregate_size2 + 1):
            aggregator_logger.debug("Seed one is {}".format(seeds1))
            aggregator_logger.debug("Seed two is {}".format(seeds2))

            if aggregation_counter1 > 1 and aggregation_counter2 == 1:
                print("Not doing", aggregation_counter1, aggregation_counter2)
            else:
                aggregate_id1 = "{:02d}".format(aggregation_counter1)
                aggregate_id2 = "{:02d}".format(aggregation_counter2)
                aggregate_home = 'a_' + aggregate_id1 + '_b_' + aggregate_id2
                file_manager.make_directories(aggregate_home)
                os.chdir(aggregate_home)

                aggregator_logger.info(" Starting aggregation cycle: {}".format(aggregation_counter1))
                print(aggregate_home)
                seeds1 = add_one_in_binary(aggregate_id1, aggregate_id2, seeds1,
                                           seed_init2, number_of_orientations, method,
                                           maximum_number_of_seeds)

                aggregator_logger.info(" Aggregation cycle: {} completed\n".format(aggregation_counter1))

            if aggregation_counter1 < aggregate_size1 and aggregation_counter2 == 1:
                os.chdir(starting_directory)
                aggregate_id1 = "{:02d}".format(aggregation_counter1 + 1)
                aggregate_id2 = "{:02d}".format(aggregation_counter2)
                aggregate_home = 'a_' + aggregate_id1 + '_b_' + aggregate_id2
                file_manager.make_directories(aggregate_home)
                os.chdir(aggregate_home)
                seed_in1 = seeds1
                print(aggregate_home)
                seed_in1 = add_one_in_binary(aggregate_id1, aggregate_id2, seed_in1,
                                             seed_init1, number_of_orientations, method,
                                             maximum_number_of_seeds)

            if aggregation_counter2 == aggregate_size2:
                seeds1 = seed_in1
            if hm_orientations == 'auto' and number_of_orientations <= 256:
                number_of_orientations += 8
            os.chdir(starting_directory)


def add_one_in_binary(aggregate_id1, aggregate_id2, seeds1, seeds2, hm_orientations, method, maximum_number_of_seeds):
    """

    """
    seeds2 = seeds2[0]
    if check_stop_signal():
        aggregator_logger.info("Function: add_one")
        return StopIteration

    aggregator_logger.info('  There are {} seed molecules'.format(len(seeds1)))
    cwd = os.getcwd()

    list_of_optimized_molecules = []
    for seed_count, each_seed in enumerate(seeds1):
        if check_stop_signal():
            aggregator_logger.info("Function: add_one")
            return
        aggregator_logger.info('   Seed: {}'.format(seed_count))
        seed_id = "{:03d}".format(seed_count)
        seeds_home = 'seed_' + seed_id
        file_manager.make_directories(seeds_home)
        os.chdir(seeds_home)
        each_seed.mol_to_xyz('seed.xyz')
        seeds2.mol_to_xyz('monomer.xyz')
        mol_id = '{0}_{1}_{2}'.format(seed_id, aggregate_id1, aggregate_id2)

        all_orientations = tabu.generate_orientations(mol_id, seeds1[seed_count], seeds2, hm_orientations)
        not_converged = all_orientations[:]
        converged = []
        for i in range(10):
            aggregator_logger.info("Round %d of block optimizations with %d molecules" % (i + 1, len(not_converged)))
            if len(not_converged) == 0:
                break
            status_list = [optimise(each_mol, method, max_cycles=700, convergence='loose') for each_mol in
                           not_converged]
            converged = [n for n, s in zip(not_converged, status_list) if s is True]
            list_of_optimized_molecules.extend(converged)
            not_converged = [n for n, s in zip(not_converged, status_list) if s == 'CycleExceeded']
            not_converged = clustering.remove_similar(not_converged)

        os.chdir(cwd)
    aggregator_logger.debug(list_of_optimized_molecules)
    if len(list_of_optimized_molecules) < 2:
        return list_of_optimized_molecules
    aggregator_logger.info("  Clustering")
    selected_seeds = clustering.choose_geometries(list_of_optimized_molecules,
                                                  maximum_number_of_seeds=maximum_number_of_seeds)
    file_manager.make_directories('selected')
    for each_file in selected_seeds:
        status = optimise(each_file, method, max_cycles=500, convergence='normal')
        if status is True:
            xyz_file = 'seed_' + each_file.name[4:7] + '/job_' + each_file.name + '/result_' + each_file.name + '.xyz'
            shutil.copy(xyz_file, 'selected/')
        else:
            selected_seeds.remove(each_file)
    return selected_seeds


def main():
    pass


if __name__ == "__main__":
    main()
