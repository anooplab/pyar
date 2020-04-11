import logging
import os
import shutil

from pyar import file_manager, tabu
from pyar.aggregator import check_stop_signal, aggregator_logger, add_one
from pyar.data_analysis import clustering
from pyar.optimiser import optimise


def new_binary_aggregate(seed_a_input, seed_b_input, aggregate_size1, aggregate_size2, hm_orientations,
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

    seed_a = seed_a_input[0]
    seed_b = seed_b_input[0]
    seed_a.name = 'A'
    seed_b.name = 'B'
    a_n = 1
    b_n = 0

    aggregate_counter = 1
    final_one = [seed_a]
    while b_n < aggregate_size2 or a_n < aggregate_size1:
        aggregate_counter += 1
        if b_n < aggregate_size2:
            b_n += 1
            monomer = seed_b
            aggregate_id1 = "{:02d}".format(a_n)
            aggregate_id2 = "{:02d}".format(b_n)

            aggregate_home = 'a_' + aggregate_id1 + '_b_' + aggregate_id2
            aggregate_id = aggregate_id1 + '_' + aggregate_id2

            file_manager.make_directories(aggregate_home)
            os.chdir(aggregate_home)

            aggregator_logger.info(" Starting aggregation cycle: {}".format(aggregate_counter))

            final_one = add_one(aggregate_id, final_one,
                                monomer, number_of_orientations, method,
                                maximum_number_of_seeds)

            aggregator_logger.info(" Aggregation cycle: {} completed\n".format(aggregate_counter))

            os.chdir(starting_directory)

        if a_n < aggregate_size1:
            a_n += 1
            monomer = seed_a
            aggregate_id1 = "{:02d}".format(a_n)
            aggregate_id2 = "{:02d}".format(b_n)

            aggregate_home = 'a_' + aggregate_id1 + '_b_' + aggregate_id2
            aggregate_id = aggregate_id1 + '_' + aggregate_id2

            file_manager.make_directories(aggregate_home)
            os.chdir(aggregate_home)

            aggregator_logger.info(" Starting aggregation cycle: {}".format(aggregate_counter))

            final_one = add_one(aggregate_id, final_one,
                                monomer, number_of_orientations, method,
                                maximum_number_of_seeds)

            aggregator_logger.info(" Aggregation cycle: {} completed\n".format(aggregate_counter))

            os.chdir(starting_directory)

            if hm_orientations == 'auto' and number_of_orientations <= 256:
                number_of_orientations += 8
            os.chdir(starting_directory)


aggregator_logger = logging.getLogger('pyar.aggregator')


def check_stop_signal():
    if os.path.exists('stop') or os.path.exists('STOP'):
        aggregator_logger.info("Found stop file, in {}".format(os.getcwd()))
        return 1


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