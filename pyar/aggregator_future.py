import os

from pyar import file_manager
from pyar.aggregator import check_stop_signal, aggregator_logger, add_one
from pyar.binary_aggregator import check_stop_signal, aggregator_logger


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

    seed = seed_a
    monomer = seed_b
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
