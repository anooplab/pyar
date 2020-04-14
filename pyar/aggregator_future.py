import os

from pyar import file_manager
from pyar.aggregator import check_stop_signal, aggregator_logger, add_one


def update_id(aid, m):
    if m == 'a':
        n = int(aid[2])
        n += 1
        aid = aid[:2] + str(n) + aid[3:]
        return aid
    if m == 'b':
        n = int(aid[6])
        n += 1
        aid = aid[:6] + str(n) + aid[7:]
        return aid


def exhaustive_binary_aggregate(seed_a_input, seed_b_input,
                                a_n_max, b_n_max,
                                hm_orientations,
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

    parent_folder = 'binary_aggregates'
    if not os.path.exists(parent_folder):
        os.mkdir(parent_folder)
    os.chdir(parent_folder)
    starting_directory = os.getcwd()

    aggregator_logger.info("Starting Aggregation in\n {}".format(starting_directory))

    seed_a = seed_a_input[0]
    seed_b = seed_b_input[0]
    seed_a.name = 'a'
    seed_b.name = 'b'
    a_seeds = [seed_a for _ in range(a_n_max - 1)]
    b_seeds = [seed_b for _ in range(b_n_max - 1)]
    a_n = 1
    b_n = 1

    seed_store = OrderedDict()

    aggregate_id = 'a_' + str(a_n) + '_b_' + str(b_n)
    file_manager.make_directories(aggregate_id)
    os.chdir(aggregate_id)
    seed_store[aggregate_id] = add_one(aggregate_id,
                                       [seed_a], seed_b,
                                       number_of_orientations,
                                       method,
                                       maximum_number_of_seeds)

    os.chdir(starting_directory)

    list_of_monomers = a_seeds + b_seeds
    start_store = copy.copy(seed_store)
    start_id = aggregate_id
    outer_counter = 1
    counter = 1
    for i in set(itertools.permutations(list_of_monomers)):
        for monomer in i:
            seed = seed_store[aggregate_id]
            aggregate_id = update_id(aggregate_id, monomer.name)
            aggregate_home = "{:03d}_{}".format(counter, aggregate_id)
            file_manager.make_directories(aggregate_home)
            os.chdir(aggregate_home)

            seed_store[aggregate_id] = add_one(aggregate_id,
                                               seed, monomer,
                                               number_of_orientations,
                                               method,
                                               maximum_number_of_seeds)
            print(monomer.name, seed_store.keys())
            os.chdir(starting_directory)
            seed_store.popitem(last=False)
            counter += 1

        outer_counter += 1
        seed_store = copy.copy(start_store)
        aggregate_id = start_id

    return


def main():
    pass
