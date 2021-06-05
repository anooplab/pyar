import copy
import logging
import os
import random
import shutil
import string
from collections import OrderedDict

from pyar import tabu, file_manager
from pyar.data_analysis import clustering
from pyar.optimiser import optimise

aggregator_logger = logging.getLogger('pyar.aggregator')


def random_permutation(iterable, r=None):
    """Random selection from itertools.permutations(iterable, r)"""
    pool = tuple(iterable)
    r = len(pool) if r is None else r
    return tuple(random.sample(pool, r))


def load_molecules(aggregate_sizes, molecules):
    seed_names = string.ascii_lowercase
    ag_id = "ag"
    monomers_to_be_added = []
    for seed_molecule, seed_name, size_of_this_seed in \
            zip(molecules, seed_names, aggregate_sizes):
        seed_molecule.name = seed_name
        for _ in range(size_of_this_seed):
            monomers_to_be_added.append(seed_molecule)
        ag_id += f"_{seed_name}_000"
    return ag_id, monomers_to_be_added


def choose_pathways(monomers_to_be_added, number_of_pathways):
    complete_pathways = set()
    for _ in range(number_of_pathways):
        new_permutation = random_permutation(monomers_to_be_added)
        if new_permutation not in complete_pathways:
            complete_pathways.add(new_permutation)
    return complete_pathways


def aggregate(molecules,
              aggregate_sizes,
              hm_orientations,
              qc_params,
              maximum_number_of_seeds,
              first_pathway,
              number_of_pathways,
              tabu_on,
              grid_on,
              site):
    """
    New aggregate module

    :param grid_on: Toggle the use of grid in generation of trial geometries.
    :param tabu_on: Toggle use of Tabu list in generation of trial geometries.
    :param site: Not used now, but needed for create_trial_molecules().
    :type number_of_pathways: int
    :param number_of_pathways: For cluster or aggregate containing different
        types of molecules or atoms, there are many pathways to explore. This
        parameter determines how many pathways to explore.
    :type first_pathway: int
    :param first_pathway: The starting pathway. This helps in restarting the broken job.
    :param molecules: molecules or atoms for aggregation or cluster formation.
    :type molecules: list(Molecules)
    :param aggregate_sizes: the number of each atoms in the final cluster.
    :type aggregate_sizes: list(int)
    :param hm_orientations: Number of trial orientations.
    :type hm_orientations: int
    :param qc_params: Parameters for Quantum Chemistry Calculations.
    :type qc_params: dict
    :param maximum_number_of_seeds: The maximum number of seeds to be selected
        for the next cycle.
    :type maximum_number_of_seeds: int
    :return: None

    """

    if hm_orientations == 'auto':
        number_of_orientations = 8
    else:
        number_of_orientations = int(hm_orientations)

    parent_folder = 'aggregates'
    file_manager.make_directories(parent_folder)
    os.chdir(parent_folder)
    starting_directory = os.getcwd()

    aggregator_logger.info(f"Starting Aggregation in\n {starting_directory}")

    ag_id, monomers_to_be_added = load_molecules(aggregate_sizes, molecules)

    if len(molecules) == 1:
        pathways_to_calculate = [monomers_to_be_added]
    else:
        pathways_to_calculate = choose_pathways(monomers_to_be_added, number_of_pathways)

    aggregator_logger.info("  The following Afbau paths will be carried out")
    for i, path in enumerate(pathways_to_calculate):
        paths_for_print = f'   {i:03d}: '
        for p in path:
            paths_for_print += p.name
        aggregator_logger.info(paths_for_print)

    seed_storage = OrderedDict()

    initial_storage = copy.deepcopy(seed_storage)
    initial_aggregate_id = ag_id

    outside_counter = first_pathway
    inside_counter = 1

    for i in pathways_to_calculate:
        for this_monomer in i:
            if len(seed_storage) < 1:
                ag_id = update_id(ag_id, this_monomer.name)
                seed_storage[ag_id] = [this_monomer]
                continue
            this_seed = seed_storage[ag_id]
            ag_id = update_id(ag_id, this_monomer.name)
            ag_home = "{}_{:03d}".format(ag_id, outside_counter)
            file_manager.make_directories(ag_home)
            os.chdir(ag_home)

            seed_storage[ag_id] = add_one_to_many_seeds(ag_id,
                                                        this_seed,
                                                        this_monomer,
                                                        number_of_orientations,
                                                        qc_params,
                                                        maximum_number_of_seeds,
                                                        tabu_on, grid_on, site)
            os.chdir(starting_directory)
            if len(seed_storage[ag_id]) == 0:
                aggregator_logger.info(f"No molecules were found from {ag_id}"
                                       f"to continue this pathway.")
                aggregator_logger.info(f"Breaking! \N{worried face}")
                break
            seed_storage.popitem(last=False)
            inside_counter += 1
        outside_counter += 1
        seed_storage = copy.copy(initial_storage)
        ag_id = initial_aggregate_id

    if hm_orientations == 'auto' and number_of_orientations <= 256:
        number_of_orientations += 8
    return


def solvate(seeds, monomer, aggregate_size, hm_orientations,
            qc_params, maximum_number_of_seeds, tabu_on=None, grid_on=None, site=None):
    """
    All monomer to seeds.

    :param seeds:
    :param monomer:
    :param aggregate_size:
    :param hm_orientations:
    :param qc_params:
    :param maximum_number_of_seeds:
    :param tabu_on:
    :param grid_on:
    :param site:
    :return:
    """

    if hm_orientations == 'auto':
        number_of_orientations = 8
    else:
        number_of_orientations = int(hm_orientations)

    starting_directory = os.getcwd()
    aggregator_logger.info(f"Starting Aggregation in\n {starting_directory}")
    for aggregation_counter in range(2, aggregate_size + 2):
        if len(seeds) == 0:
            aggregator_logger.info("No seeds to process")
            return
        aggregate_id = f"{aggregation_counter:03d}"
        aggregate_home = 'aggregate_' + aggregate_id
        file_manager.make_directories(aggregate_home)
        os.chdir(aggregate_home)

        aggregator_logger.info(" Starting aggregation cycle: {}".format(aggregation_counter))

        seeds = add_one_to_many_seeds(aggregate_id, seeds, monomer, number_of_orientations,
                                      qc_params, maximum_number_of_seeds, tabu_on, grid_on, site)

        aggregator_logger.info(f" Aggregation cycle: {aggregation_counter} completed\n")

        if hm_orientations == 'auto' and number_of_orientations <= 256:
            number_of_orientations *= 2
        os.chdir(starting_directory)
    return


def update_id(aid, the_monomer):
    f"""
        aggregate_id = f"a_{aid}:03d}_b_{bid:03d}_c_{cid:03d}"

    """
    from collections import deque
    parts_of_aid = deque(aid.split('_'))
    prefix = parts_of_aid.popleft()
    d = {}
    new_id = prefix
    while parts_of_aid:
        cid = parts_of_aid.popleft()
        n = parts_of_aid.popleft()
        d[cid] = int(n)
        if cid == the_monomer:
            d[cid] += 1
        new_id += f"_{cid}_{d[cid]:03d}"

    return new_id


def add_one_to_many_seeds(aggregate_id, seeds, monomer, hm_orientations, qc_params,
                          maximum_number_of_seeds, tabu_on, grid_on, site):
    """
    Add one monomer to all the seed molecules

    :param tabu_on: Toggle the use of Tabu list
    :type tabu_on: bool
    :param grid_on: Toggle the use of Grid
    :type grid_on: bool
    :param site: Not used
    :type maximum_number_of_seeds: int
    :param maximum_number_of_seeds: The maximum number of seeds to be
        selected for next cycle
    :param qc_params: parameters needed for calculation
    :type qc_params: dict
    :param hm_orientations: Number of orientation to be used.
    :type hm_orientations: int
    :type monomer: Molecule
    :param monomer: monomer molecule
    :type seeds: list[Molecule]
    :param seeds: seed molecules to which monomer will be added.
    :type aggregate_id: str
    :param aggregate_id: An id for the aggregate used for
        job_dir name and xyz file names
    :return: List(Molecule.Molecule)
    :rtype: List(Molecule.Molecule)
    """

    aggregator_logger.info(f'  There are {len(seeds)} seed molecules in {aggregate_id}')

    list_of_optimized_molecules = []
    for seed_count, each_seed in enumerate(seeds):
        cwd = os.getcwd()
        aggregator_logger.info(f'   Seed: {seed_count}')
        seed_id = f"{seed_count:03d}"
        seeds_home = f'seed_{seed_id}'
        file_manager.make_directories(seeds_home)
        os.chdir(seeds_home)
        each_seed.mol_to_xyz('seed.xyz')
        monomer.mol_to_xyz('monomer.xyz')
        mol_id = f'{seed_id}_{aggregate_id}'
        aggregator_logger.debug('Making orientations')
        all_orientations = tabu.create_trial_geometries(mol_id, seeds[seed_count],
                                                        monomer, hm_orientations,
                                                        tabu_on, grid_on, site)
        aggregator_logger.debug('Orientations are made.')
        optimise_all(all_orientations, list_of_optimized_molecules, qc_params)
        os.chdir(cwd)

    if len(list_of_optimized_molecules) < 2:
        selected_seeds = list_of_optimized_molecules
    elif len(list_of_optimized_molecules) <= 8:
        selected_seeds = clustering.remove_similar(list_of_optimized_molecules)
    else:
        aggregator_logger.info("  Clustering")
        selected_seeds = clustering.choose_geometries(list_of_optimized_molecules,
                                                      maximum_number_of_seeds=maximum_number_of_seeds)
    return optimise_further(qc_params, selected_seeds)


def optimise_further(qc_params, selected_seeds):
    file_manager.make_directories('selected')
    os.chdir('selected')
    qc_params["opt_threshold"] = 'normal'
    aggregator_logger.info("  Optimizing the selected molecules with higher threshold")
    less_than_ideal = []
    for each_file in selected_seeds:
        not_refined = copy.deepcopy(each_file)
        status = optimise(each_file, qc_params)
        if status is True:
            xyz_file = 'job_' + each_file.name + '/result_' + each_file.name + '.xyz'
            shutil.copy(xyz_file, '.')
        else:
            selected_seeds.remove(each_file)
            less_than_ideal.append(not_refined)
    if len(selected_seeds) != 0:
        return selected_seeds
    aggregator_logger.info("    The optimization could not be refined, \n"
                           "    so sending the loosely optimised molecules")
    return less_than_ideal


def optimise_all(all_orientations, list_of_optimized_molecules, qc_params):
    not_converged = all_orientations[:]
    status_list = []
    for i in range(10):
        if len(not_converged) != 0:
            aggregator_logger.info(
                f"    Round {i + 1:d} of block optimizations with"
                f" {len(not_converged):d} molecules")
            qc_params["opt_threshold"] = 'loose'
            status_list = [optimise(each_mol, qc_params) for each_mol in
                           not_converged]
            converged = [n for n, s in zip(not_converged, status_list) if s is True]
            not_converged = [n for n, s in zip(not_converged, status_list)
                             if s == 'CycleExceeded' and not tabu.broken(n)]
            not_converged = clustering.remove_similar(not_converged)
            list_of_optimized_molecules.extend(converged)
        else:
            aggregator_logger.info("    All molecules are processed")
            break
    else:
        aggregator_logger.info("    The following molecules are not converged"
                               "after 10 rounds")
        for n, s in zip(not_converged, status_list):
            if s == 'CycleExceeded' and not tabu.broken(n):
                aggregator_logger.info("      ", n.name)


def main():
    pass


if __name__ == "__main__":
    main()
