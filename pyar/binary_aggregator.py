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
