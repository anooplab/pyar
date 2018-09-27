import os
import shutil

import file_manager
import tabu
from data_analysis import clustering
from optimiser import optimise
import logging
aggregator_logger = logging.getLogger('pyar.aggregator')


def check_stop_signal():
    if os.path.exists('stop') or os.path.exists('STOP'):
        aggregator_logger.info("Found stop file, in {}".format(os.getcwd()))
        return 1


def aggregate(seeds, monomer, aggregate_size, hm_orientations, method):
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
        seeds = add_one(aggregate_id, seeds, monomer, number_of_orientations, method)
        aggregator_logger.info(" Aggregation cycle: {} completed\n".format(aggregation_counter))

        if hm_orientations == 'auto' and number_of_orientations <= 256:
            number_of_orientations *= 2
        os.chdir(starting_directory)
    return


def add_one(aggregate_id, seeds, monomer, hm_orientations, method):
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
        for molecule in all_orientations:
            o_status = optimise(molecule, method)
            if o_status is True and molecule.energy is not None:
                aggregator_logger.info("      E(%10s): %12.7f" % (molecule.name, molecule.energy))
                list_of_optimized_molecules.append(molecule)
            else:
                aggregator_logger.info('    Optimisation failed: {}, will '
                                       'be discarded'.format(molecule.name))
        os.chdir(cwd)

    if len(list_of_optimized_molecules) < 2:
        return list_of_optimized_molecules
    aggregator_logger.info("  Clustering")
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
