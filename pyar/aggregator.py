import copy
import logging
import os
import random
import shutil
import string
from collections import OrderedDict
import numpy as np
from pyar import tabu, file_manager
from pyar.Molecule import Molecule
from pyar.data_analysis import clustering
from pyar.old_optimiser import optimise
import re
from pyar.Molecule import atomic_data
aggregator_logger = logging.getLogger('pyar.aggregator')


def random_permutation(iterable, r=None):
    """Random selection from itertools.permutations(iterable, r)"""
    pool = tuple(iterable)
    r = len(pool) if r is None else r
    return tuple(random.sample(pool, r))


def read_old_path():
    path = []
    needed_elements = 3 if 'DEBUG' in open('pyar.log') else 2
    for line in open('pyar.log'):
        split_line = line.split(':')
        if len(split_line) == needed_elements and split_line[
            -2].strip().isnumeric():
            path.append(split_line[-1].rstrip())
    return path or None


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
    :param first_pathway: The starting pathway. This helps in restarting the
    broken job.
    :param molecules: molecules or atoms for aggregation or cluster formation.
    :type molecules: list(Molecules)
    :param aggregate_sizes: the number of each atom in the final cluster.
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

    if check_stop_signal():
        aggregator_logger.info("Function: aggregate")
        return StopIteration

    if hm_orientations == 'auto':
        number_of_orientations = 8
    else:
        number_of_orientations = int(hm_orientations)

    parent_folder = 'aggregates'
    old_path = read_old_path()
    restart = bool(old_path)

    if not restart:
        file_manager.make_directories(parent_folder)

    os.chdir(parent_folder)
    starting_directory = os.getcwd()

    if restart:
        aggregator_logger.info(
            f"Restarting Aggregation in\n {starting_directory}")
    else:
        aggregator_logger.info(
            f"Starting Aggregation in\n {starting_directory}")

    seed_names = string.ascii_lowercase
    ag_id = "ag"

    monomers_to_be_added = []
    for seed_molecule, seed_name, size_of_this_seed in zip(molecules,
                                                           seed_names,
                                                           aggregate_sizes):
        seed_molecule.name = seed_name
        monomers_to_be_added.extend(
            seed_molecule for _ in range(size_of_this_seed))
        ag_id += f"_{seed_name}_000"

    if len(molecules) == 1:
        pathways_to_calculate = [monomers_to_be_added]
    elif restart:
        pathways_to_calculate = old_path_to_new_path(monomers_to_be_added,
                                                     old_path)
    else:
        pathways_to_calculate = select_pathways(monomers_to_be_added,
                                                number_of_pathways)

        aggregator_logger.info(
            "  The following Afbau paths will be carried out")
        for i, path in enumerate(pathways_to_calculate):
            paths_for_print = f'      {i:03d}: '
            for p in path:
                paths_for_print += p.name
            aggregator_logger.info(paths_for_print)

    seed_storage = OrderedDict()

    initial_storage = copy.deepcopy(seed_storage)
    initial_aggregate_id = ag_id

    outside_counter = first_pathway
    inside_counter = 1

    for i in pathways_to_calculate:
        aggregator_logger.info(f"  Path: {i}")
        for this_monomer in i:
            if len(seed_storage) < 1:
                ag_id = update_id(ag_id, this_monomer.name)
                seed_storage[ag_id] = [this_monomer]
                continue
            this_seed = seed_storage[ag_id]
            ag_id = update_id(ag_id, this_monomer.name)
            ag_home = "{}_{:03d}".format(ag_id, outside_counter)
            if not os.path.exists(ag_home):
                file_manager.make_directories(ag_home)
            os.chdir(ag_home)

            seed_storage[ag_id] = add_one(ag_id,
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
                aggregator_logger.info('Breaking! 😟')
                break
            seed_storage.popitem(last=False)
            inside_counter += 1
        outside_counter += 1
        seed_storage = copy.copy(initial_storage)
        ag_id = initial_aggregate_id

    if hm_orientations == 'auto' and number_of_orientations <= 256:
        number_of_orientations += 8
    return


def old_path_to_new_path(monomers_to_be_added, old_path):
    complete_pathways = []
    for each in old_path:
        tmp = []
        for e in each:
            tmp.extend(a for a in monomers_to_be_added if a.name == e)
        complete_pathways.append(tuple(tmp))
    return complete_pathways


def select_pathways(monomers_to_be_added, number_of_pathways):
    complete_pathways = set()
    for _ in range(number_of_pathways):
        new_permutation = random_permutation(monomers_to_be_added)
        if new_permutation not in complete_pathways:
            complete_pathways.add(new_permutation)
    return complete_pathways


def solvate(seeds, monomer, aggregate_size, hm_orientations,
            qc_params, maximum_number_of_seeds, tabu_on=None, grid_on=None,
            site=None):
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
    if check_stop_signal():
        aggregator_logger.info("Function: solvate")
        return StopIteration

    if hm_orientations == 'auto':
        number_of_orientations = 8
    else:
        number_of_orientations = int(hm_orientations)

    starting_directory = os.getcwd()
    aggregator_logger.info(
        f"Starting Aggregation in\n {starting_directory}")
    for aggregation_counter in range(2, aggregate_size + 2):
        if len(seeds) == 0:
            aggregator_logger.info("No seeds to process")
            return
        aggregate_id = "{:03d}".format(aggregation_counter)
        aggregate_home = f'aggregate_{aggregate_id}'
        file_manager.make_directories(aggregate_home)
        os.chdir(aggregate_home)

        aggregator_logger.info(
            f" Starting aggregation cycle: {aggregation_counter}")

        seeds = add_one(aggregate_id, seeds, monomer, number_of_orientations,
                        qc_params, maximum_number_of_seeds, tabu_on, grid_on,
                        site)

        aggregator_logger.info(
            f" Aggregation cycle: {aggregation_counter} completed\n")

        if hm_orientations == 'auto' and number_of_orientations <= 256:
            number_of_orientations *= 2
        os.chdir(starting_directory)
    return


def read_orientations(molecule_id, noo):
    orientations = []
    for i in range(noo):
        xyz_file = f"trial_{i:03d}_{molecule_id}.xyz"
        new_orientation = Molecule.from_xyz(xyz_file)
        new_orientation.name = f'{i:03d}_{molecule_id}'
        orientations.append(new_orientation)
    return orientations




# def add_one(aggregate_id, seeds, monomer, hm_orientations, qc_params, maximum_number_of_seeds, tabu_on, grid_on, site):
    
#     if check_stop_signal():
#         aggregator_logger.info("Function: add_one")
#         return StopIteration
#     aggregator_logger.info(f'  There are {len(seeds)} seed molecules in {aggregate_id}')

#     cwd = os.getcwd()
#     list_of_optimized_molecules = []
#     for seed_count, each_seed in enumerate(seeds):
#         if check_stop_signal():
#             aggregator_logger.info("Function: add_one")
#             return
#         aggregator_logger.info(f'   Seed: {seed_count}')
#         seed_id = "{:03d}".format(seed_count)
#         seeds_home = f'seed_{seed_id}'
#         if not os.path.exists(seeds_home):
#             file_manager.make_directories(seeds_home)
#         os.chdir(seeds_home)
#         each_seed.mol_to_xyz('seed.xyz')
#         monomer.mol_to_xyz('monomer.xyz')
#         if len(each_seed) == 1:
#             hm_orientations = 1
#         mol_id = f'{seed_id}_{aggregate_id}'
#         aggregator_logger.debug('Making orientations')
#         if not all(os.path.exists(f"trial_{i:03d}_{mol_id}.xyz") for i in range(hm_orientations)):
#             all_orientations = tabu.create_trial_geometries(mol_id, seeds[seed_count], monomer, hm_orientations, tabu_on, grid_on, site)

#             aggregator_logger.debug('Orientations are made.')
#         else:
#             all_orientations = read_orientations(mol_id, hm_orientations)
#         not_converged = all_orientations[:]
#         status_list = [False for _ in not_converged]
#         for i in range(10):
#             if len(not_converged) > 0:
#                 aggregator_logger.info(f"    Round {i + 1:d} of block optimizations with {len(not_converged):d} molecules")

#                 qc_params["opt_threshold"] = 'loose'
#                 status_list = [optimise(each_mol, qc_params) for each_mol in not_converged]
#                 converged = [n for n, s in zip(not_converged, status_list) if s is True]
#                 list_of_optimized_molecules.extend(converged)
#                 not_converged = [n for n, s in zip(not_converged, status_list) if s == 'CycleExceeded' and not tabu.broken(n)]

#                 not_converged = clustering.remove_similar(not_converged)
#             else:
#                 aggregator_logger.info("    All molecules are processed")
#                 break
#         else:
#             aggregator_logger.info("    The following molecules are not convergedafter 10 rounds")

#             for n, s in zip(not_converged, status_list):
#                 if s == 'CycleExceeded' and not tabu.broken(n):
#                     aggregator_logger.info("      ", n.name)
#         os.chdir(cwd)
#     if os.path.exists('selected'):
#         os.chdir('selected')
#         optimized_molecules = [i.name for i in list_of_optimized_molecules]
#         job_done = []
#         selected_seeds = []
#         for i in optimized_molecules:
#             for j in os.listdir():
#                 if f'job_{i}' == j:
#                     job_done.append(j)
#                     if f'result_{i}.xyz' == j:
#                         selected_seeds.append(i)
#                         list_of_optimized_molecules.pop(j)
#         os.chdir(cwd)
#     else:
#         file_manager.make_directories('selected')
#     if len(list_of_optimized_molecules) < 2:
#         selected_seeds = list_of_optimized_molecules
#     else:
#         aggregator_logger.info("  Clustering")
#         selected_seeds = clustering.choose_geometries(list_of_optimized_molecules, maximum_number_of_seeds=maximum_number_of_seeds)

#     os.chdir('selected')
#     qc_params["opt_threshold"] = 'normal'
#     aggregator_logger.info("Optimizing the selected molecules with higher threshold")

#     less_than_ideal = []
#     for each_file in selected_seeds:
#         not_refined = copy.deepcopy(each_file)
#         status = optimise(each_file, qc_params)
#         if status is True:
#             xyz_file = f'job_{each_file.name}/result_{each_file.name}.xyz'
#             shutil.copy(xyz_file, '.')
#         else:
#             selected_seeds.remove(each_file)
#             less_than_ideal.append(not_refined)
#     if len(selected_seeds) != 0:
#         return selected_seeds
#     aggregator_logger.info("    The optimization could not be refined, \n    so sending the loosely optimised molecules")

#     return less_than_ideal

def add_one(aggregate_id, seeds, monomer, hm_orientations, qc_params, maximum_number_of_seeds, tabu_on, grid_on, site):
    if check_stop_signal():
        aggregator_logger.info("Function: add_one")
        return StopIteration
    aggregator_logger.info(f'  There are {len(seeds)} seed molecules in {aggregate_id}')

    cwd = os.getcwd()
    
    if qc_params.get('software'):
        
        list_of_optimized_molecules = []
        for seed_count, each_seed in enumerate(seeds):
            if check_stop_signal():
                aggregator_logger.info("Function: add_one")
                return
            aggregator_logger.info(f'   Seed: {seed_count}')
            seed_id = "{:03d}".format(seed_count)
            seeds_home = f'seed_{seed_id}'
            if not os.path.exists(seeds_home):
                file_manager.make_directories(seeds_home)
            os.chdir(seeds_home)
            each_seed.mol_to_xyz('seed.xyz')
            monomer.mol_to_xyz('monomer.xyz')
            if len(each_seed) == 1:
                hm_orientations = 1
            mol_id = f'{seed_id}_{aggregate_id}'
            aggregator_logger.debug('Making orientations')
            if not all(os.path.exists(f"trial_{i:03d}_{mol_id}.xyz") for i in range(hm_orientations)):
                all_orientations = tabu.create_trial_geometries(mol_id, seeds[seed_count], monomer, hm_orientations, tabu_on, grid_on, site)
                aggregator_logger.debug('Orientations are made.')
            else:
                all_orientations = read_orientations(mol_id, hm_orientations)
            not_converged = all_orientations[:]
            status_list = [False for _ in not_converged]
            for i in range(10):
                if len(not_converged) > 0:
                    aggregator_logger.info(f"    Round {i + 1:d} of block optimizations with {len(not_converged):d} molecules")
                    qc_params["opt_threshold"] = 'loose'
                    status_list = [optimise(each_mol, qc_params) for each_mol in not_converged]
                    converged = [n for n, s in zip(not_converged, status_list) if s is True]
                    list_of_optimized_molecules.extend(converged)
                    not_converged = [n for n, s in zip(not_converged, status_list) if s == 'CycleExceeded' and not tabu.broken(n)]
                    not_converged = clustering.remove_similar(not_converged)
                else:
                    aggregator_logger.info("    All molecules are processed")
                    break
            else:
                aggregator_logger.info("    The following molecules are not converged after 10 rounds")
                for n, s in zip(not_converged, status_list):
                    if s == 'CycleExceeded' and not tabu.broken(n):
                        aggregator_logger.info("      ", n.name)
            os.chdir(cwd)
        
        
        if os.path.exists('selected'):
            os.chdir('selected')
            optimized_molecules = [i.name for i in list_of_optimized_molecules]
            job_done = []
            selected_seeds = []
            for i in optimized_molecules:
                for j in os.listdir():
                    if f'job_{i}' == j:
                        job_done.append(j)
                        if f'result_{i}.xyz' == j:
                            selected_seeds.append(i)
                            list_of_optimized_molecules.pop(j)
            os.chdir(cwd)
        else:
            file_manager.make_directories('selected')
        if len(list_of_optimized_molecules) < 2:
            selected_seeds = list_of_optimized_molecules
        else:
            aggregator_logger.info("  Clustering")
            selected_seeds = clustering.choose_geometries(list_of_optimized_molecules, maximum_number_of_seeds=maximum_number_of_seeds)
            # selected_seeds = clustering.choose_geometries(list_of_optimized_molecules, 
            #                                   maximum_number_of_seeds=maximum_number_of_seeds,
            #                                   algorithm=args.clustering_algorithm)

        os.chdir('selected')
        qc_params["opt_threshold"] = 'normal'
        aggregator_logger.info("Optimizing the selected molecules with higher threshold")

        less_than_ideal = []
        for each_file in selected_seeds:
            not_refined = copy.deepcopy(each_file)
            status = optimise(each_file, qc_params)
            if status is True:
                xyz_file = f'job_{each_file.name}/result_{each_file.name}.xyz'
                shutil.copy(xyz_file, '.')
            else:
                selected_seeds.remove(each_file)
                less_than_ideal.append(not_refined)
        if len(selected_seeds) != 0:
            return selected_seeds
        aggregator_logger.info("    The optimization could not be refined, \n    so sending the loosely optimised molecules")
        return less_than_ideal
    
    else:
        # without software specified
        all_orientations = []
        for seed_count, each_seed in enumerate(seeds):
            if check_stop_signal():
                aggregator_logger.info("Function: add_one")
                return
            aggregator_logger.info(f'   Seed: {seed_count}')
            seed_id = "{:03d}".format(seed_count)
            seeds_home = f'seed_{seed_id}'
            if not os.path.exists(seeds_home):
                file_manager.make_directories(seeds_home)
            os.chdir(seeds_home)
            each_seed.mol_to_xyz('seed.xyz')
            monomer.mol_to_xyz('monomer.xyz')
            if len(each_seed) == 1:
                hm_orientations = 1
            mol_id = f'{seed_id}_{aggregate_id}'
            orientations = generate_orientations(grid_on, hm_orientations,
                                                 mol_id, monomer, seed_count,
                                                 seeds, site, tabu_on)
            all_orientations.extend(orientations)
            os.chdir(cwd)
            
        return all_orientations


def check_for_the_finished_jobs_on_restart(list_of_optimized_molecules, cwd):
    os.chdir('selected')
    optimized_molecules = [i.name for i in list_of_optimized_molecules]
    job_done = []
    result = []
    for i in optimized_molecules:
        for j in os.listdir():
            if f'job_{i}' == j:
                job_done.append(j)
                if f'result_{i}.xyz' == j:
                    result.append(i)
                    list_of_optimized_molecules.pop(j)
    os.chdir(cwd)
    return result



def generate_orientations(use_grid, num_orientations, mol_id, monomer,
                          seed_counter, seeds, site, use_tabu):
    aggregator_logger.debug('Making orientations')
    if not all(os.path.exists(f"trial_{i:03d}_{mol_id}.xyz")
               for i in range(num_orientations)):
        yield from tabu.create_trial_geometries(
            mol_id,
            seeds[seed_counter],
            monomer,
            num_orientations,
            use_tabu,
            use_grid,
            site,
        )
        aggregator_logger.debug('Orientations are made.')
    else:
        yield from read_orientations(mol_id, num_orientations)

def update_id(aid, the_monomer):
    """
        aggregate_id = "a_{:03d}_b_{:03d}_c_{:03d}".format(a_n, b_n, c_n)

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


def check_stop_signal():
    if os.path.exists('stop') or os.path.exists('STOP'):
        aggregator_logger.info(f"Found stop file, in {os.getcwd()}")
        return 1

def generate_molecule_from_formula(formula, box_size=None):
    """
    Generates a molecule object from a given formula using a Tabu-based approach within a box.

    Args:
        formula (str): The chemical formula of the molecule to generate.
        box_size (float, optional): The size of the box in Angstroms. If None, it will be auto-generated.

    Returns:
        Molecule: The generated molecule object.
    """

    # Parse the formula into a dictionary of element counts
    element_counts = {}
    for match in re.findall(r'([A-Z][a-z]*)(\d*)', formula):
        element, count = match
        count = int(count) if count else 1
        element_counts[element] = count

    # Create a list of atoms based on the formula
    atoms_list = []
    for element, count in element_counts.items():
        atoms_list.extend([element] * count)

    # Auto-generate box size if not provided
    if box_size is None:
        max_radius = max(atomic_data.vdw_radii[atomic_data.atomic_number[atom]] for atom in atoms_list)
        box_size = 2 * max_radius + 2.0  # Add some buffer space

    # Generate initial random coordinates within the box
    coordinates = np.random.uniform(0, box_size, size=(len(atoms_list), 3))

    # Create the initial molecule object
    molecule = Molecule(atoms_list, coordinates)

    # Use Tabu search to optimize the molecule's geometry
    tabu_options = {
        'number_of_orientations': 100,  # Adjust as needed
        'tabu_on': True,
        'grid_on': True,
        'd_threshold': 0.3,  # Adjust as needed
        'a_threshold': 15.0  # Adjust as needed
    }
    optimized_molecule = tabu.create_composite_molecule(molecule, molecule, tabu_options, d_scale=1.5)

    return optimized_molecule


def aggregate_from_formulas(formulas, aggregate_sizes, hm_orientations, qc_params, maximum_number_of_seeds,
                            first_pathway, number_of_pathways, tabu_on, grid_on, site):
    """
    Generates aggregates from given formulas using a Tabu-based approach within a box.

    Args:
        formulas (list): A list of chemical formulas for the molecules to aggregate.
        aggregate_sizes (list): The number of each molecule in the final aggregate.
        hm_orientations (int): Number of trial orientations.
        qc_params (dict): Parameters for Quantum Chemistry Calculations.
        maximum_number_of_seeds (int): The maximum number of seeds to be selected for the next cycle.
        first_pathway (int): The starting pathway.
        number_of_pathways (int): The number of pathways to explore.
        tabu_on (bool): Toggle the use of Tabu list.
        grid_on (bool): Toggle the use of Grid.
        site (list, optional): Not used now, but needed for create_trial_molecules().
    """

    # Generate molecule objects from formulas
    molecules = [generate_molecule_from_formula(formula) for formula in formulas]

    # Call the original aggregate function with the generated molecules
    aggregate(molecules, aggregate_sizes, hm_orientations, qc_params, maximum_number_of_seeds,
              first_pathway, number_of_pathways, tabu_on, grid_on, site)

def main():
    pass


if __name__ == "__main__":
    main()
