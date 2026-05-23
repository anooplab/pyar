import copy
import logging
import os
import random
import shutil
import string
from collections import Counter, OrderedDict
import numpy as np
from pyar import tabu, file_manager
from pyar.Molecule import Molecule
from pyar.data_analysis import clustering
from pyar.optimiser import is_cycle_exceeded, is_success, optimise
import re
from pyar.Molecule import atomic_data
aggregator_logger = logging.getLogger('pyar.aggregator')


def _format_path(path):
    return ''.join(m.name for m in path)


def _stoichiometry_label(molecule):
    """Return a compact chemical formula label for a molecule."""
    counts = Counter(molecule.atoms_list)
    parts = []

    if 'C' in counts:
        carbon = counts.pop('C')
        parts.append('C' if carbon == 1 else f'C{carbon}')
    if 'H' in counts:
        hydrogen = counts.pop('H')
        parts.append('H' if hydrogen == 1 else f'H{hydrogen}')

    for element in sorted(counts):
        count = counts[element]
        parts.append(element if count == 1 else f'{element}{count}')

    return ''.join(parts) if parts else 'unknown'


def _snapshot_selected_geometries(selected_seeds, output_root='selected'):
    """Persist the selected geometries in a stoichiometry-specific folder."""
    if not selected_seeds:
        return None

    stoichiometry = _stoichiometry_label(selected_seeds[0])
    snapshot_dir = os.path.join(output_root, f'stoichiometry_{stoichiometry}')
    if not os.path.exists(snapshot_dir):
        file_manager.make_directories(snapshot_dir)

    for molecule in selected_seeds:
        target_file = os.path.join(snapshot_dir, f'result_{molecule.name}.xyz')
        if hasattr(molecule, 'mol_to_xyz'):
            molecule.mol_to_xyz(target_file)
        else:
            with open(target_file, 'w') as fp:
                fp.write(f"{len(molecule.atoms_list):3d}\n")
                fp.write(f"{getattr(molecule, 'title', molecule.name)}: {getattr(molecule, 'energy', 0.0)}\n")
                for element_symbol, atom_coordinate in zip(molecule.atoms_list, molecule.coordinates):
                    fp.write(("%-2s%12.5f%12.5f%12.5f\n" % (
                        element_symbol, atom_coordinate[0], atom_coordinate[1], atom_coordinate[2]
                    )))

    aggregator_logger.info(
        "Saved %d selected geometries to %s",
        len(selected_seeds),
        snapshot_dir,
    )
    return snapshot_dir


def _supports_two_layer_optimization(qc_params):
    """Return True when the backend supports staged loose/normal optimization.

    The CLI now passes an explicit backend capability flag. Fall back to the
    historical xTB-family behavior only when older call sites do not provide it.
    """
    if '_two_layer_optimization' in qc_params:
        return bool(qc_params.get('_two_layer_optimization'))
    return qc_params.get('software') in {'xtb', 'xtb-aiqm1', 'xtb-aimnet2'}


def random_permutation(iterable, r=None):
    """Random selection from itertools.permutations(iterable, r)"""
    pool = tuple(iterable)
    r = len(pool) if r is None else r
    return tuple(random.sample(pool, r))


def read_old_path():
    """Read restart pathway markers from a prior log, if one exists."""
    if not os.path.exists('pyar.log'):
        return None

    path = []
    with open('pyar.log') as log_file:
        lines = log_file.readlines()
    needed_elements = 3 if any('DEBUG' in line for line in lines) else 2
    for line in lines:
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
    """Run an aggregate or cluster generation workflow.

    The workflow starts from one or more molecules, generates trial aggregate
    orientations, optimizes selected structures, and repeats pathway exploration
    until the requested aggregate sizes are reached.

    Parameters are the input molecules, target aggregate sizes, orientation
    count, quantum-chemistry settings, seed-selection limit, pathway bounds,
    and toggles for the Tabu and grid-based search steps.
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
    parent_folder_exists = os.path.isdir(parent_folder)
    restart = bool(old_path and parent_folder_exists)

    if not parent_folder_exists:
        file_manager.make_directories(parent_folder)
        if old_path:
            aggregator_logger.info(
                "Ignoring restart markers because %s was missing; starting a new aggregation.",
                parent_folder,
            )

    os.chdir(parent_folder)
    starting_directory = os.getcwd()

    if restart:
        aggregator_logger.info(f"Restarting aggregation in {starting_directory}")
    else:
        aggregator_logger.info(f"Starting aggregation in {starting_directory}")
    aggregator_logger.info(
        "Aggregate config: fragments=%d sizes=%s orientations=%s software=%s max_seeds=%s pathways=%s first_pathway=%s",
        len(molecules), aggregate_sizes, hm_orientations, qc_params.get('software'),
        maximum_number_of_seeds, number_of_pathways, first_pathway,
    )
    aggregator_logger.debug(
        "Aggregate flags: tabu=%s grid=%s site=%s restart=%s",
        tabu_on, grid_on, site, restart,
    )

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

        aggregator_logger.info("Selected build paths:")
        for i, path in enumerate(pathways_to_calculate):
            aggregator_logger.info("  %03d: %s", i, _format_path(path))

    seed_storage = OrderedDict()

    initial_storage = copy.deepcopy(seed_storage)
    initial_aggregate_id = ag_id

    outside_counter = first_pathway
    inside_counter = 1

    for i in pathways_to_calculate:
        aggregator_logger.info("Path start: %s", _format_path(i))
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
                aggregator_logger.info('Stopping this path.')
                break
            seed_storage.popitem(last=False)
            inside_counter += 1
        outside_counter += 1
        seed_storage = copy.copy(initial_storage)
        ag_id = initial_aggregate_id

    if hm_orientations == 'auto' and number_of_orientations <= 256:
        number_of_orientations += 8
    aggregator_logger.info("Aggregation workflow completed.")
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
    aggregator_logger.info(f"Starting solvation in {starting_directory}")
    aggregator_logger.info(
        "Solvation config: seeds=%d solvent=%s solvent_count=%d orientations=%s software=%s max_seeds=%s",
        len(seeds), monomer.name, aggregate_size, hm_orientations, qc_params.get('software'),
        maximum_number_of_seeds,
    )
    for aggregation_counter in range(2, aggregate_size + 2):
        if len(seeds) == 0:
            aggregator_logger.info("No seeds to process")
            return
        aggregate_id = "{:03d}".format(aggregation_counter)
        aggregate_home = f'aggregate_{aggregate_id}'
        file_manager.make_directories(aggregate_home)
        os.chdir(aggregate_home)

        aggregator_logger.info(f"Solvation cycle start: {aggregation_counter}")

        seeds = add_one(aggregate_id, seeds, monomer, number_of_orientations,
                        qc_params, maximum_number_of_seeds, tabu_on, grid_on,
                        site)

        aggregator_logger.info(f"Solvation cycle completed: {aggregation_counter}")

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

def add_one(aggregate_id, seeds, monomer, hm_orientations, qc_params, maximum_number_of_seeds, tabu_on, grid_on, site):
    if check_stop_signal():
        aggregator_logger.info("Function: add_one")
        return StopIteration
    aggregator_logger.info(f'  Seed pool for {aggregate_id}: {len(seeds)} molecules')

    cwd = os.getcwd()

    if qc_params.get('software'):
        list_of_optimized_molecules = []
        staged_optimization = _supports_two_layer_optimization(qc_params)
        aggregator_logger.info(
            "Optimization layers: %s",
            "two-stage (loose then normal)" if staged_optimization else "single-stage",
        )
        total_seeds = len(seeds)
        for seed_count, each_seed in enumerate(seeds):
            if check_stop_signal():
                aggregator_logger.info("Function: add_one")
                return
            aggregator_logger.info(
                '   Processing seed %03d/%03d',
                seed_count + 1,
                total_seeds,
            )
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
                    aggregator_logger.info(f"    Optimization round {i + 1:d}: pending={len(not_converged):d}")
                    if staged_optimization:
                        qc_params["opt_threshold"] = 'loose'
                    status_list = [optimise(each_mol, qc_params) for each_mol in not_converged]
                    converged = [n for n, s in zip(not_converged, status_list) if is_success(s)]
                    list_of_optimized_molecules.extend(converged)
                    not_converged = [n for n, s in zip(not_converged, status_list) if is_cycle_exceeded(s) and not tabu.broken(n)]
                    not_converged = clustering.remove_similar(not_converged)
                else:
                    aggregator_logger.info("    All trial molecules processed for this seed")
                    break
            else:
                aggregator_logger.info("    Molecules still unconverged after 10 rounds:")
                for n, s in zip(not_converged, status_list):
                    if is_cycle_exceeded(s) and not tabu.broken(n):
                        aggregator_logger.info("      %s", n.name)
            os.chdir(cwd)

        selected_from_restart = []
        if os.path.exists('selected'):
            selected_from_restart = check_for_the_finished_jobs_on_restart(
                list_of_optimized_molecules,
                cwd,
            )
            if selected_from_restart:
                aggregator_logger.info(
                    "Found %d already refined selected geometries from restart.",
                    len(selected_from_restart),
                )
        else:
            file_manager.make_directories('selected')
        aggregator_logger.info("  Selecting optimized pool")
        selected_seeds = clustering.choose_geometries(
            list_of_optimized_molecules,
            maximum_number_of_seeds=maximum_number_of_seeds,
            persist_basin_memory=not staged_optimization,
        )
        selected_seeds = selected_from_restart + selected_seeds
        if len(selected_seeds) > maximum_number_of_seeds:
            aggregator_logger.warning(
                "Restart plus current selection produced %d seeds; keeping lowest-energy %d.",
                len(selected_seeds),
                maximum_number_of_seeds,
            )
            selected_seeds = sorted(
                selected_seeds,
                key=lambda molecule: float(molecule.energy),
            )[:maximum_number_of_seeds]

        if not staged_optimization:
            aggregator_logger.info(
                "Backend does not support staged optimization; returning selected seeds after a single pass."
            )
            _snapshot_selected_geometries(selected_seeds)
            return selected_seeds

        os.chdir('selected')
        qc_params["opt_threshold"] = 'normal'
        aggregator_logger.info("Refining selected geometries with normal threshold")

        less_than_ideal = []
        refined_seeds = []
        for each_file in selected_seeds:
            not_refined = copy.deepcopy(each_file)
            status = optimise(each_file, qc_params)
            if is_success(status):
                xyz_file = f'job_{each_file.name}/result_{each_file.name}.xyz'
                shutil.copy(xyz_file, '.')
                refined_seeds.append(each_file)
            else:
                less_than_ideal.append(not_refined)
        selected_seeds = refined_seeds
        if len(selected_seeds) != 0:
            aggregator_logger.info("Selection result: %d refined molecules", len(selected_seeds))
            clustering.record_selected_basins(selected_seeds, output_root='.')
            _snapshot_selected_geometries(selected_seeds, output_root='.')
            return selected_seeds
        aggregator_logger.info("Selection result: no refined molecules, returning loose set (%d)", len(less_than_ideal))
        clustering.record_selected_basins(less_than_ideal, output_root='.')
        _snapshot_selected_geometries(less_than_ideal, output_root='.')
        return less_than_ideal
    
    else:
        # without software specified
        all_orientations = []
        total_seeds = len(seeds)
        for seed_count, each_seed in enumerate(seeds):
            if check_stop_signal():
                aggregator_logger.info("Function: add_one")
                return
            aggregator_logger.info(
                '   Processing seed %03d/%03d',
                seed_count + 1,
                total_seeds,
            )
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

        aggregator_logger.info("Generated %d trial geometries without QM optimization", len(all_orientations))
        return all_orientations


def check_for_the_finished_jobs_on_restart(list_of_optimized_molecules, cwd):
    try:
        os.chdir('selected')
        selected_names = set(os.listdir())
        result = []
        remaining = []
        for molecule in list_of_optimized_molecules:
            job_name = f'job_{molecule.name}'
            result_name = f'result_{molecule.name}.xyz'
            if job_name in selected_names and result_name in selected_names:
                result.append(molecule)
            else:
                remaining.append(molecule)
        list_of_optimized_molecules[:] = remaining
        return result
    finally:
        os.chdir(cwd)


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

_FORMULA_SYMBOLS = set(atomic_data.atomic_number)


def _normalize_symbol(symbol):
    """Normalize an element token to standard capitalization."""
    if len(symbol) == 1:
        return symbol.upper()
    return symbol[0].upper() + symbol[1:].lower()


def _read_formula_symbol(formula, position):
    """Read one element symbol from a formula string."""
    char = formula[position]
    if not char.isalpha():
        raise ValueError(f"Invalid chemical formula: {formula!r}")

    next_char = formula[position + 1] if position + 1 < len(formula) else ""
    if char.isupper() and next_char.islower():
        symbol = _normalize_symbol(formula[position:position + 2])
        if symbol in _FORMULA_SYMBOLS:
            return symbol, position + 2
        raise ValueError(f"Unknown element symbol: {formula[position:position + 2]}")

    if char.islower() and next_char.islower():
        symbol = _normalize_symbol(formula[position:position + 2])
        if symbol in _FORMULA_SYMBOLS:
            return symbol, position + 2

    symbol = _normalize_symbol(char)
    if symbol in _FORMULA_SYMBOLS:
        return symbol, position + 1

    if next_char and next_char.islower():
        raise ValueError(f"Unknown element symbol: {formula[position:position + 2]}")
    raise ValueError(f"Unknown element symbol: {char}")


def _parse_formula(formula):
    """Expand a chemical formula string into a list of atom symbols."""
    formula = formula.strip()
    atoms_list = []
    position = 0
    while position < len(formula):
        if formula[position].isspace():
            position += 1
            continue

        symbol, position = _read_formula_symbol(formula, position)

        digit_start = position
        while position < len(formula) and formula[position].isdigit():
            position += 1
        count = int(formula[digit_start:position]) if position > digit_start else 1
        atoms_list.extend([symbol] * count)

    if not atoms_list:
        raise ValueError(f"Invalid chemical formula: {formula!r}")
    return atoms_list


def expand_formula_to_aggregate_inputs(formula):
    """Convert a formula into aggregate input symbols and multiplicities."""
    atoms_list = _parse_formula(formula)
    fragments = []
    aggregate_sizes = []
    for atom in atoms_list:
        if fragments and fragments[-1] == atom:
            aggregate_sizes[-1] += 1
        else:
            fragments.append(atom)
            aggregate_sizes.append(1)
    return fragments, aggregate_sizes


def _atom_radius(atom):
    """Return a packing radius for an element symbol."""
    radius = atomic_data.vdw_radius.get(atom)
    if radius is None or np.isnan(radius):
        radius = atomic_data.covalent_radius[atom]
    return float(radius)


def _estimate_formula_box_size(atoms_list, box_size=None):
    """Estimate a coordinate box large enough for an initial formula geometry."""
    if box_size is not None:
        if box_size <= 0:
            raise ValueError("box_size must be positive")
        return float(box_size)

    max_radius = max(_atom_radius(atom) for atom in atoms_list)
    volume_scale = len(atoms_list) ** (1.0 / 3.0)
    return max(2.0 * max_radius + 2.0, 2.5 * volume_scale * max_radius + 2.0)


def _pack_formula_coordinates(atoms_list, box_size, rng=None, minimum_distance_factor=0.85,
                              max_attempts=500):
    """Pack atoms into a box while avoiding obvious overlaps."""
    rng = np.random.default_rng() if rng is None else rng
    half_box = box_size / 2.0
    coordinates = np.zeros((len(atoms_list), 3), dtype=float)
    placed_indices = []
    placement_order = sorted(range(len(atoms_list)),
                             key=lambda index: (-_atom_radius(atoms_list[index]), index))

    for index in placement_order:
        atom = atoms_list[index]
        attempts = 0
        while True:
            candidate = rng.uniform(-half_box, half_box, size=3)
            if all(
                np.linalg.norm(candidate - coordinates[other_index]) >=
                minimum_distance_factor * (_atom_radius(atom) + _atom_radius(atoms_list[other_index]))
                for other_index in placed_indices
            ):
                coordinates[index] = candidate
                placed_indices.append(index)
                break

            attempts += 1
            if attempts >= max_attempts:
                box_size *= 1.15
                half_box = box_size / 2.0
                attempts = 0
                aggregator_logger.debug(
                    "Expanded formula packing box to %s for atom %s", box_size, atom
                )

    return coordinates


def generate_molecule_from_formula(formula, box_size=None, rng=None):
    """
    Generate a molecule object from a chemical formula.

    Args:
        formula (str): The chemical formula of the molecule to generate.
        box_size (float, optional): The size of the box in Angstroms. If None, it will be auto-generated.
        rng (numpy.random.Generator, optional): Random generator used for packing.

    Returns:
        Molecule: The generated molecule object.
    """
    atoms_list = _parse_formula(formula)
    box_size = _estimate_formula_box_size(atoms_list, box_size=box_size)
    coordinates = _pack_formula_coordinates(atoms_list, box_size, rng=rng)
    return Molecule(atoms_list, coordinates, name=formula, title=formula)


def aggregate_from_formulas(formulas, aggregate_sizes, hm_orientations, qc_params, maximum_number_of_seeds,
                            first_pathway, number_of_pathways, tabu_on, grid_on, site):
    """
    Generate initial molecules from formulas and run the aggregate workflow.

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

    molecules = [generate_molecule_from_formula(formula) for formula in formulas]

    return aggregate(molecules, aggregate_sizes, hm_orientations, qc_params, maximum_number_of_seeds,
                     first_pathway, number_of_pathways, tabu_on, grid_on, site)

def main():
    """Module entry point reserved for future aggregate-specific commands."""
    return None


if __name__ == "__main__":
    main()
