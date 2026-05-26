"""Shared growth-workflow helpers for aggregate and solvation paths."""

from __future__ import annotations

import copy
from contextlib import contextmanager
import logging
import os
import random
import shutil
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

from pyar import file_manager, molecule_io, trial_generation
from pyar.Molecule import Molecule, atomic_data
from pyar.data_analysis import clustering
from pyar.optimiser import is_cycle_exceeded, is_success, optimise

aggregator_logger = logging.getLogger("pyar.aggregator")


@contextmanager
def _working_directory(path):
    """Temporarily change into a working directory and restore the caller cwd."""
    previous_cwd = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(previous_cwd)


def _format_path(path):
    return "".join(m.name for m in path)


def _stoichiometry_label(molecule):
    """Return a compact chemical formula label for a molecule."""
    counts = Counter(molecule.atoms_list)
    parts = []

    if "C" in counts:
        carbon = counts.pop("C")
        parts.append("C" if carbon == 1 else f"C{carbon}")
    if "H" in counts:
        hydrogen = counts.pop("H")
        parts.append("H" if hydrogen == 1 else f"H{hydrogen}")

    for element in sorted(counts):
        count = counts[element]
        parts.append(element if count == 1 else f"{element}{count}")

    return "".join(parts) if parts else "unknown"


def _snapshot_selected_geometries(
    selected_seeds,
    output_root="selected",
    summary_lines=None,
    group_by_stoichiometry=False,
):
    """Persist the selected geometries to disk."""
    if not selected_seeds:
        return None

    snapshot_dir = output_root
    if group_by_stoichiometry:
        stoichiometry = _stoichiometry_label(selected_seeds[0])
        snapshot_dir = os.path.join(output_root, f"stoichiometry_{stoichiometry}")
    else:
        stoichiometry = _stoichiometry_label(selected_seeds[0])
    if not os.path.exists(snapshot_dir):
        file_manager.make_directories(snapshot_dir)
    else:
        for old_result in Path(snapshot_dir).glob("result_*.xyz"):
            old_result.unlink()

    for molecule in selected_seeds:
        target_file = os.path.join(snapshot_dir, f"result_{molecule.name}.xyz")
        molecule_io.write_xyz(molecule, target_file)

    if summary_lines:
        summary_path = os.path.join(snapshot_dir, "README.txt")
        with open(summary_path, "w") as fp:
            fp.write(f"Selected geometries for stoichiometry {stoichiometry}\n")
            fp.write("=" * 60 + "\n")
            for line in summary_lines:
                fp.write(f"{line}\n")

    aggregator_logger.info(
        "Saved %d selected geometries to %s",
        len(selected_seeds),
        snapshot_dir,
    )
    return snapshot_dir


def _resolve_orientation_count(hm_orientations):
    """Return the integer orientation count requested by the caller."""
    if hm_orientations == "auto":
        return 8
    return int(hm_orientations)


def _seed_specific_orientation_count(seed, requested_count):
    """Return the trial count for a seed, honoring monoatomic special cases."""
    return 1 if len(seed) == 1 else requested_count


def _discover_selected_result_files(aggregate_root):
    """Return result files collected from all pathway-level selected folders."""
    root_path = Path(aggregate_root)
    result_files = set()
    for ag_dir in root_path.glob("ag_*"):
        if not ag_dir.is_dir():
            continue
        result_files.update(ag_dir.rglob("selected/result_*.xyz"))
    return sorted(result_files)


def _selection_name_from_path(result_file, aggregate_root):
    """Build a filesystem-safe name from a collected result path."""
    result_path = Path(result_file)
    try:
        relative_path = result_path.relative_to(Path(aggregate_root))
    except ValueError:
        relative_path = result_path
    return str(relative_path.with_suffix("")).replace(os.sep, "_")


def _finalize_selected_geometries(aggregate_root=".", maximum_number_of_seeds=12, algorithm="hybrid"):
    """Cluster pathway-level selected results into the final stoichiometry groups."""
    result_files = _discover_selected_result_files(aggregate_root)
    if not result_files:
        aggregator_logger.info("No pathway-level selected results were found for final clustering.")
        return []

    grouped = defaultdict(list)
    for each_file in result_files:
        molecule = Molecule.from_xyz(str(each_file))
        molecule.energy = clustering.read_energy_from_xyz_file(str(each_file))
        molecule.name = _selection_name_from_path(each_file, aggregate_root)
        grouped[_stoichiometry_label(molecule)].append(molecule)

    final_selected = []
    selected_root = os.path.join(aggregate_root, "selected")
    if not os.path.exists(selected_root):
        file_manager.make_directories(selected_root)
    for stoichiometry, molecules in sorted(grouped.items()):
        aggregator_logger.info(
            "Final clustering for stoichiometry %s from %d pathway results",
            stoichiometry,
            len(molecules),
        )
        selected = clustering.choose_geometries(
            molecules,
            maximum_number_of_seeds=maximum_number_of_seeds,
            apply_basin_memory=False,
            algorithm=algorithm,
        )
        final_selected.extend(selected)
        _snapshot_selected_geometries(
            selected,
            output_root=selected_root,
            summary_lines=[
                f"Stoichiometry: {stoichiometry}",
                "Selection mode: final cross-path clustering",
                f"Source geometries: {len(molecules)}",
                f"Selected geometries: {len(selected)}",
                f"Algorithm: {algorithm}",
            ],
            group_by_stoichiometry=True,
        )

    return final_selected


def _supports_two_layer_optimization(qc_params):
    """Return True when the backend supports staged loose/normal optimization."""
    if "_two_layer_optimization" in qc_params:
        return bool(qc_params.get("_two_layer_optimization"))
    return qc_params.get("software") in {"xtb", "xtb-aiqm1", "xtb-aimnet2"}


def random_permutation(iterable, r=None):
    """Random selection from itertools.permutations(iterable, r)."""
    pool = tuple(iterable)
    r = len(pool) if r is None else r
    return tuple(random.sample(pool, r))


def read_old_path():
    """Read restart pathway markers from a prior log, if one exists."""
    if not os.path.exists("pyar.log"):
        return None

    path = []
    with open("pyar.log") as log_file:
        lines = log_file.readlines()
    needed_elements = 3 if any("DEBUG" in line for line in lines) else 2
    for line in lines:
        split_line = line.split(":")
        if len(split_line) == needed_elements and split_line[-2].strip().isnumeric():
            path.append(split_line[-1].rstrip())
    return path or None


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


def read_orientations(molecule_id, noo):
    orientations = []
    for i in range(noo):
        xyz_file = f"trial_{i:03d}_{molecule_id}.xyz"
        new_orientation = Molecule.from_xyz(xyz_file)
        new_orientation.name = f"{i:03d}_{molecule_id}"
        orientations.append(new_orientation)
    return orientations


def add_one(aggregate_id, seeds, monomer, hm_orientations, qc_params, maximum_number_of_seeds, site):
    if check_stop_signal():
        aggregator_logger.info("Function: add_one")
        return StopIteration
    aggregator_logger.info(f"  Seed pool for {aggregate_id}: {len(seeds)} molecules")

    cwd = os.getcwd()

    if qc_params.get("software"):
        list_of_optimized_molecules = []
        staged_optimization = _supports_two_layer_optimization(qc_params)
        preselection_qc_params = dict(qc_params)
        if staged_optimization:
            preselection_qc_params["opt_threshold"] = "loose"
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
                "   Processing seed %03d/%03d",
                seed_count + 1,
                total_seeds,
            )
            seed_id = "{:03d}".format(seed_count)
            seeds_home = f"seed_{seed_id}"
            if not os.path.exists(seeds_home):
                file_manager.make_directories(seeds_home)
            with _working_directory(seeds_home):
                each_seed.mol_to_xyz("seed.xyz")
                monomer.mol_to_xyz("monomer.xyz")
                seed_orientations = _seed_specific_orientation_count(each_seed, hm_orientations)
                mol_id = f"{seed_id}_{aggregate_id}"
                aggregator_logger.debug("Making orientations")
                if not all(os.path.exists(f"trial_{i:03d}_{mol_id}.xyz") for i in range(seed_orientations)):
                    all_orientations = trial_generation.create_trial_geometries(
                        mol_id, seeds[seed_count], monomer, seed_orientations, site
                    )
                    aggregator_logger.debug("Orientations are made.")
                else:
                    all_orientations = read_orientations(mol_id, seed_orientations)
                not_converged = all_orientations[:]
                status_list = [False for _ in not_converged]
                for i in range(10):
                    if len(not_converged) > 0:
                        aggregator_logger.info(
                            f"    Optimization round {i + 1:d}: pending={len(not_converged):d}"
                        )
                        status_list = [
                            optimise(each_mol, preselection_qc_params)
                            for each_mol in not_converged
                        ]
                        converged = [n for n, s in zip(not_converged, status_list) if is_success(s)]
                        list_of_optimized_molecules.extend(converged)
                        not_converged = [
                            n
                            for n, s in zip(not_converged, status_list)
                            if is_cycle_exceeded(s) and not trial_generation.broken(n)
                        ]
                        not_converged = clustering.remove_similar(not_converged)
                    else:
                        aggregator_logger.info("    All trial molecules processed for this seed")
                        break
                else:
                    aggregator_logger.info("    Molecules still unconverged after 10 rounds:")
                    for n, s in zip(not_converged, status_list):
                        if is_cycle_exceeded(s) and not trial_generation.broken(n):
                            aggregator_logger.info("      %s", n.name)

        selected_from_restart = []
        if os.path.exists("selected"):
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
            file_manager.make_directories("selected")
        aggregator_logger.info("  Selecting optimized pool")
        selected_seeds = clustering.choose_geometries(
            list_of_optimized_molecules,
            maximum_number_of_seeds=maximum_number_of_seeds,
            persist_basin_memory=not staged_optimization,
            group_basin_by_stoichiometry=False,
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
            _snapshot_selected_geometries(
                selected_seeds,
                summary_lines=[
                    "Selection mode: single-stage",
                    f"Backend: {qc_params.get('software')}",
                    f"Selected geometries: {len(selected_seeds)}",
                ],
                group_by_stoichiometry=False,
            )
            return selected_seeds

        refinement_qc_params = dict(qc_params)
        refinement_qc_params["opt_threshold"] = "normal"
        with _working_directory("selected"):
            aggregator_logger.info("Refining selected geometries with normal threshold")

            less_than_ideal = []
            refined_seeds = []
            for each_file in selected_seeds:
                not_refined = copy.deepcopy(each_file)
                status = optimise(each_file, refinement_qc_params)
                if is_success(status):
                    xyz_file = f"job_{each_file.name}/result_{each_file.name}.xyz"
                    shutil.copy(xyz_file, ".")
                    refined_seeds.append(each_file)
                else:
                    less_than_ideal.append(not_refined)
            selected_seeds = refined_seeds
            if len(selected_seeds) != 0:
                aggregator_logger.info("Selection result: %d refined molecules", len(selected_seeds))
                _snapshot_selected_geometries(
                    selected_seeds,
                    output_root=".",
                    summary_lines=[
                        "Selection mode: two-stage",
                        "Stage 1: loose preselection",
                        "Stage 2: normal-threshold refinement",
                        f"Backend: {qc_params.get('software')}",
                        f"Selected geometries: {len(selected_seeds)}",
                    ],
                    group_by_stoichiometry=False,
                )
                return selected_seeds
            aggregator_logger.info(
                "Selection result: no refined molecules, returning loose set (%d)",
                len(less_than_ideal),
            )
            _snapshot_selected_geometries(
                less_than_ideal,
                output_root=".",
                summary_lines=[
                    "Selection mode: two-stage",
                    "Stage 1: loose preselection",
                    "Stage 2: normal-threshold refinement failed for all selected seeds",
                    f"Backend: {qc_params.get('software')}",
                    f"Selected geometries: {len(less_than_ideal)}",
                ],
                group_by_stoichiometry=False,
            )
            return less_than_ideal

    # without software specified
    # TODO: geometry-only growth currently carries every trial orientation
    # forward to the next stoichiometry because there is no energy-based
    # seed selection step in this branch. Add a bounded diversity-based
    # selection path here so large builds do not expand combinatorially.
    all_orientations = []
    total_seeds = len(seeds)
    for seed_count, each_seed in enumerate(seeds):
        if check_stop_signal():
            aggregator_logger.info("Function: add_one")
            return
        aggregator_logger.info(
            "   Processing seed %03d/%03d",
            seed_count + 1,
            total_seeds,
        )
        seed_id = "{:03d}".format(seed_count)
        seeds_home = f"seed_{seed_id}"
        if not os.path.exists(seeds_home):
            file_manager.make_directories(seeds_home)
        with _working_directory(seeds_home):
            each_seed.mol_to_xyz("seed.xyz")
            monomer.mol_to_xyz("monomer.xyz")
            seed_orientations = _seed_specific_orientation_count(each_seed, hm_orientations)
            mol_id = f"{seed_id}_{aggregate_id}"
            orientations = generate_orientations(
                seed_orientations, mol_id, monomer, seed_count, seeds, site
            )
            all_orientations.extend(orientations)

    aggregator_logger.info(
        "Generated %d trial geometries without QM optimization", len(all_orientations)
    )
    return all_orientations


def check_for_the_finished_jobs_on_restart(list_of_optimized_molecules, cwd):
    try:
        os.chdir("selected")
        selected_names = set(os.listdir())
        result = []
        remaining = []
        for molecule in list_of_optimized_molecules:
            job_name = f"job_{molecule.name}"
            result_name = f"result_{molecule.name}.xyz"
            if job_name in selected_names and result_name in selected_names:
                result.append(molecule)
            else:
                remaining.append(molecule)
        list_of_optimized_molecules[:] = remaining
        return result
    finally:
        os.chdir(cwd)


def generate_orientations(num_orientations, mol_id, monomer, seed_counter, seeds, site):
    aggregator_logger.debug("Making orientations")
    if not all(os.path.exists(f"trial_{i:03d}_{mol_id}.xyz") for i in range(num_orientations)):
        yield from trial_generation.create_trial_geometries(
            mol_id,
            seeds[seed_counter],
            monomer,
            num_orientations,
            site,
        )
        aggregator_logger.debug("Orientations are made.")
    else:
        yield from read_orientations(mol_id, num_orientations)


def update_id(aid, the_monomer):
    """Increment one monomer count inside an aggregate identifier."""
    from collections import deque

    parts_of_aid = deque(aid.split("_"))
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
    if os.path.exists("stop") or os.path.exists("STOP"):
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
        symbol = _normalize_symbol(formula[position : position + 2])
        if symbol in _FORMULA_SYMBOLS:
            return symbol, position + 2
        raise ValueError(f"Unknown element symbol: {formula[position:position + 2]}")

    if char.islower() and next_char.islower():
        symbol = _normalize_symbol(formula[position : position + 2])
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


def _pack_formula_coordinates(atoms_list, box_size, rng=None, minimum_distance_factor=0.85, max_attempts=500):
    """Pack atoms into a box while avoiding obvious overlaps."""
    rng = np.random.default_rng() if rng is None else rng
    half_box = box_size / 2.0
    coordinates = np.zeros((len(atoms_list), 3), dtype=float)
    placed_indices = []
    placement_order = sorted(
        range(len(atoms_list)),
        key=lambda index: (-_atom_radius(atoms_list[index]), index),
    )

    for index in placement_order:
        atom = atoms_list[index]
        attempts = 0
        while True:
            candidate = rng.uniform(-half_box, half_box, size=3)
            if all(
                np.linalg.norm(candidate - coordinates[other_index])
                >= minimum_distance_factor
                * (_atom_radius(atom) + _atom_radius(atoms_list[other_index]))
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
    """Generate a molecule object from a chemical formula."""
    atoms_list = _parse_formula(formula)
    box_size = _estimate_formula_box_size(atoms_list, box_size=box_size)
    coordinates = _pack_formula_coordinates(atoms_list, box_size, rng=rng)
    return Molecule(atoms_list, coordinates, name=formula, title=formula)
