"""Aggregate workflow entrypoint."""

from __future__ import annotations

import copy
import os
import string
from collections import OrderedDict

from pyar import file_manager
from pyar.workflows._growth import (
    add_one,
    aggregator_logger,
    check_stop_signal,
    expand_formula_to_aggregate_inputs,
    generate_molecule_from_formula,
    old_path_to_new_path,
    read_old_path,
    select_pathways,
    update_id,
    _finalize_selected_geometries,
    _format_path,
    _resolve_orientation_count,
    _working_directory,
)

__all__ = [
    "aggregate",
    "aggregate_from_formulas",
    "expand_formula_to_aggregate_inputs",
    "generate_molecule_from_formula",
]


def aggregate(
    molecules,
    aggregate_sizes,
    hm_orientations,
    qc_params,
    maximum_number_of_seeds,
    first_pathway,
    number_of_pathways,
    site,
):
    """Run an aggregate or cluster generation workflow."""
    if check_stop_signal():
        aggregator_logger.info("Function: aggregate")
        return StopIteration

    number_of_orientations = _resolve_orientation_count(hm_orientations)

    parent_folder = "aggregates"
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

    with _working_directory(parent_folder):
        starting_directory = os.getcwd()

        if restart:
            aggregator_logger.info(f"Restarting aggregation in {starting_directory}")
        else:
            aggregator_logger.info(f"Starting aggregation in {starting_directory}")
        aggregator_logger.info(
            "Aggregate config: fragments=%d sizes=%s orientations=%s software=%s max_seeds=%s pathways=%s first_pathway=%s",
            len(molecules),
            aggregate_sizes,
            hm_orientations,
            qc_params.get("software"),
            maximum_number_of_seeds,
            number_of_pathways,
            first_pathway,
        )
        aggregator_logger.debug("Aggregate site=%s restart=%s", site, restart)

        seed_names = string.ascii_lowercase
        ag_id = "ag"

        monomers_to_be_added = []
        for seed_molecule, seed_name, size_of_this_seed in zip(
            molecules, seed_names, aggregate_sizes
        ):
            seed_molecule.name = seed_name
            monomers_to_be_added.extend(seed_molecule for _ in range(size_of_this_seed))
            ag_id += f"_{seed_name}_000"

        aggregator_logger.info("Aggregate identifier: %s", ag_id)
        aggregator_logger.info("Aggregate output root: %s/%s", starting_directory, ag_id)

        if len(molecules) == 1:
            pathways_to_calculate = [monomers_to_be_added]
        elif restart:
            pathways_to_calculate = old_path_to_new_path(monomers_to_be_added, old_path)
        else:
            pathways_to_calculate = select_pathways(
                monomers_to_be_added,
                number_of_pathways,
            )
            aggregator_logger.info("Selected build paths:")
            for i, path in enumerate(pathways_to_calculate):
                aggregator_logger.info("  %03d: %s", i, _format_path(path))

        seed_storage = OrderedDict()
        initial_storage = copy.deepcopy(seed_storage)
        initial_aggregate_id = ag_id
        outside_counter = first_pathway

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
                with _working_directory(ag_home):
                    seed_storage[ag_id] = add_one(
                        ag_id,
                        this_seed,
                        this_monomer,
                        number_of_orientations,
                        qc_params,
                        maximum_number_of_seeds,
                        site,
                    )
                if len(seed_storage[ag_id]) == 0:
                    aggregator_logger.info(
                        f"No molecules were found from {ag_id} to continue this pathway."
                    )
                    aggregator_logger.info("Stopping this path.")
                    break
                seed_storage.popitem(last=False)
            outside_counter += 1
            seed_storage = copy.copy(initial_storage)
            ag_id = initial_aggregate_id

        final_selected = _finalize_selected_geometries(
            aggregate_root=".",
            maximum_number_of_seeds=maximum_number_of_seeds,
            algorithm="hybrid",
        )
        if final_selected:
            aggregator_logger.info(
                "Final cross-path selected geometries written: %d",
                len(final_selected),
            )

        if hm_orientations == "auto" and number_of_orientations <= 256:
            number_of_orientations += 8
        aggregator_logger.info("Aggregation workflow completed.")
        return


def aggregate_from_formulas(
    formulas,
    aggregate_sizes,
    hm_orientations,
    qc_params,
    maximum_number_of_seeds,
    first_pathway,
    number_of_pathways,
    site,
):
    """Generate initial molecules from formulas and run the aggregate workflow."""
    molecules = [generate_molecule_from_formula(formula) for formula in formulas]
    return aggregate(
        molecules,
        aggregate_sizes,
        hm_orientations,
        qc_params,
        maximum_number_of_seeds,
        first_pathway,
        number_of_pathways,
        site,
    )
