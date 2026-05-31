"""Aggregate workflow orchestration for PyAR.

This module owns the aggregation and cluster-generation workflow used by
``pyar-cli -a``. It is
responsible for:

* validating the aggregate request and restart state
* creating or resuming the ``aggregates/`` run directory
* selecting pathways and orientation sets
* generating trial candidates and running backend optimization
* selecting survivors from pathway-level and final clustering stages
* persisting ``aggregates/state.json`` plus selected-geometry snapshots
* returning a structured :class:`~pyar.workflow_results.AggregateResult`
* running the backend optimization steps for each pathway
* collecting the selected geometries and workflow result metadata

The public entry points are :func:`aggregate` and
:func:`aggregate_from_formulas`.
"""

from __future__ import annotations

import copy
import os
import string
from collections import OrderedDict
from pathlib import Path

import numpy as np

from pyar.state.aggregate import AggregateRunState, AggregateStateError
from pyar import file_manager
from pyar.workflow_results import AggregateResult
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
    sampling_configuration,
    _resolve_orientation_count,
    workflow_run_directory,
    workflow_state_path,
    stopped_workflow_result,
    _working_directory,
)

__all__ = [
    "aggregate",
    "aggregate_from_formulas",
    "expand_formula_to_aggregate_inputs",
    "generate_molecule_from_formula",
]


def _molecule_signature(molecule):
    """Return stable input geometry metadata used to validate restarts."""
    return {
        "atoms": list(molecule.atoms_list),
        "coordinates": np.asarray(molecule.coordinates, dtype=float).tolist(),
        "charge": molecule.charge,
        "multiplicity": molecule.multiplicity,
        "scftype": molecule.scftype,
        "fragment_definition": getattr(molecule, "fragments", []),
    }


def _build_aggregate_request(
    molecules,
    aggregate_sizes,
    hm_orientations,
    qc_params,
    maximum_number_of_seeds,
    first_pathway,
    number_of_pathways,
    site,
):
    """Build the scientific request persisted with aggregation state.

    The persisted request is the restart contract for aggregation. It captures
    the fragment composition, orientation count, backend configuration, and
    pathway-selection controls so a resumed run can reject incompatible input
    changes.
    """
    return {
        "aggregate_sizes": [int(size) for size in aggregate_sizes],
        "orientations": hm_orientations,
        "backend_parameters": dict(qc_params),
        "maximum_number_of_seeds": int(maximum_number_of_seeds),
        "first_pathway": int(first_pathway),
        "number_of_pathways": int(number_of_pathways),
        "site": None if site is None else list(site),
        "fragments": [_molecule_signature(molecule) for molecule in molecules],
    }


def _pathway_from_label(monomers_to_be_added, label):
    """Restore one persisted pathway label to its input molecule objects.

    Pathways are persisted as labels so they can survive restarts and legacy
    checkpoint migration. This helper maps a stored label back to the concrete
    molecule objects used in the current run.
    """
    by_name = {}
    for molecule in monomers_to_be_added:
        by_name.setdefault(molecule.name, molecule)
    try:
        return tuple(by_name[name] for name in label)
    except KeyError as exc:
        raise AggregateStateError(
            f"Saved aggregation pathway {label!r} contains unknown fragment {exc.args[0]!r}."
        ) from exc


def _selected_result_paths(pattern):
    """Return result paths relative to the aggregate workflow directory."""
    return sorted(str(path) for path in Path(".").glob(pattern))


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
    """Run an aggregate or cluster-generation workflow.

    The function creates or resumes the aggregation state, iterates over the
    selected pathways, runs the per-orientation optimization for each path,
    and returns an :class:`~pyar.workflow_results.AggregateResult` describing
    the final run state.
    """
    if check_stop_signal():
        aggregator_logger.info("Function: aggregate")
        return stopped_workflow_result("aggregate", os.getcwd())

    number_of_orientations = _resolve_orientation_count(hm_orientations)
    sampling = sampling_configuration(
        number_of_orientations=number_of_orientations,
        use_angles=None,
    )

    root_directory = os.getcwd()
    parent_folder = "aggregates"
    request = _build_aggregate_request(
        molecules,
        aggregate_sizes,
        hm_orientations,
        qc_params,
        maximum_number_of_seeds,
        first_pathway,
        number_of_pathways,
        site,
    )
    run_state = AggregateRunState.load(root_directory, request)
    old_path = read_old_path()
    parent_folder_exists = os.path.isdir(parent_folder)
    restart = run_state is not None or bool(old_path and parent_folder_exists)

    if not parent_folder_exists:
        file_manager.make_directories(parent_folder)
        if old_path:
            aggregator_logger.info(
                "Ignoring restart markers because %s was missing; starting a new aggregation.",
                parent_folder,
            )
    elif run_state is None and not old_path and any(Path(parent_folder).iterdir()):
        raise AggregateStateError(
            "Existing aggregates directory has no resumable state; "
            "preserve it and start in a new directory, or remove it explicitly."
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

        if run_state is not None:
            aggregator_logger.info(
                "Aggregation state detected: resuming %d remaining pathways",
                len(run_state.remaining_pathway_labels()),
            )
        elif old_path and parent_folder_exists:
            pathways_to_calculate = old_path_to_new_path(monomers_to_be_added, old_path)
            run_state = AggregateRunState.create(
                root_directory,
                request,
                [_format_path(path) for path in pathways_to_calculate],
                legacy_import="pyar.log",
                sampling=sampling,
            )
            aggregator_logger.warning(
                "Imported legacy pyar.log pathway markers into aggregates/state.json."
            )
        else:
            if len(molecules) == 1:
                pathways_to_calculate = [monomers_to_be_added]
            else:
                pathways_to_calculate = select_pathways(
                    monomers_to_be_added,
                    number_of_pathways,
                )
            aggregator_logger.info("Selected build paths:")
            for i, path in enumerate(pathways_to_calculate):
                aggregator_logger.info("  %03d: %s", i, _format_path(path))
            run_state = AggregateRunState.create(
                root_directory,
                request,
                [_format_path(path) for path in pathways_to_calculate],
                sampling=sampling,
            )

        pathways_to_calculate = [
            (label, _pathway_from_label(monomers_to_be_added, label))
            for label in run_state.remaining_pathway_labels()
        ]

        seed_storage = OrderedDict()
        initial_storage = copy.deepcopy(seed_storage)
        initial_aggregate_id = ag_id
        outside_counter = first_pathway + run_state.completed_pathway_count()

        for label, pathway in pathways_to_calculate:
            aggregator_logger.info("Path start: %s", label)
            for this_monomer in pathway:
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
            run_state.complete_pathway(
                label,
                _selected_result_paths(f"ag_*_{outside_counter:03d}/selected/result_*.xyz"),
            )
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
        run_state.finish(_selected_result_paths("selected/**/result_*.xyz"))
        result = AggregateResult(
            workflow="aggregate",
            status=run_state.data["status"],
            run_directory=workflow_run_directory(root_directory, "aggregates"),
            state_path=workflow_state_path(workflow_run_directory(root_directory, "aggregates")),
            selected_paths=tuple(run_state.data.get("final_selected_results", [])),
            metadata={
                "pathways": tuple(run_state.data.get("pathways", [])),
                "completed_pathways": tuple(
                    record["label"] for record in run_state.data.get("completed_pathways", [])
                ),
                "legacy_import": run_state.data.get("legacy_import"),
                "sampling": run_state.data.get("sampling", sampling),
            },
        )

        if hm_orientations == "auto" and number_of_orientations <= 256:
            number_of_orientations += 8
        aggregator_logger.info("Aggregation workflow completed.")
        return result


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
    """Generate initial molecules from formulas and run the aggregate workflow.

    This is the convenience entry point used by ``--formula`` aggregation
    runs. It converts the provided formulas into molecule objects and then
    delegates to :func:`aggregate`.
    """
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
