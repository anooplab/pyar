"""Solvation workflow orchestration for PyAR.

This module owns the growth workflow behind ``pyar-cli --solvate``. It
validates the request, creates or resumes the ``solvation/`` run directory,
generates trial orientations for each cycle, evaluates candidates through the
selected backend, selects the surviving seed geometries, persists
``solvation/state.json`` plus cycle snapshots, and returns a structured
:class:`~pyar.workflow_results.SolvationResult`.
"""

from __future__ import annotations

import os

from pyar.solvation import SolvationRequest
from pyar.solvation.request import normalize_connectivity_policy
from pyar.state.solvation import SolvationRunState, SolvationStateError
from pyar import file_manager
from pyar.workflow_results import SolvationResult
from pyar.workflows._growth import (
    add_one,
    aggregator_logger,
    check_stop_signal,
    sampling_configuration,
    stopped_workflow_result,
    _resolve_orientation_count,
    _working_directory,
    workflow_run_directory,
    workflow_state_path,
)

__all__ = ["solvate"]


def _resolve_solvation_connectivity_policy(connectivity_policy):
    """Return the enforced connectivity policy for solvation workflows."""
    return normalize_connectivity_policy(connectivity_policy)


def _build_solvation_result(run_state, root_directory, sampling):
    """Return the final structured result for a solvation run."""
    return SolvationResult(
        workflow="solvation",
        status=run_state.data["status"],
        run_directory=workflow_run_directory(root_directory, "solvation"),
        state_path=workflow_state_path(workflow_run_directory(root_directory, "solvation")),
        selected_paths=tuple(run_state.data.get("final_seeds", [])),
        metadata={
            "completed_cycles": tuple(
                record["cycle"] for record in run_state.data.get("completed_cycles", [])
            ),
            "next_cycle": run_state.data.get("next_cycle"),
            "sampling": run_state.data.get("sampling", sampling),
            "connectivity_policy": "off",
        },
    )


def solvate(
    seeds,
    monomer,
    aggregate_size,
    hm_orientations,
    qc_params,
    maximum_number_of_seeds,
    site=None,
    connectivity_policy="off",
):
    """Run the solvation workflow for the provided seeds and solvent.

    The workflow either resumes an existing ``solvation/`` state tree or
    creates a new one, then iterates over the requested growth cycles while
    recording the seed geometries selected for the next cycle. The returned
    :class:`~pyar.workflow_results.SolvationResult` describes the final run
    state and the output directory.
    """
    if check_stop_signal():
        aggregator_logger.info("Function: solvate")
        return stopped_workflow_result("solvation", os.getcwd())

    solvation_request = SolvationRequest.from_options(
        seeds,
        monomer,
        aggregate_size,
        hm_orientations,
        qc_params,
        maximum_number_of_seeds,
        site,
        connectivity_policy,
    )
    seeds = list(solvation_request.seeds)
    monomer = solvation_request.monomer
    aggregate_size = solvation_request.aggregate_size
    qc_params = dict(solvation_request.backend_parameters)
    maximum_number_of_seeds = solvation_request.maximum_number_of_seeds
    site = None if solvation_request.site is None else list(solvation_request.site)
    resolved_connectivity_policy = solvation_request.connectivity_policy

    number_of_orientations = _resolve_orientation_count(hm_orientations)

    root_directory = os.getcwd()
    request = solvation_request.to_state_dict()
    sampling = sampling_configuration(
        number_of_orientations=number_of_orientations,
        use_angles=monomer.number_of_atoms > 1,
    )
    run_state = SolvationRunState.load(root_directory, request)
    if run_state is None:
        if os.path.isdir("solvation") and any(os.scandir("solvation")):
            raise SolvationStateError(
                "Existing solvation directory has no resumable state; "
                "preserve it and start in a new directory, or remove it explicitly."
            )
        file_manager.make_directories("solvation")
        run_state = SolvationRunState.create(root_directory, request, seeds, sampling=sampling)
    else:
        seeds = run_state.current_seed_molecules()

    starting_directory = os.getcwd()
    aggregator_logger.info(f"Starting solvation in {starting_directory}")
    aggregator_logger.info("Connectivity policy: off (solvation workflow)")
    aggregator_logger.info(
        "Solvation state detected: %s; next cycle=%d",
        "resuming" if run_state.data.get("completed_cycles") else "starting fresh",
        run_state.data["next_cycle"],
    )
    aggregator_logger.info(
        "Solvation config: seeds=%d solvent=%s solvent_count=%d orientations=%s software=%s max_seeds=%s",
        len(seeds),
        monomer.name,
        aggregate_size,
        hm_orientations,
        qc_params.get("software"),
        maximum_number_of_seeds,
    )
    for aggregation_counter in range(run_state.data["next_cycle"], aggregate_size + 2):
        if len(seeds) == 0:
            aggregator_logger.info("No seeds to process")
            run_state.finish("completed_no_seeds")
            return _build_solvation_result(run_state, root_directory, sampling)
        aggregate_id = "{:03d}".format(aggregation_counter)
        aggregate_home = f"aggregate_{aggregate_id}"
        file_manager.make_directories(aggregate_home)
        with _working_directory(aggregate_home):
            aggregator_logger.info("Solvation cycle path: %s", aggregate_home)
            aggregator_logger.info(f"Solvation cycle start: {aggregation_counter}")
            seeds = add_one(
                aggregate_id,
                seeds,
                monomer,
                number_of_orientations,
                qc_params,
                maximum_number_of_seeds,
                site,
                connectivity_policy=resolved_connectivity_policy,
            )
            aggregator_logger.info(f"Solvation cycle completed: {aggregation_counter}")
        run_state.complete_cycle(aggregation_counter, seeds)

        if hm_orientations == "auto" and number_of_orientations <= 256:
            number_of_orientations *= 2
    run_state.finish()
    return _build_solvation_result(run_state, root_directory, sampling)
