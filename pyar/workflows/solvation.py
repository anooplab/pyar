"""Solvation workflow orchestration for PyAR.

This module owns the growth workflow behind ``pyar-cli --solvate``. It
persists a restartable solvation state, generates trial orientations for each
cycle, and returns a structured result that describes the final selected
seed geometries.
"""

from __future__ import annotations

import os
from pathlib import Path

from pyar.state.solvation import SolvationRunState, SolvationStateError
from pyar import file_manager
from pyar.workflow_results import SolvationResult
from pyar.workflows._growth import (
    add_one,
    aggregator_logger,
    check_stop_signal,
    _resolve_orientation_count,
    _working_directory,
)

__all__ = ["solvate"]


def _molecule_signature(molecule):
    """Return stable molecule metadata used to validate solvation restarts."""
    return {
        "name": molecule.name,
        "atoms": list(molecule.atoms_list),
        "coordinates": [[float(coord) for coord in row] for row in molecule.coordinates],
        "charge": molecule.charge,
        "multiplicity": molecule.multiplicity,
        "scftype": molecule.scftype,
        "fragment_definition": getattr(molecule, "fragments", []),
    }


def _build_solvation_request(seeds, monomer, aggregate_size, hm_orientations, qc_params, maximum_number_of_seeds, site):
    """Build the persisted request used to validate solvation restarts.

    The request captures the seed geometries, solvent monomer, cycle size,
    orientation count, backend configuration, and site constraint so a future
    invocation can only resume if the scientific inputs are unchanged.
    """
    return {
        "aggregate_size": int(aggregate_size),
        "orientations": hm_orientations,
        "backend_parameters": dict(qc_params),
        "maximum_number_of_seeds": int(maximum_number_of_seeds),
        "site": None if site is None else list(site),
        "seeds": [_molecule_signature(seed) for seed in seeds],
        "monomer": _molecule_signature(monomer),
    }


def solvate(seeds, monomer, aggregate_size, hm_orientations, qc_params, maximum_number_of_seeds, site=None):
    """Run the solvation workflow for the provided seeds and solvent.

    The workflow either resumes an existing ``solvation/`` state tree or
    creates a new one, then iterates over the requested growth cycles while
    recording the seed geometries selected for the next cycle. The returned
    :class:`~pyar.workflow_results.SolvationResult` describes the final run
    state and the output directory.
    """
    if check_stop_signal():
        aggregator_logger.info("Function: solvate")
        run_directory = str(Path.cwd().resolve() / "solvation")
        return SolvationResult(
            workflow="solvation",
            status="stopped",
            run_directory=run_directory,
            state_path=str(Path(run_directory) / "state.json"),
        )

    number_of_orientations = _resolve_orientation_count(hm_orientations)

    root_directory = os.getcwd()
    request = _build_solvation_request(
        seeds,
        monomer,
        aggregate_size,
        hm_orientations,
        qc_params,
        maximum_number_of_seeds,
        site,
    )
    run_state = SolvationRunState.load(root_directory, request)
    if run_state is None:
        if os.path.isdir("solvation") and any(os.scandir("solvation")):
            raise SolvationStateError(
                "Existing solvation directory has no resumable state; "
                "preserve it and start in a new directory, or remove it explicitly."
            )
        file_manager.make_directories("solvation")
        run_state = SolvationRunState.create(root_directory, request, seeds)
    else:
        seeds = run_state.current_seed_molecules()

    starting_directory = os.getcwd()
    aggregator_logger.info(f"Starting solvation in {starting_directory}")
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
            return SolvationResult(
                workflow="solvation",
                status=run_state.data["status"],
                run_directory=str(Path(root_directory).resolve() / "solvation"),
                state_path=str(Path(root_directory).resolve() / "solvation" / "state.json"),
                selected_paths=tuple(run_state.data.get("final_seeds", [])),
                metadata={
                    "completed_cycles": tuple(
                        record["cycle"] for record in run_state.data.get("completed_cycles", [])
                    ),
                    "next_cycle": run_state.data.get("next_cycle"),
                },
            )
        aggregate_id = "{:03d}".format(aggregation_counter)
        aggregate_home = f"aggregate_{aggregate_id}"
        file_manager.make_directories(aggregate_home)
        with _working_directory(aggregate_home):
            aggregator_logger.info(f"Solvation cycle start: {aggregation_counter}")
            seeds = add_one(
                aggregate_id,
                seeds,
                monomer,
                number_of_orientations,
                qc_params,
                maximum_number_of_seeds,
                site,
            )
            aggregator_logger.info(f"Solvation cycle completed: {aggregation_counter}")
        run_state.complete_cycle(aggregation_counter, seeds)

        if hm_orientations == "auto" and number_of_orientations <= 256:
            number_of_orientations *= 2
    run_state.finish()
    return SolvationResult(
        workflow="solvation",
        status=run_state.data["status"],
        run_directory=str(Path(root_directory).resolve() / "solvation"),
        state_path=str(Path(root_directory).resolve() / "solvation" / "state.json"),
        selected_paths=tuple(run_state.data.get("final_seeds", [])),
        metadata={
            "completed_cycles": tuple(
                record["cycle"] for record in run_state.data.get("completed_cycles", [])
            ),
            "next_cycle": run_state.data.get("next_cycle"),
        },
    )
