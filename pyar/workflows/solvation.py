"""Solvation workflow entrypoint."""

from __future__ import annotations

import os

from pyar import file_manager
from pyar.workflows._growth import (
    add_one,
    aggregator_logger,
    check_stop_signal,
    _resolve_orientation_count,
    _working_directory,
)

__all__ = ["solvate"]


def solvate(seeds, monomer, aggregate_size, hm_orientations, qc_params, maximum_number_of_seeds, site=None):
    """Add solvent molecules to the current seeds."""
    if check_stop_signal():
        aggregator_logger.info("Function: solvate")
        return StopIteration

    number_of_orientations = _resolve_orientation_count(hm_orientations)

    starting_directory = os.getcwd()
    aggregator_logger.info(f"Starting solvation in {starting_directory}")
    aggregator_logger.info(
        "Solvation config: seeds=%d solvent=%s solvent_count=%d orientations=%s software=%s max_seeds=%s",
        len(seeds),
        monomer.name,
        aggregate_size,
        hm_orientations,
        qc_params.get("software"),
        maximum_number_of_seeds,
    )
    for aggregation_counter in range(2, aggregate_size + 2):
        if len(seeds) == 0:
            aggregator_logger.info("No seeds to process")
            return
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

        if hm_orientations == "auto" and number_of_orientations <= 256:
            number_of_orientations *= 2
