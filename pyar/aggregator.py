"""Legacy aggregator module.

This module is retained as a compatibility shim. New code should import the
workflow entrypoints from :mod:`pyar.workflows.aggregate` and
:mod:`pyar.workflows.solvation`.
"""

from __future__ import annotations

import importlib
import warnings

aggregate_workflow = importlib.import_module("pyar.workflows.aggregate")
solvation_workflow = importlib.import_module("pyar.workflows.solvation")
from pyar.workflows._growth import (
    add_one,
    aggregator_logger,
    check_for_the_finished_jobs_on_restart,
    check_stop_signal,
    expand_formula_to_aggregate_inputs,
    generate_molecule_from_formula,
    generate_orientations,
    old_path_to_new_path,
    read_old_path,
    read_orientations,
    select_pathways,
    update_id,
    _finalize_selected_geometries,
    _format_path,
    _resolve_orientation_count,
    _seed_specific_orientation_count,
    _snapshot_selected_geometries,
    _supports_two_layer_optimization,
    _working_directory,
)

__all__ = [
    "aggregate",
    "aggregate_from_formulas",
    "solvate",
    "add_one",
    "aggregator_logger",
    "check_for_the_finished_jobs_on_restart",
    "check_stop_signal",
    "expand_formula_to_aggregate_inputs",
    "generate_molecule_from_formula",
    "generate_orientations",
    "old_path_to_new_path",
    "read_old_path",
    "read_orientations",
    "select_pathways",
    "update_id",
    "_finalize_selected_geometries",
    "_format_path",
    "_resolve_orientation_count",
    "_seed_specific_orientation_count",
    "_snapshot_selected_geometries",
    "_supports_two_layer_optimization",
    "_working_directory",
]


def _warn_legacy(name):
    warnings.warn(
        f"pyar.aggregator.{name} is deprecated; import from pyar.workflows instead.",
        DeprecationWarning,
        stacklevel=2,
    )


def aggregate(*args, **kwargs):
    """Deprecated alias for :func:`pyar.workflows.aggregate.aggregate`."""
    _warn_legacy("aggregate")
    return aggregate_workflow.aggregate(*args, **kwargs)


def aggregate_from_formulas(*args, **kwargs):
    """Deprecated alias for :func:`pyar.workflows.aggregate.aggregate_from_formulas`."""
    _warn_legacy("aggregate_from_formulas")
    return aggregate_workflow.aggregate_from_formulas(*args, **kwargs)


def solvate(*args, **kwargs):
    """Deprecated alias for :func:`pyar.workflows.solvation.solvate`."""
    _warn_legacy("solvate")
    return solvation_workflow.solvate(*args, **kwargs)
