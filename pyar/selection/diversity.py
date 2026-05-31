"""Diversity selection helpers."""

from __future__ import annotations

import numpy as np

__all__ = [
    "_finalize_selection",
    "_limit_seed_count",
    "_log_seed_shortfall",
    "_max_min_diversity_select",
]


def _log_seed_shortfall(requested, available, context):
    """Log that fewer unique geometries exist than the requested seed count."""
    from pyar.selection import clustering

    clustering.cluster_logger.info(
        "Requested %d seeds, but only %d unique geometries were available after %s.",
        requested,
        available,
        context,
    )


def _finalize_selection(selected_molecules, basin_registry_path, existing_entries=None):
    """Persist basin memory and return the selected molecules."""
    if basin_registry_path:
        from pyar.selection.basin_memory import _persist_basin_registry

        _persist_basin_registry(basin_registry_path, selected_molecules, existing_entries=existing_entries)
    return selected_molecules


def _limit_seed_count(molecules, maximum_number_of_seeds, reason):
    """Enforce a hard upper bound on selected molecules."""
    from pyar.selection import clustering

    if len(molecules) <= maximum_number_of_seeds:
        if len(molecules) < maximum_number_of_seeds:
            _log_seed_shortfall(maximum_number_of_seeds, len(molecules), reason)
        return molecules
    clustering.cluster_logger.warning(
        "Seed count exceeded limit after %s: keeping lowest-energy %d of %d",
        reason,
        maximum_number_of_seeds,
        len(molecules),
    )
    limited = sorted(molecules, key=lambda m: float(m.energy))[:maximum_number_of_seeds]
    print_energy_table(limited)
    return limited


def _max_min_diversity_select(features, molecules, maximum_number_of_seeds, initial_selected_indices=None):
    """Select diverse seeds via greedy max-min (farthest-point) selection.

    When initial_selected_indices is provided, those points are treated as fixed
    anchors and only the additional selected molecules are returned.
    """
    from pyar.selection import clustering

    if initial_selected_indices is None:
        initial_selected_indices = []

    if len(molecules) <= maximum_number_of_seeds and not initial_selected_indices:
        return molecules

    selected_indices = list(initial_selected_indices)
    returned_indices = []
    remaining_indices = set(range(len(molecules))) - set(selected_indices)

    if not selected_indices:
        energies = np.array([float(m.energy) for m in molecules])
        start_index = int(np.argmin(energies))
        selected_indices.append(start_index)
        returned_indices.append(start_index)
        remaining_indices.remove(start_index)

    while len(selected_indices) < maximum_number_of_seeds and remaining_indices:
        best_idx = None
        best_score = -np.inf
        for idx in remaining_indices:
            min_distance = min(
                np.linalg.norm(features[idx] - features[sidx]) for sidx in selected_indices
            )
            if min_distance > best_score and not np.isclose(min_distance, best_score):
                best_score = min_distance
                best_idx = idx
            elif np.isclose(min_distance, best_score) and best_idx is not None:
                if float(molecules[idx].energy) < float(molecules[best_idx].energy):
                    best_idx = idx
        selected_indices.append(best_idx)
        returned_indices.append(best_idx)
        remaining_indices.remove(best_idx)

    selected = [molecules[i] for i in returned_indices]
    if initial_selected_indices:
        clustering.cluster_logger.info(
            "Added %d filler seeds by max-min strategy.",
            len(selected),
        )
    else:
        clustering.cluster_logger.info(
            "Selected %d diverse seeds by max-min strategy.",
            len(selected),
        )
    print_energy_table(selected)
    return selected


def print_energy_table(molecules, stream=None, title=None):
    """Report energies with relative values against the global minimum."""
    from pyar.selection import reports

    return reports.print_energy_table(molecules, stream=stream, title=title)
