"""Deduplication and connectivity filtering for selected geometries."""

from __future__ import annotations

import numpy as np

import pyar.representations
from pyar.structure_comparison import rmsd as _rmsd

__all__ = [
    "calc_fingerprint_distance",
    "remove_similar",
    "deduplicate_structures",
]


def calc_fingerprint_distance(a, b):
    """Calculate the distance between two fingerprints."""
    return np.linalg.norm(
        pyar.representations.fingerprint(a.atoms_list, a.coordinates)
        - pyar.representations.fingerprint(b.atoms_list, b.coordinates)
    )


def _structure_is_similar(candidate, kept):
    """Quick prefilter using fingerprint distance."""
    # Resolve through the canonical clustering module so tests and extensions
    # can patch one selection seam.
    from pyar.selection import clustering

    fingerprint_distance = clustering.calc_fingerprint_distance(candidate, kept)
    return abs(fingerprint_distance) < 1.0


def _kabsch_rotation(reference, mobile):
    """Compatibility wrapper for the comparison-layer Kabsch primitive."""
    return _rmsd.kabsch_rotation(reference, mobile)


def _kabsch_rmsd(reference, mobile):
    """Compatibility wrapper for the comparison-layer Kabsch RMSD."""
    return _rmsd.kabsch_rmsd(reference, mobile)


def _equivalent_atom_groups(atoms):
    """Compatibility wrapper for element grouping."""
    return _rmsd.equivalent_atom_groups(atoms)


def _exact_element_orders(reference_atoms, mobile_atoms, maximum_permutations=720):
    """Compatibility wrapper for exact element-preserving permutations."""
    return _rmsd.exact_element_orders(
        reference_atoms, mobile_atoms, maximum_permutations=maximum_permutations,
    )


def _assigned_element_order(reference_atoms, reference, mobile_atoms, mobile):
    """Compatibility wrapper for Hungarian element-constrained assignment."""
    return _rmsd.assigned_element_order(reference_atoms, reference, mobile_atoms, mobile)


def _iterative_assigned_rmsd(reference_atoms, reference, mobile_atoms, mobile):
    """Compatibility wrapper for iterative element-assignment RMSD."""
    return _rmsd.iterative_assigned_rmsd(
        reference_atoms, reference, mobile_atoms, mobile,
    )


def _rmsd_after_alignment(candidate, kept):
    """Compatibility wrapper for permutation-aware aligned RMSD."""
    return _rmsd.rmsd_after_alignment(candidate, kept)


def _adaptive_duplicate_rmsd_threshold(molecules):
    """Estimate a duplicate RMSD threshold from the current candidate pool."""
    if len(molecules) < 2:
        return 0.10

    sampled_rmsd = []
    sample_limit = 40
    for left_index, left in enumerate(molecules):
        if len(sampled_rmsd) >= sample_limit:
            break
        for right in molecules[left_index + 1:]:
            if len(sampled_rmsd) >= sample_limit:
                break
            if len(left.atoms_list) != len(right.atoms_list):
                continue
            try:
                if not _structure_is_similar(left, right):
                    continue
                sampled_rmsd.append(_rmsd_after_alignment(left, right))
            except Exception:
                continue

    if not sampled_rmsd:
        return 0.10

    finite_rmsd = [value for value in sampled_rmsd if np.isfinite(value)]
    if not finite_rmsd:
        return 0.10

    baseline = float(np.percentile(finite_rmsd, 5))
    threshold = max(0.05, baseline * 0.75)
    return float(np.clip(threshold, 0.05, 0.15))


def _deduplicate_ordered(ordered_molecules, comparator, *, logger=None, skip_single_atoms=False):
    """Shared conservative elimination loop; diagnostics never imply duplicates."""
    final_list, removed_duplicates, diagnostics = [], [], []

    class Reporter:
        def warning(self, message, *args):
            diagnostics.append(message % args)
            if logger is not None:
                logger.warning(message, *args)

        def debug(self, message, *args):
            if logger is not None:
                logger.debug(message, *args)

    reporter = Reporter()
    for candidate in ordered_molecules:
        duplicate = False
        for kept in final_list:
            if skip_single_atoms and (len(candidate.atoms_list) < 2 or len(kept.atoms_list) < 2):
                continue
            try:
                comparison = comparator.compare(candidate, kept)
            except Exception as exc:
                reporter.warning(
                    "Retaining %s and %s: structural comparison failed (%s: %s).",
                    candidate.name, kept.name, type(exc).__name__, exc,
                )
                continue
            if comparison.equivalent is True:
                aligned_rmsd = comparison.distance
                duplicate = True
                removed_duplicates.append((candidate.name, kept.name, aligned_rmsd))
                reporter.debug(
                    'Removing {} as a near-duplicate of {}'.format(candidate.name, kept.name)
                )
                break
            if (comparison.metadata.get("comparison_complete") is False
                    and comparison.metadata.get("fallback_status") != "ok"):
                reporter.warning(
                    "Retaining %s and %s: graph RMSD was incomplete after %d mappings; iRMSD check status=%s.",
                    candidate.name,
                    kept.name,
                    comparison.metadata.get("isomorphisms_evaluated", 0),
                    comparison.metadata.get("fallback_status", "not_run"),
                )
        if not duplicate:
            final_list.append(candidate)
    return final_list, removed_duplicates, diagnostics


def deduplicate_structures(molecules, *, threshold=None, ordering="auto"):
    """Deduplicate with the existing adaptive threshold and graph-first policy.

    Auto ordering prefers lowest energy only when every energy is available;
    otherwise input order is preserved. Workflows retain their historical
    ordering, single-atom exemption and reporting through ``remove_similar``.
    """
    from pyar.structure_comparison import GraphFirstDeduplicationComparator

    molecules = list(molecules)
    if ordering not in {"auto", "input", "energy"}:
        raise ValueError('ordering must be auto, input, or energy')
    all_energies = all(m.energy is not None and np.isfinite(float(m.energy)) for m in molecules)
    if ordering == "energy" and not all_energies:
        raise ValueError('Energy ordering requires every structure to have an energy')
    if ordering == "energy" or (ordering == "auto" and all_energies):
        molecules = sorted(molecules, key=lambda molecule: float(molecule.energy))
    resolved = _adaptive_duplicate_rmsd_threshold(molecules) if threshold is None else threshold
    comparator = GraphFirstDeduplicationComparator(threshold=resolved)
    kept, removed, diagnostics = _deduplicate_ordered(molecules, comparator)
    return {"kept": kept, "removed": removed, "diagnostics": diagnostics,
            "rmsd_threshold": resolved, "input_count": len(molecules)}


def remove_similar(list_of_molecules):
    """Remove geometrical duplicates under the graph-first comparison policy.

    iRMSD is used only when graph identity matched but graph mapping enumeration
    was incomplete. Any failed, asymmetric-threshold, or diagnostic fallback
    keeps both candidates.
    """
    from pyar.selection import clustering

    ordered_molecules = sorted(list_of_molecules, key=lambda molecule: (float(molecule.energy), molecule.name))
    rmsd_threshold = _adaptive_duplicate_rmsd_threshold(ordered_molecules)
    from pyar.structure_comparison import GraphFirstDeduplicationComparator

    comparator = GraphFirstDeduplicationComparator(threshold=rmsd_threshold)
    clustering.cluster_logger.debug('Number of molecules before similarity elimination,  {}'.format(len(ordered_molecules)))
    final_list, removed_duplicates, _ = _deduplicate_ordered(
        ordered_molecules, comparator, logger=clustering.cluster_logger, skip_single_atoms=True)
    clustering.cluster_logger.debug('Number of molecules after similarity elimination,  {}'.format(len(final_list)))
    if removed_duplicates:
        clustering.cluster_logger.info(
            "Similarity pruning retained %d geometries and removed %d near-duplicates (RMSD threshold %.4f Ang).",
            len(final_list),
            len(removed_duplicates),
            rmsd_threshold,
        )
        detail_limit = 20
        for candidate_name, kept_name, aligned_rmsd in removed_duplicates[:detail_limit]:
            clustering.cluster_logger.info(
                "Similarity pruning removed %s as a near-duplicate of %s (RMSD %.6f < %.6f).",
                candidate_name,
                kept_name,
                aligned_rmsd,
                rmsd_threshold,
            )
        if len(removed_duplicates) > detail_limit:
            clustering.cluster_logger.info(
                "Similarity pruning omitted %d additional near-duplicate details.",
                len(removed_duplicates) - detail_limit,
            )
    from pyar.selection.reports import print_energy_table

    print_energy_table(final_list)
    return final_list


def _prefer_connected_structures(molecules, policy="prefer"):
    """Prefer geometries that remain connected under a covalent-radius graph."""
    if not molecules:
        return molecules

    try:
        from pyar.sampling import trial_generator as trial_generation
    except Exception as exc:
        from pyar.selection import clustering

        clustering.cluster_logger.debug(
            "Connectivity filter unavailable, keeping original candidate pool: %s",
            exc,
        )
        return molecules

    connected = []
    disconnected = []
    for molecule in molecules:
        try:
            if trial_generation.broken(molecule):
                disconnected.append(molecule)
            else:
                connected.append(molecule)
        except Exception as exc:
            from pyar.selection import clustering

            clustering.cluster_logger.debug(
                "Connectivity check failed for %s, keeping it in the connected pool: %s",
                getattr(molecule, 'name', 'unknown'),
                exc,
            )
            connected.append(molecule)

    if not connected:
        from pyar.selection import clustering

        if policy == "strict":
            clustering.cluster_logger.warning(
                "Connectivity policy is strict, but all candidates are disconnected under the covalent graph; returning an empty pool."
            )
            return []

        clustering.cluster_logger.info(
            "Connectivity preference requested, but all candidates are disconnected under the covalent graph; retaining the original pool."
        )
        return molecules

    if disconnected:
        from pyar.selection import clustering

        clustering.cluster_logger.info(
            "Connectivity preference retained %d connected candidates and discarded %d disconnected candidates.",
            len(connected),
            len(disconnected),
        )

    return connected
