"""Deduplication and connectivity filtering for selected geometries."""

from __future__ import annotations

import itertools
import math
from collections import Counter

import numpy as np
from scipy.optimize import linear_sum_assignment

import pyar.representations


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
            if not _structure_is_similar(left, right):
                continue
            try:
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


def _kabsch_rotation(reference, mobile):
    """Return the proper rotation matrix aligning ``mobile`` onto ``reference``."""
    covariance = mobile.T @ reference
    left_vectors, _, right_vectors = np.linalg.svd(covariance)
    correction = np.eye(3)
    if np.linalg.det(left_vectors @ right_vectors) < 0:
        correction[-1, -1] = -1.0
    return left_vectors @ correction @ right_vectors


def _kabsch_rmsd(reference, mobile):
    """Return RMSD after optimal translation and proper rotation."""
    reference_centered = reference - np.mean(reference, axis=0)
    mobile_centered = mobile - np.mean(mobile, axis=0)
    rotation = _kabsch_rotation(reference_centered, mobile_centered)
    aligned = mobile_centered @ rotation
    diff = aligned - reference_centered
    return float(np.sqrt(np.mean(np.sum(diff * diff, axis=1))))


def _equivalent_atom_groups(atoms):
    """Return indices grouped by element, in deterministic element order."""
    return {
        element: [index for index, atom in enumerate(atoms) if atom == element]
        for element in sorted(set(atoms))
    }


def _exact_element_orders(reference_atoms, mobile_atoms, maximum_permutations=720):
    """Yield element-preserving mobile atom orders when exact matching is tractable."""
    reference_groups = _equivalent_atom_groups(reference_atoms)
    mobile_groups = _equivalent_atom_groups(mobile_atoms)
    permutation_count = math.prod(
        math.factorial(len(indices)) for indices in mobile_groups.values()
    )
    if permutation_count > maximum_permutations:
        return None

    elements = list(reference_groups)
    permutations = [
        list(itertools.permutations(mobile_groups[element])) for element in elements
    ]

    def orders():
        for element_orders in itertools.product(*permutations):
            mobile_order = [None] * len(reference_atoms)
            for element, element_order in zip(elements, element_orders):
                for reference_index, mobile_index in zip(
                    reference_groups[element],
                    element_order,
                ):
                    mobile_order[reference_index] = mobile_index
            yield mobile_order

    return orders()


def _assigned_element_order(reference_atoms, reference, mobile_atoms, mobile):
    """Match mobile atoms to reference atoms with element-constrained assignment."""
    order = [None] * len(reference_atoms)
    reference_groups = _equivalent_atom_groups(reference_atoms)
    mobile_groups = _equivalent_atom_groups(mobile_atoms)
    for element, reference_indices in reference_groups.items():
        mobile_indices = mobile_groups[element]
        distance_matrix = np.sum(
            (
                reference[np.asarray(reference_indices)][:, None, :]
                - mobile[np.asarray(mobile_indices)][None, :, :]
            )
            ** 2,
            axis=2,
        )
        row_indices, column_indices = linear_sum_assignment(distance_matrix)
        for row_index, column_index in zip(row_indices, column_indices):
            order[reference_indices[row_index]] = mobile_indices[column_index]
    return order


def _iterative_assigned_rmsd(reference_atoms, reference, mobile_atoms, mobile):
    """Estimate permutation-aware RMSD for systems too large for exact matching."""
    reference_centered = reference - np.mean(reference, axis=0)
    mobile_centered = mobile - np.mean(mobile, axis=0)
    order = _assigned_element_order(
        reference_atoms,
        reference_centered,
        mobile_atoms,
        mobile_centered,
    )
    for _ in range(8):
        rotation = _kabsch_rotation(reference_centered, mobile_centered[order])
        aligned_mobile = mobile_centered @ rotation
        next_order = _assigned_element_order(
            reference_atoms,
            reference_centered,
            mobile_atoms,
            aligned_mobile,
        )
        if next_order == order:
            break
        order = next_order
    return _kabsch_rmsd(reference, mobile[order])


def _rmsd_after_alignment(candidate, kept):
    """Return translation-, rotation-, and element-order-invariant RMSD."""
    candidate_coords = np.asarray(candidate.coordinates, dtype=float)
    kept_coords = np.asarray(kept.coordinates, dtype=float)
    if (
        candidate_coords.shape != kept_coords.shape
        or Counter(candidate.atoms_list) != Counter(kept.atoms_list)
    ):
        return float('inf')

    candidate_orders = _exact_element_orders(kept.atoms_list, candidate.atoms_list)
    if candidate_orders is not None:
        return min(
            _kabsch_rmsd(kept_coords, candidate_coords[order])
            for order in candidate_orders
        )

    return _iterative_assigned_rmsd(
        kept.atoms_list,
        kept_coords,
        candidate.atoms_list,
        candidate_coords,
    )


def remove_similar(list_of_molecules):
    """Remove near-duplicate geometries from a candidate pool."""
    from pyar.selection import clustering

    ordered_molecules = sorted(list_of_molecules, key=lambda molecule: (float(molecule.energy), molecule.name))
    final_list = []
    rmsd_threshold = _adaptive_duplicate_rmsd_threshold(ordered_molecules)
    clustering.cluster_logger.debug('Number of molecules before similarity elimination,  {}'.format(len(ordered_molecules)))
    for candidate in ordered_molecules:
        duplicate = False
        for kept in final_list:
            if len(candidate.atoms_list) < 2 or len(kept.atoms_list) < 2:
                continue
            if not _structure_is_similar(candidate, kept):
                continue
            if _rmsd_after_alignment(candidate, kept) < rmsd_threshold:
                duplicate = True
                clustering.cluster_logger.debug(
                    'Removing {} as a near-duplicate of {}'.format(candidate.name, kept.name)
                )
                break
        if not duplicate:
            final_list.append(candidate)
    clustering.cluster_logger.debug('Number of molecules after similarity elimination,  {}'.format(len(final_list)))
    from pyar.selection.clustering import print_energy_table

    print_energy_table(final_list)
    return final_list


def _prefer_connected_structures(molecules):
    """Prefer geometries that remain connected under a covalent-radius graph."""
    if len(molecules) < 2:
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

        clustering.cluster_logger.warning(
            "All candidate geometries are disconnected; continuing with the original pool."
        )
        return molecules

    if disconnected:
        from pyar.selection import clustering

        clustering.cluster_logger.info(
            "Discarded %d disconnected candidate geometries before clustering.",
            len(disconnected),
        )

    return connected
