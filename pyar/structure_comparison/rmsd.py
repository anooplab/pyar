"""Legacy permutation-aware Kabsch RMSD primitives."""

from __future__ import annotations

import itertools
import math
from collections import Counter

import numpy as np
from scipy.optimize import linear_sum_assignment


def kabsch_rotation(reference, mobile):
    """Return the proper rotation matrix aligning ``mobile`` onto ``reference``."""
    covariance = mobile.T @ reference
    left_vectors, _, right_vectors = np.linalg.svd(covariance)
    correction = np.eye(3)
    if np.linalg.det(left_vectors @ right_vectors) < 0:
        correction[-1, -1] = -1.0
    return left_vectors @ correction @ right_vectors


def kabsch_rmsd(reference, mobile):
    """Return RMSD after optimal translation and proper rotation."""
    reference_centered = reference - np.mean(reference, axis=0)
    mobile_centered = mobile - np.mean(mobile, axis=0)
    rotation = kabsch_rotation(reference_centered, mobile_centered)
    aligned = mobile_centered @ rotation
    diff = aligned - reference_centered
    return float(np.sqrt(np.mean(np.sum(diff * diff, axis=1))))


def equivalent_atom_groups(atoms):
    """Return indices grouped by element, in deterministic element order."""
    return {
        element: [index for index, atom in enumerate(atoms) if atom == element]
        for element in sorted(set(atoms))
    }


def exact_element_orders(reference_atoms, mobile_atoms, maximum_permutations=720):
    """Yield element-preserving mobile atom orders when exact matching is tractable."""
    reference_groups = equivalent_atom_groups(reference_atoms)
    mobile_groups = equivalent_atom_groups(mobile_atoms)
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
                    reference_groups[element], element_order,
                ):
                    mobile_order[reference_index] = mobile_index
            yield mobile_order

    return orders()


def assigned_element_order(reference_atoms, reference, mobile_atoms, mobile):
    """Match mobile atoms to reference atoms with element-constrained assignment."""
    order = [None] * len(reference_atoms)
    reference_groups = equivalent_atom_groups(reference_atoms)
    mobile_groups = equivalent_atom_groups(mobile_atoms)
    for element, reference_indices in reference_groups.items():
        mobile_indices = mobile_groups[element]
        distance_matrix = np.sum(
            (
                reference[np.asarray(reference_indices)][:, None, :]
                - mobile[np.asarray(mobile_indices)][None, :, :]
            ) ** 2,
            axis=2,
        )
        row_indices, column_indices = linear_sum_assignment(distance_matrix)
        for row_index, column_index in zip(row_indices, column_indices):
            order[reference_indices[row_index]] = mobile_indices[column_index]
    return order


def iterative_assigned_rmsd(reference_atoms, reference, mobile_atoms, mobile):
    """Estimate permutation-aware RMSD for systems too large for exact matching."""
    reference_centered = reference - np.mean(reference, axis=0)
    mobile_centered = mobile - np.mean(mobile, axis=0)
    order = assigned_element_order(
        reference_atoms, reference_centered, mobile_atoms, mobile_centered,
    )
    for _ in range(8):
        rotation = kabsch_rotation(reference_centered, mobile_centered[order])
        aligned_mobile = mobile_centered @ rotation
        next_order = assigned_element_order(
            reference_atoms, reference_centered, mobile_atoms, aligned_mobile,
        )
        if next_order == order:
            break
        order = next_order
    return kabsch_rmsd(reference, mobile[order])


def rmsd_after_alignment(candidate, kept):
    """Return translation-, rotation-, and element-order-invariant RMSD."""
    candidate_coords = np.asarray(candidate.coordinates, dtype=float)
    kept_coords = np.asarray(kept.coordinates, dtype=float)
    if (
        candidate_coords.shape != kept_coords.shape
        or Counter(candidate.atoms_list) != Counter(kept.atoms_list)
    ):
        return float("inf")

    candidate_orders = exact_element_orders(kept.atoms_list, candidate.atoms_list)
    if candidate_orders is not None:
        return min(
            kabsch_rmsd(kept_coords, candidate_coords[order])
            for order in candidate_orders
        )

    return iterative_assigned_rmsd(
        kept.atoms_list, kept_coords, candidate.atoms_list, candidate_coords,
    )


__all__ = [
    "assigned_element_order", "equivalent_atom_groups", "exact_element_orders",
    "iterative_assigned_rmsd", "kabsch_rotation", "kabsch_rmsd",
    "rmsd_after_alignment",
]
