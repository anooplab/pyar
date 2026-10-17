"""Versioned, bounded rigid-body moves for fixed-composition clusters."""

from __future__ import annotations

import numpy as np


def _rotation(rng):
    axis = rng.normal(size=3)
    axis /= np.linalg.norm(axis)
    angle = rng.uniform(-np.pi, np.pi)
    cross = np.array([[0., -axis[2], axis[1]], [axis[2], 0., -axis[0]], [-axis[1], axis[0], 0.]])
    return np.eye(3) + np.sin(angle) * cross + (1. - np.cos(angle)) * (cross @ cross)


def propose_rigid_move(molecule, rng):
    """Return a copied cluster after translating/rotating one or several waters."""
    candidate = molecule.copy()
    candidate.solute_atom_indices = tuple(molecule.solute_atom_indices)
    candidate.solvent_fragments = tuple(tuple(part) for part in molecule.solvent_fragments)
    candidate.fragments = [list(candidate.solute_atom_indices), *[list(part) for part in candidate.solvent_fragments]]
    for attribute in ("microsolvation_coverage", "microsolvation_target_atom", "microsolvation_target_kind",
                      "microsolvation_target_fragment"):
        if hasattr(molecule, attribute):
            setattr(candidate, attribute, getattr(molecule, attribute))
    fragments = [tuple(map(int, part)) for part in molecule.solvent_fragments]
    if not fragments:
        raise ValueError("whole-cluster perturbation requires at least one solvent fragment")
    if len(fragments) >= 2 and rng.random() < 0.15:
        left, right = sorted(map(int, rng.choice(len(fragments), size=2, replace=False)))
        a, b = np.asarray(fragments[left], dtype=int), np.asarray(fragments[right], dtype=int)
        if len(a) != len(b):
            raise ValueError("solvent-site exchange requires equal fragment atom counts")
        coords = np.asarray(candidate.coordinates, dtype=float).copy()
        masses = np.asarray(candidate.atomic_mass, dtype=float)
        center_a = np.average(coords[a], axis=0, weights=masses[a])
        center_b = np.average(coords[b], axis=0, weights=masses[b])
        shape_a, shape_b = coords[a] - center_a, coords[b] - center_b
        coords[a], coords[b] = shape_b + center_a, shape_a + center_b
        candidate.coordinates = coords
        candidate.fragments_coordinates = candidate.split_coordinates()
        candidate.fragments_atoms_list = candidate.split_atoms_lists()
        return candidate, {"policy": "whole-cluster-perturbation-v1", "move": "solvent-site-exchange",
                           "fragments": [left, right], "mode": "exchange-centered-rigid-fragment-geometries"}
    nmove = 1 if len(fragments) == 1 or rng.random() < 0.75 else min(3, len(fragments))
    chosen = sorted(map(int, rng.choice(len(fragments), size=nmove, replace=False)))
    coords = np.asarray(candidate.coordinates, dtype=float).copy()
    records = []
    for index in chosen:
        ids = np.asarray(fragments[index], dtype=int)
        center = coords[ids].mean(axis=0)
        rotation = _rotation(rng)
        distance = rng.uniform(0.5, 1.5)
        direction = rng.normal(size=3); direction /= np.linalg.norm(direction)
        coords[ids] = (coords[ids] - center) @ rotation.T + center + distance * direction
        records.append({"fragment_index": index, "rotation": rotation.tolist(),
                        "translation_angstrom": (distance * direction).tolist()})
    candidate.coordinates = coords
    candidate.fragments_coordinates = candidate.split_coordinates()
    candidate.fragments_atoms_list = candidate.split_atoms_lists()
    return candidate, {"policy": "whole-cluster-perturbation-v1", "move": "single-fragment-rigid" if nmove == 1 else "multi-fragment-rigid",
                       "fragments": records}


def has_interfragment_clash(molecule, scale=0.85):
    """Reject overlaps between distinct fragments using vdW radii (placement-compatible scale)."""
    from pyar.data import new_atomic_data as atomic_data
    groups = [tuple(map(int, molecule.solute_atom_indices)), *[tuple(map(int, x)) for x in molecule.solvent_fragments]]
    xyz = np.asarray(molecule.coordinates, dtype=float)
    symbols = molecule.atoms_list
    radii = np.asarray([atomic_data.vdw_radius[atom] for atom in symbols], dtype=float)
    for left in range(len(groups)):
        for right in range(left + 1, len(groups)):
            a, b = np.asarray(groups[left], int), np.asarray(groups[right], int)
            if len(a) and len(b) and np.any(np.linalg.norm(xyz[a, None] - xyz[None, b], axis=2) < scale * (radii[a, None] + radii[None, b])):
                return True
    return False


def preserves_cluster_identity(molecule, request, expected_multiplicity):
    """Check atom order, fragment partition, electronic state and solute graph."""
    solute_count = len(request.solute)
    solvent_count = len(request.solvent)
    expected_atoms = list(request.solute.atoms_list) + list(request.solvent.atoms_list) * request.count
    expected_fragments = [tuple(range(solute_count + i * solvent_count,
                                      solute_count + (i + 1) * solvent_count))
                          for i in range(request.count)]
    if (list(molecule.atoms_list) != expected_atoms
            or tuple(molecule.solute_atom_indices) != tuple(range(solute_count))
            or tuple(tuple(part) for part in molecule.solvent_fragments) != tuple(expected_fragments)
            or int(molecule.charge) != int(request.solute.charge) + request.count * int(request.solvent.charge)
            or int(molecule.multiplicity) != int(expected_multiplicity)):
        return False
    from pyar.property import get_bond_matrix
    original = get_bond_matrix(np.asarray(request.solute.coordinates), np.asarray(request.solute.covalent_radius))
    observed = get_bond_matrix(np.asarray(molecule.coordinates[:solute_count]),
                               np.asarray(request.solute.covalent_radius))
    return np.array_equal(original, observed)
