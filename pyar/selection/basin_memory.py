"""Persistence helpers for basin selection memory."""

from __future__ import annotations

import json
import os
from collections import Counter

import numpy as np

import pyar.representations


def _stoichiometry_label(molecule):
    """Return a compact stoichiometry label for basin snapshots."""
    counts = Counter(molecule.atoms_list)
    parts = []

    if 'C' in counts:
        carbon = counts.pop('C')
        parts.append('C' if carbon == 1 else f'C{carbon}')
    if 'H' in counts:
        hydrogen = counts.pop('H')
        parts.append('H' if hydrogen == 1 else f'H{hydrogen}')

    for element in sorted(counts):
        count = counts[element]
        parts.append(element if count == 1 else f'{element}{count}')

    return ''.join(parts) if parts else 'unknown'


def _basin_registry_path(molecule, output_root='selected', group_by_stoichiometry=True):
    """Return the persistence path for basin memory."""
    if group_by_stoichiometry:
        return os.path.join(
            output_root,
            f'stoichiometry_{_stoichiometry_label(molecule)}',
            'basin_registry.json',
        )
    return os.path.join(output_root, 'basin_registry.json')


def _fingerprint_signature(molecule):
    """Return a stable, normalized fingerprint signature for a molecule."""
    coordinates = np.asarray(molecule.coordinates, dtype=float)
    signature = pyar.representations.fingerprint(molecule.atoms_list, coordinates)
    signature = np.real_if_close(signature)
    signature = np.asarray(signature, dtype=float).ravel()
    norm = np.linalg.norm(signature)
    if norm > 0:
        signature = signature / norm
    return signature.tolist()


def _entry_fingerprint(entry):
    """Return a registry fingerprint as an array, or None if it is invalid."""
    try:
        return np.asarray(entry['fingerprint'], dtype=float).ravel()
    except (KeyError, TypeError, ValueError):
        return None


def _load_basin_registry(registry_path):
    """Load previously selected basin signatures from disk."""
    if not registry_path or not os.path.exists(registry_path):
        return []

    try:
        with open(registry_path, 'r') as fp:
            payload = json.load(fp)
    except Exception as exc:
        from pyar.selection import clustering

        clustering.cluster_logger.warning("Could not load basin registry %s: %s", registry_path, exc)
        return []

    entries = payload.get('entries', []) if isinstance(payload, dict) else payload
    return [entry for entry in entries if isinstance(entry, dict) and 'fingerprint' in entry]


def _persist_basin_registry(registry_path, selected_molecules, existing_entries=None, max_entries=200):
    """Persist selected basins so later runs can avoid rediscovering them."""
    if not registry_path or not selected_molecules:
        return

    if existing_entries is None:
        existing_entries = _load_basin_registry(registry_path)

    seen_signatures = {
        tuple(np.round(fingerprint, 6))
        for fingerprint in (_entry_fingerprint(entry) for entry in existing_entries)
        if fingerprint is not None
    }
    updated_entries = list(existing_entries)
    for molecule in selected_molecules:
        signature = _fingerprint_signature(molecule)
        signature_key = tuple(np.round(np.asarray(signature, dtype=float), 6))
        if signature_key in seen_signatures:
            continue
        seen_signatures.add(signature_key)
        updated_entries.append(
            {
                'name': molecule.name,
                'energy': float(molecule.energy),
                'fingerprint': signature,
            }
        )

    if not updated_entries:
        return

    target_directory = os.path.dirname(registry_path)
    if target_directory and not os.path.exists(target_directory):
        os.makedirs(target_directory, exist_ok=True)

    payload = {
        'stoichiometry': os.path.basename(os.path.dirname(registry_path)).replace('stoichiometry_', ''),
        'entries': updated_entries[-max_entries:],
    }
    with open(registry_path, 'w') as fp:
        json.dump(payload, fp, indent=2)

    from pyar.selection import clustering

    clustering.cluster_logger.info(
        "Updated basin registry with %d selected geometries at %s",
        len(selected_molecules),
        registry_path,
    )


def _basin_novelty_scores(molecules, basin_entries):
    """Rank molecules by novelty against previous basin signatures."""
    basin_signatures = [
        fingerprint
        for fingerprint in (_entry_fingerprint(entry) for entry in basin_entries)
        if fingerprint is not None
    ]
    scores = []
    for molecule in molecules:
        signature = np.asarray(_fingerprint_signature(molecule), dtype=float)
        comparable_signatures = [
            basin for basin in basin_signatures if basin.shape == signature.shape
        ]
        if comparable_signatures:
            novelty = min(np.linalg.norm(signature - basin) for basin in comparable_signatures)
        else:
            novelty = np.inf
        scores.append((novelty, float(molecule.energy), molecule))
    return scores


def _apply_basin_memory(molecules, maximum_number_of_seeds, basin_entries):
    """Reduce the MBTR candidate pool using previous basin discoveries."""
    if not basin_entries or len(molecules) <= maximum_number_of_seeds:
        return molecules

    pool_cap = min(len(molecules), maximum_number_of_seeds * 2)
    if len(molecules) <= pool_cap:
        return molecules

    scored_molecules = _basin_novelty_scores(molecules, basin_entries)
    novelty_ranked = sorted(scored_molecules, key=lambda item: (-item[0], item[1], item[2].name))
    energy_ranked = sorted(scored_molecules, key=lambda item: (item[1], item[2].name))

    selected = []
    selected_ids = set()
    for _, _, molecule in energy_ranked[:maximum_number_of_seeds]:
        selected.append(molecule)
        selected_ids.add(id(molecule))
    for _, _, molecule in novelty_ranked:
        if len(selected) >= pool_cap:
            break
        if id(molecule) in selected_ids:
            continue
        selected.append(molecule)
        selected_ids.add(id(molecule))

    from pyar.selection import clustering

    clustering.cluster_logger.info(
        "Basin memory reduced candidate pool from %d to %d before MBTR selection.",
        len(molecules),
        len(selected),
    )
    return selected


def record_selected_basins(selected_molecules, output_root='selected'):
    """Persist final selected geometries in the per-stoichiometry basin registry."""
    if not selected_molecules:
        return None
    registry_path = _basin_registry_path(selected_molecules[0], output_root=output_root)
    _persist_basin_registry(registry_path, selected_molecules)
    return registry_path
