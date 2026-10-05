"""Energy-ranked filtering followed by optional shared graph-first pruning."""

import numpy as np

from pyar.selection.deduplication import deduplicate_structures
from pyar.selection.reports import HARTREE_TO_KCAL_MOL


def select_structures(molecules, *, within=None, top=None, unique=False, threshold=None):
    if within is None and top is None:
        raise ValueError('Specify --within or --top; --unique is a modifier')
    if within is not None and (not np.isfinite(within) or within < 0):
        raise ValueError('--within must be finite and nonnegative')
    if top is not None and (not isinstance(top, int) or top < 1):
        raise ValueError('--top must be a positive integer')
    molecules = list(molecules)
    for molecule in molecules:
        if molecule.energy is None or not np.isfinite(float(molecule.energy)):
            raise ValueError(f'Missing finite energy for {molecule.name}')
    ranked = sorted(molecules, key=lambda molecule: float(molecule.energy))
    if not ranked:
        return {'kept': [], 'minimum_energy_hartree': None, 'deduplication': None}
    minimum = float(ranked[0].energy)
    if within is not None:
        # Compare against the Hartree boundary directly to avoid cancellation
        # turning an exact kcal/mol boundary into a slightly larger difference.
        cutoff = minimum + within / HARTREE_TO_KCAL_MOL
        ranked = [molecule for molecule in ranked if float(molecule.energy) <= cutoff]
    if top is not None:
        ranked = ranked[:top]
    pruned = deduplicate_structures(ranked, threshold=threshold, ordering='input') if unique else None
    return {'kept': ranked if pruned is None else pruned['kept'],
            'minimum_energy_hartree': minimum, 'deduplication': pruned}
