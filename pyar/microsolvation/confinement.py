"""Post-optimization first-shell validity checks.

This module intentionally does not claim to apply an optimization restraint.
It detects solvent migration outside the requested surface shell after
geometry optimization and provides the same geometric validation in
geometry-only mode.
"""

from __future__ import annotations

import numpy as np


def shell_postfilter(molecule, target_surface, *, tolerance):
    """Evaluate whether every solvent-fragment centre remains near the target surface."""
    solute_count = len(molecule.solute_atom_indices)
    solvent_centers = [
        np.mean(molecule.coordinates[np.asarray(fragment, dtype=int)], axis=0)
        for fragment in molecule.solvent_fragments
    ]
    if len(target_surface.points):
        distances = [
            float(np.min(np.linalg.norm(target_surface.points - center, axis=1)))
            for center in solvent_centers
        ]
    else:
        distances = [float("inf") for _ in solvent_centers]
    maximum = max(distances, default=0.0)
    return {
        "shell_valid": bool(distances) and maximum <= float(tolerance),
        "max_solvent_center_to_surface_angstrom": maximum,
        "shell_tolerance_angstrom": float(tolerance),
        "confinement_mode": "post-optimization-shell-validity-filter",
        "optimization_restraint_applied": False,
        "solute_atom_count": solute_count,
    }
