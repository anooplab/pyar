"""Differentiable, opt-in outer-sphere confinement for shell searches."""

from __future__ import annotations

import numpy as np


def outer_sphere_energy_gradient(coordinates_angstrom, masses, solute_indices, *, radius_angstrom,
                                 force_constant_hartree_per_angstrom2, solvent_fragments=None):
    """Return wall energy (Eh) and Cartesian gradient (Eh/Angstrom).

    Solvent-fragment centers outside the sphere centered on the instantaneous
    solute COM receive a flat-bottom harmonic penalty. Fragment and solute COM
    derivatives are distributed over their atoms by mass fractions.
    """
    xyz = np.asarray(coordinates_angstrom, dtype=float)
    masses = np.asarray(masses, dtype=float)
    solute = np.asarray(solute_indices, dtype=int)
    if xyz.ndim != 2 or xyz.shape[1] != 3 or masses.shape != (len(xyz),):
        raise ValueError("coordinates and masses must have shapes (N, 3) and (N,)")
    if not np.isfinite(xyz).all() or not np.isfinite(masses).all() or np.any(masses <= 0):
        raise ValueError("wall coordinates and positive masses must be finite")
    if not len(solute) or len(set(solute.tolist())) != len(solute) or np.any(solute < 0) or np.any(solute >= len(xyz)):
        raise ValueError("solute atom indices must be distinct valid indices")
    radius = float(radius_angstrom)
    strength = float(force_constant_hartree_per_angstrom2)
    if not np.isfinite(radius) or radius <= 0 or not np.isfinite(strength) or strength < 0:
        raise ValueError("wall radius must be positive and strength finite and non-negative")
    solvent = np.setdiff1d(np.arange(len(xyz)), solute, assume_unique=False)
    gradient = np.zeros_like(xyz)
    if not len(solvent) or strength == 0:
        return 0.0, gradient
    fragments = [np.asarray(fragment, dtype=int) for fragment in (solvent_fragments or [solvent])]
    flattened = np.concatenate(fragments) if fragments else np.asarray([], dtype=int)
    if (any(not len(fragment) or np.any(fragment < 0) or np.any(fragment >= len(xyz)) for fragment in fragments)
            or len(np.unique(flattened)) != len(flattened)
            or not np.array_equal(np.sort(flattened), solvent)):
        raise ValueError("solvent fragments must partition all non-solute atoms")
    mass_fraction = masses[solute] / masses[solute].sum()
    center = np.average(xyz[solute], axis=0, weights=masses[solute])
    fragment_centers = np.asarray([
        np.average(xyz[fragment], axis=0, weights=masses[fragment]) for fragment in fragments
    ])
    displacement = fragment_centers - center
    distances = np.linalg.norm(displacement, axis=1)
    excess = np.maximum(distances - radius, 0.0)
    active = distances > radius
    fragment_gradient = np.zeros_like(displacement)
    fragment_gradient[active] = (strength * excess[active] / distances[active])[:, None] * displacement[active]
    for fragment, fragment_force in zip(fragments, fragment_gradient):
        gradient[fragment] += masses[fragment, None] / masses[fragment].sum() * fragment_force
    gradient[solute] -= mass_fraction[:, None] * fragment_gradient.sum(axis=0)
    energy = 0.5 * strength * float(np.dot(excess, excess))
    if not np.isfinite(energy) or not np.isfinite(gradient).all():
        raise ValueError("outer-sphere wall produced non-finite energy or gradient")
    return energy, gradient
