"""Deterministic accessible surface points and solvent coverage metrics.

The surface uses atom-centred van der Waals spheres expanded by a spherical
probe radius. It is a discretized placement surface, not an analytical SASA.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from pyar.data import new_atomic_data as atomic_data
from pyar.sampling.sphere import fibonacci_sphere


@dataclass(frozen=True)
class SoluteSurface:
    points: np.ndarray
    normals: np.ndarray
    parent_atoms: np.ndarray
    weights: np.ndarray

    def restrict(self, atom_indices):
        if atom_indices is None:
            return self
        mask = np.isin(self.parent_atoms, np.asarray(atom_indices, dtype=int))
        return SoluteSurface(self.points[mask], self.normals[mask], self.parent_atoms[mask], self.weights[mask])


def _vdw_radius(symbol):
    value = atomic_data.vdw_radius.get(symbol)
    if value is None or not np.isfinite(value):
        value = atomic_data.covalent_radius[symbol]
    return float(value)


def build_solute_surface(solute, *, points_per_atom=96, probe_radius=1.4):
    """Return exposed Fibonacci samples on probe-expanded solute atom spheres."""
    points_per_atom = int(points_per_atom)
    if points_per_atom < 1 or not np.isfinite(probe_radius) or probe_radius < 0:
        raise ValueError("surface sampling requires positive points_per_atom and non-negative probe_radius")
    directions = fibonacci_sphere(points_per_atom)
    coordinates = np.asarray(solute.coordinates, dtype=float)
    radii = np.asarray([_vdw_radius(atom) + float(probe_radius) for atom in solute.atoms_list])
    positions, normals, parents, weights = [], [], [], []
    for atom_index, (center, radius) in enumerate(zip(coordinates, radii)):
        candidates = center + radius * directions
        exposed = np.ones(points_per_atom, dtype=bool)
        for neighbor_index, (neighbor, neighbor_radius) in enumerate(zip(coordinates, radii)):
            if neighbor_index == atom_index:
                continue
            exposed &= np.linalg.norm(candidates - neighbor, axis=1) >= neighbor_radius
        for direction, point in zip(directions[exposed], candidates[exposed]):
            positions.append(point)
            normals.append(direction)
            parents.append(atom_index)
            weights.append(4.0 * np.pi * radius * radius / points_per_atom)
    if not positions:
        raise ValueError("No accessible solute surface points were generated")
    return SoluteSurface(np.asarray(positions), np.asarray(normals), np.asarray(parents, dtype=int), np.asarray(weights))


def surface_coverage(surface, solvent_atoms, solvent_coordinates, *, probe_radius=0.0):
    """Return covered/open surface counts and area-weighted coverage fraction."""
    atom_radii = np.asarray([_vdw_radius(atom) + float(probe_radius) for atom in solvent_atoms])
    coordinates = np.asarray(solvent_coordinates, dtype=float)
    distances = np.linalg.norm(surface.points[:, None, :] - coordinates[None, :, :], axis=2)
    covered = np.any(distances <= atom_radii[None, :], axis=1)
    total_weight = float(np.sum(surface.weights))
    covered_weight = float(np.sum(surface.weights[covered]))
    return {
        "surface_points": int(len(surface.points)),
        "covered_points": int(np.count_nonzero(covered)),
        "open_points": int(len(surface.points) - np.count_nonzero(covered)),
        "coverage_fraction": 0.0 if total_weight == 0 else covered_weight / total_weight,
        "covered_mask": covered,
    }
