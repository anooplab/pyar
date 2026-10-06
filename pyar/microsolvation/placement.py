"""Surface-targeted solvent placement with full-cluster clash checking."""

from __future__ import annotations

import numpy as np

from pyar.core.molecule import Molecule
from pyar.sampling.rotation import halton_quaternions, quaternions_to_euler_zxz
from pyar.microsolvation.surface import build_solute_surface, surface_coverage
from pyar.microsolvation.confinement import shell_postfilter


def _surface_for_seed(seed, solute_atom_count, request):
    solute = Molecule(
        seed.atoms_list[:solute_atom_count], seed.coordinates[:solute_atom_count],
        name=request.solute.name, title=request.solute.title,
    )
    return build_solute_surface(
        solute, points_per_atom=request.surface_points_per_atom,
        probe_radius=request.probe_radius,
    ).restrict(request.site)


def _existing_solvent_coordinates(seed):
    start = len(seed.solute_atom_indices)
    return seed.coordinates[start:]


def _target_indices(surface, seed, solvent_count):
    """Order eligible original-solute points by deterministic max-min coverage."""
    if not len(surface.points):
        raise ValueError("No accessible surface points remain for the selected solute region")
    solvent_coordinates = _existing_solvent_coordinates(seed)
    coverage = surface_coverage(
        surface, seed.atoms_list[len(seed.solute_atom_indices):], solvent_coordinates,
        probe_radius=0.0,
    ) if len(solvent_coordinates) else {
        "covered_mask": np.zeros(len(surface.points), dtype=bool),
        "coverage_fraction": 0.0,
    }
    open_indices = np.flatnonzero(~coverage["covered_mask"])
    saturated = len(open_indices) == 0
    candidates = open_indices if len(open_indices) else np.arange(len(surface.points))
    centers = []
    for fragment in getattr(seed, "solvent_fragments", ()):
        centers.append(np.mean(seed.coordinates[np.asarray(fragment, dtype=int)], axis=0))
    chosen = []
    available = candidates.copy()
    while len(chosen) < min(solvent_count, len(candidates)):
        points = surface.points[available]
        if centers or chosen:
            references = centers + [surface.points[index] for index in chosen]
            distances = np.linalg.norm(points[:, None, :] - np.asarray(references)[None, :, :], axis=2)
            scores = distances.min(axis=1)
        elif saturated:
            # Once occupied, prefer the least-covered points, deterministically.
            center = np.mean(surface.points, axis=0)
            scores = np.linalg.norm(points - center, axis=1)
        else:
            center = np.mean(surface.points, axis=0)
            scores = np.linalg.norm(points - center, axis=1)
        selected_position = int(np.argmax(scores))
        chosen.append(int(available[selected_position]))
        available = np.delete(available, selected_position)
    if not chosen:
        chosen = [int(candidates[0])]
    return chosen, saturated


def _has_clash(cluster, incoming, *, scale=0.85):
    distances = np.linalg.norm(
        cluster.coordinates[:, None, :] - incoming.coordinates[None, :, :], axis=2
    )
    cluster_radii = np.asarray(cluster.vdw_radius, dtype=float)
    incoming_radii = np.asarray(incoming.vdw_radius, dtype=float)
    cluster_radii = np.where(np.isfinite(cluster_radii), cluster_radii, cluster.covalent_radius)
    incoming_radii = np.where(np.isfinite(incoming_radii), incoming_radii, incoming.covalent_radius)
    radii = cluster_radii[:, None] + incoming_radii[None, :]
    return bool(np.any(distances < scale * radii))


def _place_at_target(solvent, target, normal, collision_environment, *, maximum_push=5.0):
    """Place a solvent COM at a target point, pushing outward to clear all atoms."""
    candidate = solvent.copy()
    candidate.move_to_origin()
    candidate.translate(target - candidate.centroid)
    displacement = 0.0
    while _has_clash(collision_environment, candidate):
        displacement += 0.2
        if displacement > maximum_push:
            return None
        candidate.translate(0.2 * normal)
    return candidate


def _combine(seed, solvent, request, name):
    combined = seed.merged_with(solvent)
    combined.name = name
    combined.title = f"{request.solute.name} + {len(seed.solvent_fragments) + 1} solvent molecule(s)"
    solute_indices = tuple(seed.solute_atom_indices)
    previous_fragments = tuple(tuple(int(i) for i in part) for part in getattr(seed, "solvent_fragments", ()))
    start = len(seed.atoms_list)
    combined.solute_atom_indices = solute_indices
    combined.solvent_fragments = previous_fragments + (tuple(range(start, start + len(solvent))),)
    combined.fragments = [list(solute_indices), *[list(x) for x in combined.solvent_fragments]]
    combined.fragments_coordinates = combined.split_coordinates()
    combined.fragments_atoms_list = combined.split_atoms_lists()
    from pyar.molecule_merge import combine_multiplicity
    combined.charge = seed.charge + request.solvent.charge
    combined.multiplicity = combine_multiplicity(seed.multiplicity, request.solvent.multiplicity)
    combined.scftype = "rhf" if seed.scftype == request.solvent.scftype == "rhf" else "uhf"
    return combined


def generate_microsolvation_candidates(seed, request, *, step, solute_atom_count):
    """Generate a bounded set targeting only original-solute surface samples.

    Previously accepted solvent molecules are used by ``_place_at_target`` as
    collision geometry and by the coverage model. They never generate target
    points; all targets come from a surface built from the fixed solute atom
    prefix of the current seed.
    """
    if not hasattr(seed, "solute_atom_indices"):
        seed.solute_atom_indices = tuple(range(solute_atom_count))
        seed.solvent_fragments = ()
    surface = _surface_for_seed(seed, solute_atom_count, request)
    target_ids, saturated = _target_indices(surface, seed, request.number_of_orientations)
    quaternions = halton_quaternions(request.number_of_orientations, seed=step)
    angles = quaternions_to_euler_zxz(quaternions)
    candidates = []
    for orientation_index in range(request.number_of_orientations):
        target_index = target_ids[orientation_index % len(target_ids)]
        solvent = request.solvent.copy()
        solvent.move_to_origin()
        if len(solvent) > 1:
            solvent.rotate_3d(angles[orientation_index])
        placed = _place_at_target(
            solvent, surface.points[target_index], surface.normals[target_index], seed,
        )
        if placed is None:
            continue
        molecule = _combine(seed, placed, request, f"micro_{step:03d}_{len(candidates):03d}")
        molecule.microsolvation_target_atom = int(surface.parent_atoms[target_index])
        molecule.microsolvation_target_point = tuple(float(x) for x in surface.points[target_index])
        molecule.surface_saturated_at_placement = saturated
        molecule.microsolvation_surface = surface
        candidates.append(molecule)
    return candidates


def assess_shell(molecule, request):
    """Return coverage and first-shell validity using the original solute prefix."""
    solute_count = len(molecule.solute_atom_indices)
    solute = Molecule(
        molecule.atoms_list[:solute_count], molecule.coordinates[:solute_count],
        name=request.solute.name, title=request.solute.title,
    )
    full_surface = build_solute_surface(
        solute, points_per_atom=request.surface_points_per_atom,
        probe_radius=request.probe_radius,
    )
    target_surface = full_surface.restrict(request.site) if request.site is not None else full_surface
    solvent_coordinates = molecule.coordinates[solute_count:]
    solvent_atoms = molecule.atoms_list[solute_count:]
    coverage = surface_coverage(target_surface, solvent_atoms, solvent_coordinates)
    full_coverage = surface_coverage(full_surface, solvent_atoms, solvent_coordinates)
    return {
        key: value for key, value in coverage.items() if key != "covered_mask"
    } | {
        "whole_solute_coverage_fraction": full_coverage["coverage_fraction"],
        "whole_solute_surface_points": full_coverage["surface_points"],
        "whole_solute_covered_points": full_coverage["covered_points"],
        **shell_postfilter(molecule, target_surface, tolerance=request.shell_tolerance),
        "site_atoms": None if request.site is None else list(request.site),
    }
