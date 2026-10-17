"""Geometry-derived and case-triggered microsolvation wall policy."""

from __future__ import annotations

import numpy as np


WALL_POLICY_VERSION = "case-dependent-outer-sphere-v1"
POLICY_TEMPERATURE_K = 298.15
POLICY_PENALTY_KBT = 3.0
POLICY_OVERRUN_ANGSTROM = 0.5
POLICY_RADIAL_MARGIN_ANGSTROM = 0.5
_KB_HARTREE_PER_K = 3.166811563e-6


def derive_outer_sphere_parameters(request):
    """Derive a first-shell envelope and a soft, thermal-scale wall strength.

    The radius encloses the farthest sampled accessible solute surface point,
    the solvent fragment COM-to-atom extent, and a 0.5 A geometric allowance.
    The harmonic strength makes a 0.5 A overrun cost 3 kBT at 298.15 K. The
    thermal reference calibrates search stiffness only; it is not a solvent
    thermodynamics model or a sampling temperature.
    """
    from ase import Atoms
    from pyar.microsolvation.surface import build_solute_surface

    solute = request.solute
    solvent = request.solvent
    solute_atoms = Atoms(symbols=solute.atoms_list, positions=solute.coordinates)
    solvent_atoms = Atoms(symbols=solvent.atoms_list, positions=solvent.coordinates)
    solute_center = np.average(solute.coordinates, axis=0, weights=solute_atoms.get_masses())
    solvent_center = np.average(solvent.coordinates, axis=0, weights=solvent_atoms.get_masses())
    solvent_extent = float(np.max(np.linalg.norm(solvent.coordinates - solvent_center, axis=1)))
    surface = build_solute_surface(
        solute, points_per_atom=request.surface_points_per_atom,
        probe_radius=request.probe_radius,
    ).restrict(request.site)
    if not len(surface.points):
        raise ValueError("automatic wall policy needs accessible points on the requested solute surface")
    surface_radius = float(np.max(np.linalg.norm(surface.points - solute_center, axis=1)))
    radius = surface_radius + solvent_extent + POLICY_RADIAL_MARGIN_ANGSTROM
    target_overrun_energy = POLICY_PENALTY_KBT * _KB_HARTREE_PER_K * POLICY_TEMPERATURE_K
    force_constant = 2.0 * target_overrun_energy / POLICY_OVERRUN_ANGSTROM**2
    return {
        "policy": WALL_POLICY_VERSION,
        "radius_angstrom": radius,
        "force_constant_hartree_per_angstrom2": force_constant,
        "radius_derivation": {
            "model": "maximum-accessible-surface-radius-plus-solvent-fragment-com-extent-plus-margin",
            "solute_surface_radius_angstrom": surface_radius,
            "solvent_fragment_com_extent_angstrom": solvent_extent,
            "radial_margin_angstrom": POLICY_RADIAL_MARGIN_ANGSTROM,
            "solute_center": "mass-weighted-center-of-mass",
        },
        "strength_derivation": {
            "model": "three-kbt-at-half-angstrom-overrun-v1",
            "reference_temperature_kelvin": POLICY_TEMPERATURE_K,
            "penalty_kbt": POLICY_PENALTY_KBT,
            "overrun_angstrom": POLICY_OVERRUN_ANGSTROM,
            "target_penalty_hartree": target_overrun_energy,
        },
    }


def should_apply_automatic_wall(request, shell_assessment):
    """Activate confinement only when an unbiased shell optimization escaped."""
    return bool(
        request.wall == "auto"
        and request.objective == "shell"
        and not shell_assessment.get("shell_valid", False)
    )
