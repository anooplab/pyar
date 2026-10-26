"""Validated request model for solute-centred microsolvation."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping

import numpy as np

from pyar.aggregation.request import _molecule_signature


class MicrosolvationRequestError(ValueError):
    """Raised when a microsolvation request is scientifically invalid."""


def _validate_wall_mode_options(wall, wall_radius, wall_force_constant):
    """Validate wall mode and which modes accept explicit wall parameters."""
    if wall not in {"none", "outer-sphere", "auto"}:
        raise MicrosolvationRequestError("wall must be 'none', 'outer-sphere', or 'auto'")
    if wall == "none" and (wall_radius is not None or wall_force_constant is not None):
        raise MicrosolvationRequestError("--wall-radius and --wall-force-constant require --wall outer-sphere")
    if wall == "auto" and (wall_radius is not None or wall_force_constant is not None):
        raise MicrosolvationRequestError("automatic wall policy derives its own radius and strength")


@dataclass(frozen=True)
class MicrosolvationRequest:
    """Resolved settings for explicit shell or finite-cluster sampling."""

    solute: Any
    solvent: Any
    count: int
    number_of_orientations: int = 8
    backend_parameters: Mapping[str, Any] = field(default_factory=dict)
    maximum_number_of_seeds: int = 12
    site: tuple[int, ...] | None = None
    surface_points_per_atom: int = 96
    probe_radius: float = 1.4
    shell_tolerance: float = 3.5
    objective: str = "shell"
    placement_score: str = "geometric"
    site_groups: tuple[tuple[str, tuple[int, ...]], ...] = ()
    contact_distance: float | None = None
    wall: str = "none"
    wall_radius: float | None = None
    wall_force_constant: float | None = None
    resampling_moves: int = 0
    resampling_evaluation_budget: int = 0
    resampling_time_limit: float = 0.0
    resampling_seed: int = 0

    def __post_init__(self):
        for label, value in (
            ("count", self.count),
            ("orientations", self.number_of_orientations),
            ("maximum_number_of_seeds", self.maximum_number_of_seeds),
            ("surface_points_per_atom", self.surface_points_per_atom),
        ):
            if isinstance(value, bool) or not isinstance(value, int) or value < 1:
                raise MicrosolvationRequestError(f"{label} must be a positive integer")
        for label, molecule in (("solute", self.solute), ("solvent", self.solvent)):
            coordinates = np.asarray(getattr(molecule, "coordinates", None), dtype=float)
            atoms = getattr(molecule, "atoms_list", ())
            if not atoms or coordinates.shape != (len(atoms), 3) or not np.isfinite(coordinates).all():
                raise MicrosolvationRequestError(f"{label} must contain atoms with finite coordinates")
            charge = getattr(molecule, "charge", None)
            multiplicity = getattr(molecule, "multiplicity", None)
            if charge is None or multiplicity is None:
                raise MicrosolvationRequestError(f"{label} charge and multiplicity must be resolved")
            electrons = sum(molecule.atomic_number) - int(charge)
            if (electrons < 0 or int(multiplicity) < 1
                    or int(multiplicity) - 1 > electrons
                    or (electrons - (int(multiplicity) - 1)) % 2):
                raise MicrosolvationRequestError(f"Invalid charge/multiplicity for {label} {molecule.name}")
        if not np.isfinite(float(self.probe_radius)) or self.probe_radius < 0:
            raise MicrosolvationRequestError("probe radius must be finite and non-negative")
        if not np.isfinite(float(self.shell_tolerance)) or self.shell_tolerance <= 0:
            raise MicrosolvationRequestError("shell tolerance must be finite and positive")
        object.__setattr__(self, "backend_parameters", dict(self.backend_parameters or {}))
        if self.objective not in {"shell", "cluster"}:
            raise MicrosolvationRequestError("objective must be 'shell' or 'cluster'")
        _validate_wall_mode_options(self.wall, self.wall_radius, self.wall_force_constant)
        if self.wall == "outer-sphere":
            if self.objective != "shell":
                raise MicrosolvationRequestError("the outer-sphere wall is only supported for shell objective")
            if (isinstance(self.wall_radius, bool) or self.wall_radius is None
                    or not np.isfinite(float(self.wall_radius)) or float(self.wall_radius) <= 0):
                raise MicrosolvationRequestError("outer-sphere wall requires a positive --wall-radius in Angstrom")
            if (isinstance(self.wall_force_constant, bool) or self.wall_force_constant is None
                    or not np.isfinite(float(self.wall_force_constant))
                    or float(self.wall_force_constant) <= 0):
                raise MicrosolvationRequestError("outer-sphere wall requires a positive --wall-force-constant in hartree/angstrom^2")
            from pyar.backend_capabilities import backend_supports_geometry_optimization, normalize_backend_name
            software = self.backend_parameters.get("software")
            if not software:
                raise MicrosolvationRequestError("outer-sphere wall requires an energy-gradient backend")
            backend = normalize_backend_name(software)
            if not backend_supports_geometry_optimization(backend):
                raise MicrosolvationRequestError(f"outer-sphere wall requires a registered energy-gradient provider; {backend!r} has none")
            if self.backend_parameters.get("geometry_optimizer") != "geometric":
                raise MicrosolvationRequestError("outer-sphere wall requires --geometry-optimizer geometric")
        if self.wall == "auto" and self.objective == "shell":
            from pyar.backend_capabilities import backend_supports_geometry_optimization, normalize_backend_name
            software = self.backend_parameters.get("software")
            if not software:
                raise MicrosolvationRequestError("automatic wall policy requires an energy-gradient backend")
            backend = normalize_backend_name(software)
            if not backend_supports_geometry_optimization(backend):
                raise MicrosolvationRequestError(f"automatic wall policy requires a registered energy-gradient provider; {backend!r} has none")
            if self.backend_parameters.get("geometry_optimizer") != "geometric":
                raise MicrosolvationRequestError("automatic wall policy requires --geometry-optimizer geometric")
        if self.placement_score not in {"geometric", "energy"}:
            raise MicrosolvationRequestError("placement score must be 'geometric' or 'energy'")
        if self.placement_score == "energy" and not self.backend_parameters.get("software"):
            raise MicrosolvationRequestError("energy placement scoring requires --backend")
        if self.placement_score == "energy":
            from pyar.backend_capabilities import backend_supports_geometry_optimization, normalize_backend_name
            backend = normalize_backend_name(self.backend_parameters["software"])
            if not backend_supports_geometry_optimization(backend):
                raise MicrosolvationRequestError(
                    f"energy placement scoring requires a registered energy-gradient provider; {backend!r} has none"
                )
        if self.objective == "cluster" and self.placement_score != "energy":
            raise MicrosolvationRequestError("cluster objective requires --placement-score energy")
        if (isinstance(self.resampling_moves, bool) or not isinstance(self.resampling_moves, int)
                or self.resampling_moves < 0):
            raise MicrosolvationRequestError("resampling moves must be a non-negative integer")
        if self.resampling_moves:
            if self.objective != "cluster":
                raise MicrosolvationRequestError("whole-cluster resampling requires --objective cluster")
            if not self.backend_parameters.get("software") or self.backend_parameters.get("geometry_optimizer") != "geometric":
                raise MicrosolvationRequestError("whole-cluster resampling requires an energy-gradient backend and geometric optimizer")
            if (isinstance(self.resampling_evaluation_budget, bool)
                    or not isinstance(self.resampling_evaluation_budget, int) or self.resampling_evaluation_budget < 1
                    or isinstance(self.resampling_time_limit, bool)
                    or not np.isfinite(float(self.resampling_time_limit)) or self.resampling_time_limit <= 0
                    or isinstance(self.resampling_seed, bool) or not isinstance(self.resampling_seed, int)):
                raise MicrosolvationRequestError("resampling requires positive evaluation/time budgets and an integer seed")
        elif (self.resampling_evaluation_budget or self.resampling_time_limit or self.resampling_seed != 0):
            raise MicrosolvationRequestError("resampling budgets require --resampling-moves")
        if self.site is not None:
            site = tuple(int(index) for index in self.site)
            if not site or len(set(site)) != len(site):
                raise MicrosolvationRequestError("--site requires one or more distinct solute atom indices")
            for index in site:
                if index < 0 or index >= len(self.solute):
                    raise MicrosolvationRequestError(
                        f"Solute site index {index} is out of range (valid range: 0..{len(self.solute) - 1})"
                    )
            object.__setattr__(self, "site", site)
        groups = tuple((str(name).strip(), tuple(int(index) for index in indices))
                       for name, indices in self.site_groups)
        if groups:
            names = [name for name, _ in groups]
            if (any(not name.strip() for name in names) or len(names) != len(set(names))
                    or any(not indices or len(indices) != len(set(indices)) for _, indices in groups)):
                raise MicrosolvationRequestError("site groups require unique names and distinct non-empty atom indices")
            for name, indices in groups:
                if any(index < 0 or index >= len(self.solute) for index in indices):
                    raise MicrosolvationRequestError(f"Site group {name!r} contains an out-of-range solute atom index")
            if (isinstance(self.contact_distance, bool) or self.contact_distance is None
                    or not np.isfinite(float(self.contact_distance)) or float(self.contact_distance) <= 0):
                raise MicrosolvationRequestError("site contact reporting requires a positive --contact-distance in Angstrom")
        elif self.contact_distance is not None:
            raise MicrosolvationRequestError("--contact-distance requires at least one --site-group")
        object.__setattr__(self, "site_groups", groups)
        from pyar.molecule_merge import combine_multiplicity
        electrons = sum(self.solute.atomic_number)
        charge = int(self.solute.charge)
        multiplicity = int(self.solute.multiplicity)
        backend = self.backend_parameters.get("software")
        capabilities = None
        if backend:
            from pyar.backend_capabilities import get_backend_capabilities, normalize_backend_name
            backend = normalize_backend_name(backend)
            capabilities = get_backend_capabilities(backend)
        for step in range(self.count + 1):
            if (electrons - charge < multiplicity - 1
                    or (electrons - charge - multiplicity + 1) % 2):
                raise MicrosolvationRequestError(
                    f"Fragment spin-combination rule yields an invalid cluster state after {step} solvent(s)"
                )
            if capabilities and (
                (charge and not capabilities.supports_charge)
                or (multiplicity != 1 and not capabilities.supports_multiplicity)
            ):
                raise MicrosolvationRequestError(
                    f"Backend {backend!r} does not support the combined electronic state after {step} solvent(s)"
                )
            if step < self.count:
                electrons += sum(self.solvent.atomic_number)
                charge += int(self.solvent.charge)
                multiplicity = combine_multiplicity(multiplicity, int(self.solvent.multiplicity))

    def to_state_dict(self):
        result = {
            "solute": _molecule_signature(self.solute),
            "solvent": _molecule_signature(self.solvent),
            "count": self.count,
            "orientations": self.number_of_orientations,
            "backend_parameters": dict(self.backend_parameters),
            "maximum_number_of_seeds": self.maximum_number_of_seeds,
            "solute_site_atoms": None if self.site is None else list(self.site),
            "surface": {
                "model": "vdw-probe-fibonacci-v1",
                "points_per_atom": self.surface_points_per_atom,
                "probe_radius_angstrom": float(self.probe_radius),
            },
            "sampling": {"directions": "fibonacci", "rotations": "halton"},
            "shell_policy": "original-solute-first-shell-postfilter-v1",
            "shell_tolerance_angstrom": float(self.shell_tolerance),
            "selection_policy": {
                "name": "surface-coverage-energy-diversity-v1",
                "deduplication": "adaptive-graph-first-rmsd",
                "coverage_frontier_fraction_tolerance": 0.03,
            },
            "confinement": {"mode": "placement-and-post-optimization-shell-filter"},
        }
        # Keep the historical default request byte-for-byte compatible so
        # existing v1/v2 runs retain their established continuation owner.
        if self.objective != "shell" or self.placement_score != "geometric" or self.wall != "none":
            if self.wall == "auto":
                result["selection_policy"] = {
                    "name": "objective-aware-v2",
                    "objective": self.objective,
                    "placement_score": self.placement_score,
                    "energy_unit": "hartree",
                    "ranking": ("fixed-composition-energy" if self.objective == "cluster"
                                else "solute-surface-coverage-then-physical-energy"),
                    "coverage_frontier_fraction_tolerance": 0.03 if self.objective == "shell" else None,
                    "cluster_energy_frontier_hartree": 0.005 if self.objective == "cluster" else None,
                    "unbiased_energy_frontier_hartree": None,
                    "prescreen_pool_multiplier": 4,
                }
            elif self.wall == "outer-sphere":
                # Preserve the opt-in MS3 selector and fingerprint so interrupted
                # explicit-wall runs retain their existing continuation contract.
                result["selection_policy"] = {
                    "name": "objective-aware-v1",
                    "objective": self.objective,
                    "placement_score": self.placement_score,
                    "energy_unit": "hartree",
                    "ranking": "unbiased-post-relaxation-energy",
                    "cluster_energy_frontier_hartree": None,
                    "unbiased_energy_frontier_hartree": 0.005,
                    "prescreen_pool_multiplier": 4,
                }
            else:
                result["selection_policy"] = {
                    "name": "objective-aware-v1",
                    "objective": self.objective,
                    "placement_score": self.placement_score,
                    "energy_unit": "hartree",
                    "ranking": None,
                    "cluster_energy_frontier_hartree": 0.005 if self.objective == "cluster" else None,
                    "unbiased_energy_frontier_hartree": None,
                    "prescreen_pool_multiplier": 4,
                }
            result["placement_policy"] = {
                "name": "solvent-cluster-growth-v1" if self.objective == "cluster" else "solute-surface-v1",
            }
            if self.objective == "cluster":
                result["shell_policy"] = "reported-not-filtered-cluster-growth-v1"
                result["confinement"] = {"mode": "none; finite-cluster objective"}
        if self.site_groups:
            result["site_contact_policy"] = {
                "name": "explicit-atom-distance-v1",
                "groups": [{"name": name, "solute_atom_indices": list(indices)}
                           for name, indices in self.site_groups],
                "distance_angstrom": float(self.contact_distance),
                "shell_assignment": "direct-solute-contact; then solvent-contact-to-first-shell-v1",
            }
        if self.resampling_moves:
            result["resampling_policy"] = {
                "name": "whole-cluster-perturbation-v1", "moves": int(self.resampling_moves),
                "evaluation_budget": int(self.resampling_evaluation_budget),
                "time_limit_seconds": float(self.resampling_time_limit), "seed": int(self.resampling_seed),
                "rng": "numpy-PCG64", "collision_scale": 0.85,
            }
        if self.wall != "none":
            if self.wall == "outer-sphere":
                result["confinement"] = {
                    "mode": "opt-in-outer-sphere-harmonic-v1",
                    "objective": "shell-search-only",
                    "center": "instantaneous-solute-center-of-mass",
                    "radius_angstrom": float(self.wall_radius),
                    "force_constant_hartree_per_angstrom2": float(self.wall_force_constant),
                    "final_stage": "unbiased-relaxation-and-physical-energy-ranking",
                }
            elif self.objective == "cluster":
                result["confinement"] = {
                    "mode": "none; finite-cluster objective",
                    "requested_policy": "auto-disabled-for-cluster-objective",
                }
            else:
                from pyar.microsolvation.wall_policy import derive_outer_sphere_parameters
                result["confinement"] = {
                    "mode": "case-dependent-outer-sphere-v1",
                    "trigger": "apply-only-after-unbiased-optimization-fails-shell-validity",
                    "objective": "shell-search-only",
                    "center": "instantaneous-solute-center-of-mass",
                    **derive_outer_sphere_parameters(self),
                    "final_stage": "unbiased-relaxation-and-physical-energy-ranking",
                }
        return result
