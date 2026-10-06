"""Validated request model for solute-centred microsolvation."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping

import numpy as np

from pyar.aggregation.request import _molecule_signature


class MicrosolvationRequestError(ValueError):
    """Raised when a microsolvation request is scientifically invalid."""


@dataclass(frozen=True)
class MicrosolvationRequest:
    """Resolved settings for a finite solvent shell around one fixed solute."""

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
        object.__setattr__(self, "backend_parameters", dict(self.backend_parameters or {}))
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
        return {
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
