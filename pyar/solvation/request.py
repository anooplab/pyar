"""Typed request model for solvation workflows."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping

import numpy as np


class SolvationRequestError(ValueError):
    """Raised when a solvation request is invalid."""


def _molecule_signature(molecule):
    """Return stable molecule metadata used to validate restarts."""
    return {
        "name": molecule.name,
        "atoms": list(molecule.atoms_list),
        "coordinates": np.asarray(molecule.coordinates, dtype=float).tolist(),
        "charge": molecule.charge,
        "multiplicity": molecule.multiplicity,
        "scftype": molecule.scftype,
        "fragment_definition": list(getattr(molecule, "fragments", [])),
    }


def normalize_connectivity_policy(connectivity_policy):
    """Return the enforced connectivity policy for solvation workflows."""
    normalized_policy = "auto" if connectivity_policy is None else str(connectivity_policy).lower()
    if normalized_policy in {"auto", "off"}:
        return "off"
    if normalized_policy in {"prefer", "strict"}:
        raise SolvationRequestError(
            "Solvation workflows require connectivity policy 'off'; "
            f"received {connectivity_policy!r}."
        )
    raise SolvationRequestError(
        f"Unknown connectivity policy: {connectivity_policy!r}. "
        "Expected one of 'auto', 'off', 'prefer', or 'strict'."
    )


@dataclass(frozen=True)
class SolvationRequest:
    """Validated options for one solvation workflow run."""

    seeds: tuple[Any, ...]
    monomer: Any
    aggregate_size: int
    orientations: Any
    backend_parameters: Mapping[str, Any] = field(default_factory=dict)
    maximum_number_of_seeds: int = 1
    site: tuple[Any, ...] | None = None
    connectivity_policy: str = "off"
    seed_signatures: tuple[Mapping[str, Any], ...] = field(default_factory=tuple)
    monomer_signature: Mapping[str, Any] = field(default_factory=dict)

    @classmethod
    def from_options(
        cls,
        seeds,
        monomer,
        aggregate_size,
        hm_orientations,
        qc_params,
        maximum_number_of_seeds,
        site,
        connectivity_policy,
    ):
        """Build a normalized request from public workflow arguments."""
        seeds = tuple(seeds or ())
        if not seeds:
            raise SolvationRequestError("Solvation requires at least one seed molecule")
        if monomer is None:
            raise SolvationRequestError("Solvation requires one solvent molecule")

        aggregate_size = int(aggregate_size)
        if aggregate_size < 1:
            raise SolvationRequestError("--solvation-size must be at least 1")

        maximum_number_of_seeds = int(maximum_number_of_seeds)
        if maximum_number_of_seeds < 1:
            raise SolvationRequestError("--maximum-number-of-seeds must be at least 1")

        normalized_site = None if site is None else tuple(site)
        return cls(
            seeds=seeds,
            monomer=monomer,
            aggregate_size=aggregate_size,
            orientations=hm_orientations,
            backend_parameters=dict(qc_params or {}),
            maximum_number_of_seeds=maximum_number_of_seeds,
            site=normalized_site,
            connectivity_policy=normalize_connectivity_policy(connectivity_policy),
            seed_signatures=tuple(_molecule_signature(seed) for seed in seeds),
            monomer_signature=_molecule_signature(monomer),
        )

    def to_state_dict(self):
        """Return the JSON-serializable representation stored in state files."""
        return {
            "aggregate_size": self.aggregate_size,
            "orientations": self.orientations,
            "backend_parameters": dict(self.backend_parameters),
            "maximum_number_of_seeds": self.maximum_number_of_seeds,
            "site": None if self.site is None else list(self.site),
            "connectivity_policy": self.connectivity_policy,
            "seeds": [dict(seed) for seed in self.seed_signatures],
            "monomer": dict(self.monomer_signature),
        }
