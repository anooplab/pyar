"""Typed request model for aggregation workflows."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping

import numpy as np


class AggregateRequestError(ValueError):
    """Raised when an aggregation request is invalid."""


def _molecule_signature(molecule):
    """Return stable input geometry metadata used to validate restarts."""
    return {
        "atoms": list(molecule.atoms_list),
        "coordinates": np.asarray(molecule.coordinates, dtype=float).tolist(),
        "charge": molecule.charge,
        "multiplicity": molecule.multiplicity,
        "scftype": molecule.scftype,
        "fragment_definition": list(getattr(molecule, "fragments", [])),
    }


def _normalize_connectivity_policy(connectivity_policy):
    """Return the persisted connectivity-policy request value."""
    normalized = "auto" if connectivity_policy is None else str(connectivity_policy).lower()
    if normalized not in {"auto", "off", "prefer", "strict"}:
        raise AggregateRequestError(
            f"Unknown connectivity policy: {connectivity_policy!r}. "
            "Expected one of 'auto', 'off', 'prefer', or 'strict'."
        )
    return normalized


@dataclass(frozen=True)
class AggregateRequest:
    """Validated options for one aggregation workflow run."""

    molecules: tuple[Any, ...]
    aggregate_sizes: tuple[int, ...]
    orientations: Any
    backend_parameters: Mapping[str, Any] = field(default_factory=dict)
    maximum_number_of_seeds: int = 1
    first_pathway: int = 0
    number_of_pathways: int = 1
    site: tuple[Any, ...] | None = None
    connectivity_policy: str = "auto"
    fragments: tuple[Mapping[str, Any], ...] = field(default_factory=tuple)

    @classmethod
    def from_options(
        cls,
        molecules,
        aggregate_sizes,
        hm_orientations,
        qc_params,
        maximum_number_of_seeds,
        first_pathway,
        number_of_pathways,
        site,
        connectivity_policy,
    ):
        """Build a normalized request from public workflow arguments."""
        molecules = tuple(molecules or ())
        if not molecules:
            raise AggregateRequestError("Aggregation requires at least one molecule")

        aggregate_sizes = tuple(int(size) for size in (aggregate_sizes or ()))
        if len(aggregate_sizes) != len(molecules):
            raise AggregateRequestError("Aggregate sizes must be specified for every molecule")
        if any(size < 1 for size in aggregate_sizes):
            raise AggregateRequestError("Aggregate sizes must be positive integers")

        maximum_number_of_seeds = int(maximum_number_of_seeds)
        if maximum_number_of_seeds < 1:
            raise AggregateRequestError("--maximum-number-of-seeds must be at least 1")

        first_pathway = int(first_pathway)
        if first_pathway < 0:
            raise AggregateRequestError("--first-pathway must be non-negative")

        number_of_pathways = int(number_of_pathways)
        if number_of_pathways < 1:
            raise AggregateRequestError("--number-of-pathways must be at least 1")

        normalized_site = None if site is None else tuple(site)
        return cls(
            molecules=molecules,
            aggregate_sizes=aggregate_sizes,
            orientations=hm_orientations,
            backend_parameters=dict(qc_params or {}),
            maximum_number_of_seeds=maximum_number_of_seeds,
            first_pathway=first_pathway,
            number_of_pathways=number_of_pathways,
            site=normalized_site,
            connectivity_policy=_normalize_connectivity_policy(connectivity_policy),
            fragments=tuple(_molecule_signature(molecule) for molecule in molecules),
        )

    def to_state_dict(self):
        """Return the JSON-serializable representation stored in state files."""
        return {
            "aggregate_sizes": list(self.aggregate_sizes),
            "orientations": self.orientations,
            "backend_parameters": dict(self.backend_parameters),
            "maximum_number_of_seeds": self.maximum_number_of_seeds,
            "first_pathway": self.first_pathway,
            "number_of_pathways": self.number_of_pathways,
            "site": None if self.site is None else list(self.site),
            "connectivity_policy": self.connectivity_policy,
            "fragments": [dict(fragment) for fragment in self.fragments],
        }
