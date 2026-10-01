"""Typed request model for aggregation workflows."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping

import numpy as np
from pyar.selection.distances import DISTANCE_METRICS
from pyar.selection.policy import SELECTION_POLICY_VERSION, normalize_system_type


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


def _normalize_selection_feature(feature):
    normalized = str(feature or "auto").strip().lower().replace("_", "-")
    if normalized == "histogram":
        normalized = "distance-histogram"
    if normalized not in {"auto", "mbtr", "soap", "distance-histogram"}:
        raise AggregateRequestError("Unknown selection feature. Choose auto, mbtr, soap, or distance-histogram.")
    return normalized


def _normalize_selection_system_type(system_type):
    try:
        return normalize_system_type(system_type)
    except ValueError as exc:
        raise AggregateRequestError(str(exc)) from exc


def _normalize_selection_algorithm(algorithm):
    normalized = str(algorithm or "auto").strip().lower()
    if normalized not in {"auto", "hybrid", "hdbscan", "agglomerative", "dbscan", "optics", "maxmin", "max-min", "max_min"}:
        raise AggregateRequestError("Unknown selection algorithm.")
    return normalized


def _normalize_distance_metric(metric):
    normalized = str(metric or "euclidean").strip().lower()
    if normalized not in DISTANCE_METRICS:
        raise AggregateRequestError(
            f"Unknown selection distance {metric!r}. Choose one of: {', '.join(DISTANCE_METRICS)}"
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
    selection_feature: str = "auto"
    selection_algorithm: str = "auto"
    selection_distance: str = "euclidean"
    selection_system_type: str = "auto"
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
        selection_feature="auto",
        selection_algorithm="auto",
        selection_distance="euclidean",
        selection_system_type="auto",
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
            selection_feature=_normalize_selection_feature(selection_feature),
            selection_algorithm=_normalize_selection_algorithm(selection_algorithm),
            selection_distance=_normalize_distance_metric(selection_distance),
            selection_system_type=_normalize_selection_system_type(selection_system_type),
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
            "selection_feature": self.selection_feature,
            "selection_algorithm": self.selection_algorithm,
            "selection_distance": self.selection_distance,
            "selection_system_type": self.selection_system_type,
            "selection_policy_version": SELECTION_POLICY_VERSION,
            "fragments": [dict(fragment) for fragment in self.fragments],
        }
