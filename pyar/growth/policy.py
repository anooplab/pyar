"""Immutable operation settings shared by aggregation and fixed-addend growth."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from types import MappingProxyType
from typing import Any, Mapping


@dataclass(frozen=True)
class SelectionOptions:
    feature: str = "auto"
    algorithm: str = "auto"
    distance: str = "euclidean"
    system_type: str = "auto"

    @classmethod
    def from_request(cls, request):
        return cls(
            request.selection_feature,
            request.selection_algorithm,
            request.selection_distance,
            request.selection_system_type,
        )


@dataclass(frozen=True)
class GrowthPolicy:
    orientations: int = 1
    qc_params: Mapping[str, Any] = field(default_factory=dict)
    maximum_number_of_seeds: int = 12
    site: object = None
    connectivity_policy: str | None = None
    selection: SelectionOptions = field(default_factory=SelectionOptions)
    bound_geometry_only: bool = True

    def __post_init__(self):
        # Freeze ownership of the outer settings mapping; nested backend values
        # retain their established types. Backends receive ordinary fresh dicts.
        object.__setattr__(self, "qc_params", MappingProxyType(dict(self.qc_params)))

    @classmethod
    def from_aggregate_request(cls, request, orientations):
        return cls(
            orientations,
            request.backend_parameters,
            request.maximum_number_of_seeds,
            None if request.site is None else list(request.site),
            request.connectivity_policy,
            SelectionOptions.from_request(request),
        )

    @classmethod
    def from_grow_request(cls, request, qc_params, site):
        return cls(
            request.number_of_orientations,
            qc_params,
            request.maximum_number_of_seeds,
            site,
            request.connectivity_policy,
            SelectionOptions.from_request(request),
        )


@dataclass(frozen=True)
class OptimizationStages:
    loose: dict
    refine: dict
    staged: bool


@dataclass(frozen=True)
class GrowthRun:
    aggregate_id: str
    seeds: list
    monomer: object
    policy: GrowthPolicy
    directory: Path
    connectivity_policy: str
    stages: OptimizationStages


@dataclass(frozen=True)
class SeedRun:
    run: GrowthRun
    index: int
    molecule: object

    @property
    def directory(self):
        return self.run.directory / f"seed_{self.index:03d}"

    @property
    def molecule_id(self):
        return f"{self.index:03d}_{self.run.aggregate_id}"


@dataclass(frozen=True)
class SeedOptimization:
    accepted: list
    retry: list
    statuses: list


@dataclass(frozen=True)
class FinalSelectionRun:
    directory: Path = Path(".")
    source_records: Mapping | None = None
    energy_ranked: bool = True


@dataclass(frozen=True)
class FinalGroupSelection:
    stoichiometry: str
    sources: list
    selected: list
    diagnostics: dict
