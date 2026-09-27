"""Shared immutable result models for structure comparison services."""

from __future__ import annotations

from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Any, Mapping


@dataclass(frozen=True)
class IdentityResult:
    """Chemical identity data, independent of any one identifier format."""

    canonical_key: str
    method: str
    inchi: str | None = None
    smiles: str | None = None
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        object.__setattr__(self, "metadata", MappingProxyType(dict(self.metadata)))


@dataclass(frozen=True)
class ComparisonResult:
    """Outcome of a geometrical comparison between two structures."""

    compatible: bool
    distance: float | None
    equivalent: bool | None
    method: str
    threshold: float | None = None
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        object.__setattr__(self, "metadata", MappingProxyType(dict(self.metadata)))
