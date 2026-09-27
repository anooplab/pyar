"""Geometrical comparison service interfaces."""

from __future__ import annotations

from typing import Protocol, runtime_checkable

from pyar.structure_comparison.models import ComparisonResult


@runtime_checkable
class StructureComparator(Protocol):
    """Service that compares geometries when their chemistry is compatible."""

    def compare(self, first, second) -> ComparisonResult:
        """Return compatibility, distance, and optional equivalence decision."""
