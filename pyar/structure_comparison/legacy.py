"""Adapters around PyAR's current identity and geometry algorithms."""

from __future__ import annotations

from collections import Counter

from pyar.structure_comparison.models import ComparisonResult


class LegacyRMSDComparator:
    """Use the existing Coulomb prefilter and permutation/Kabsch RMSD."""

    method = "coulomb-fingerprint-permutation-kabsch-rmsd"

    def __init__(self, threshold=None):
        self.threshold = threshold

    def compare(self, first, second) -> ComparisonResult:
        from pyar.selection import deduplication

        compatible = (
            len(first.atoms_list) == len(second.atoms_list)
            and Counter(first.atoms_list) == Counter(second.atoms_list)
        )
        if not compatible:
            return ComparisonResult(False, None, None, self.method, self.threshold)
        from pyar.selection import clustering

        fingerprint_distance = clustering.calc_fingerprint_distance(first, second)
        if not abs(fingerprint_distance) < 1.0:
            return ComparisonResult(
                True, None, False, self.method, self.threshold,
                {"fingerprint_distance": float(fingerprint_distance), "prefiltered": True},
            )
        distance = deduplication._rmsd_after_alignment(first, second)
        equivalent = None if self.threshold is None else distance < self.threshold
        return ComparisonResult(
            True, distance, equivalent, self.method, self.threshold,
            {"fingerprint_distance": float(fingerprint_distance)},
        )


__all__ = ["LegacyRMSDComparator"]
