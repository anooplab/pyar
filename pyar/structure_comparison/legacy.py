"""Adapters around PyAR's current identity and geometry algorithms."""

from __future__ import annotations

from collections import Counter

import numpy as np

import pyar.representations
from pyar.structure_comparison.models import ComparisonResult
from pyar.structure_comparison.rmsd import rmsd_after_alignment


class LegacyRMSDComparator:
    """Use the existing Coulomb prefilter and permutation/Kabsch RMSD."""

    method = "coulomb-fingerprint-permutation-kabsch-rmsd"

    def __init__(self, threshold=None, fingerprint_distance=None):
        self.threshold = threshold
        self._fingerprint_distance = fingerprint_distance or self._legacy_fingerprint_distance

    @staticmethod
    def _legacy_fingerprint_distance(first, second):
        return np.linalg.norm(
            pyar.representations.fingerprint(first.atoms_list, first.coordinates)
            - pyar.representations.fingerprint(second.atoms_list, second.coordinates)
        )

    def compare(self, first, second) -> ComparisonResult:
        compatible = (
            len(first.atoms_list) == len(second.atoms_list)
            and Counter(first.atoms_list) == Counter(second.atoms_list)
        )
        if not compatible:
            return ComparisonResult(False, None, None, self.method, self.threshold)
        fingerprint_distance = self._fingerprint_distance(first, second)
        if not abs(fingerprint_distance) < 1.0:
            return ComparisonResult(
                True, None, False, self.method, self.threshold,
                {"fingerprint_distance": float(fingerprint_distance), "prefiltered": True},
            )
        distance = rmsd_after_alignment(first, second)
        equivalent = None if self.threshold is None else distance < self.threshold
        return ComparisonResult(
            True, distance, equivalent, self.method, self.threshold,
            {"fingerprint_distance": float(fingerprint_distance)},
        )


__all__ = ["LegacyRMSDComparator"]
