"""Coulomb eigenvalue prefilter with permutation-aware Kabsch RMSD."""

from __future__ import annotations

from collections import Counter

import numpy as np

import pyar.representations
from pyar.structure_comparison.models import ComparisonResult
from pyar.structure_comparison.rmsd import rmsd_after_alignment


class CoulombEigenvalueRMSDComparator:
    """Use sorted Coulomb eigenvalues as a prefilter before aligned RMSD."""

    method = "coulomb-eigenvalue-prefilter-permutation-kabsch-rmsd"

    def __init__(self, threshold=None, coulomb_eigenvalue_distance=None):
        self.threshold = threshold
        self._eigenvalue_distance = coulomb_eigenvalue_distance or self._coulomb_eigenvalue_distance

    @staticmethod
    def _coulomb_eigenvalue_distance(first, second):
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
        eigenvalue_distance = self._eigenvalue_distance(first, second)
        if not abs(eigenvalue_distance) < 1.0:
            return ComparisonResult(
                True, None, False, self.method, self.threshold,
                {"coulomb_eigenvalue_distance": float(eigenvalue_distance), "prefiltered": True},
            )
        distance = rmsd_after_alignment(first, second)
        equivalent = None if self.threshold is None else distance < self.threshold
        return ComparisonResult(
            True, distance, equivalent, self.method, self.threshold,
            {"coulomb_eigenvalue_distance": float(eigenvalue_distance)},
        )


__all__ = ["CoulombEigenvalueRMSDComparator"]
