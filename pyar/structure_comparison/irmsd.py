"""Optional adapter for the permutation-invariant iRMSD package."""

from __future__ import annotations

from collections import Counter

import numpy as np

from pyar.structure_comparison.models import ComparisonResult


class IRMSDComparator:
    """Compare equal-formula geometries with the optional ``irmsd`` package.

    iRMSD handles atom permutation and alignment but does not prove molecular
    graph identity. Use :class:`GraphRMSDComparator` when connectivity must be
    part of the compatibility gate.
    """

    method = "irmsd-element-permutation-rmsd"

    def __init__(self, threshold=None, inversion=2):
        if inversion not in (0, 1, 2):
            raise ValueError("inversion must be 0 (auto), 1 (on), or 2 (off)")
        self.threshold = threshold
        self.inversion = inversion

    def compare(self, first, second) -> ComparisonResult:
        compatible = (
            len(first.atoms_list) == len(second.atoms_list)
            and Counter(first.atoms_list) == Counter(second.atoms_list)
        )
        if not compatible:
            return ComparisonResult(False, None, None, self.method, self.threshold)

        try:
            import irmsd
        except ImportError as exc:
            raise ImportError(
                "IRMSDComparator requires the optional dependency; install "
                "pyar-chem[structure-comparison]"
            ) from exc

        molecule_type = irmsd.Molecule
        first_ir = molecule_type(
            symbols=list(first.atoms_list), positions=np.asarray(first.coordinates, dtype=float),
        )
        second_ir = molecule_type(
            symbols=list(second.atoms_list), positions=np.asarray(second.coordinates, dtype=float),
        )
        distance, _, _ = irmsd.get_irmsd_molecule(
            first_ir, second_ir, iinversion=self.inversion,
        )
        distance = float(distance)
        equivalent = None if self.threshold is None else distance < self.threshold
        return ComparisonResult(
            True, distance, equivalent, self.method, self.threshold,
            {"inversion": self.inversion, "connectivity_checked": False},
        )


__all__ = ["IRMSDComparator"]
