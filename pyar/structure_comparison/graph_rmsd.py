"""Connectivity-constrained permutation and Kabsch RMSD comparison."""

from __future__ import annotations

from collections import Counter

import networkx as nx
import numpy as np

from pyar.structure_comparison.models import ComparisonResult
from pyar.structure_comparison.rmsd import kabsch_rmsd
from pyar.structure_comparison.coordinate_graph import infer_coordinate_graph


def infer_molecular_graph(molecule, bond_scale=1.15):
    """Compatibility wrapper for the shared coordinate-only graph inference.

    This deliberately exposes the bond scale: XYZ has no bond orders, and
    distance-based connectivity is a model choice. Nodes retain element labels.
    """
    return infer_coordinate_graph(molecule, model="covalent-radii", scale=bond_scale)


class GraphRMSDComparator:
    """Only calculate RMSD when inferred molecular graphs are isomorphic."""

    method = "element-labeled-graph-isomorphism-kabsch-rmsd"

    def __init__(self, threshold=None, bond_scale=1.15, max_isomorphisms=10000):
        if bond_scale <= 0:
            raise ValueError("bond_scale must be positive")
        if not isinstance(max_isomorphisms, int) or max_isomorphisms <= 0:
            raise ValueError("max_isomorphisms must be a positive integer")
        self.threshold = threshold
        self.bond_scale = float(bond_scale)
        self.max_isomorphisms = max_isomorphisms

    def compare(self, first, second) -> ComparisonResult:
        if (
            len(first.atoms_list) != len(second.atoms_list)
            or Counter(first.atoms_list) != Counter(second.atoms_list)
        ):
            return ComparisonResult(False, None, None, self.method, self.threshold)

        first_graph = infer_molecular_graph(first, self.bond_scale)
        second_graph = infer_molecular_graph(second, self.bond_scale)
        matcher = nx.algorithms.isomorphism.GraphMatcher(
            first_graph, second_graph,
            node_match=lambda left, right: left["element"] == right["element"],
        )
        best = None
        isomorphisms_evaluated = 0
        first_coords = np.asarray(first.coordinates, dtype=float)
        second_coords = np.asarray(second.coordinates, dtype=float)
        for mapping in matcher.isomorphisms_iter():
            isomorphisms_evaluated += 1
            if isomorphisms_evaluated > self.max_isomorphisms:
                return ComparisonResult(
                    True, None, None, self.method, self.threshold,
                    {
                        "bond_scale": self.bond_scale,
                        "connectivity_checked": True,
                        "connectivity_match": True,
                        "comparison_complete": False,
                        "isomorphism_limit": self.max_isomorphisms,
                        "isomorphisms_evaluated": self.max_isomorphisms,
                    },
                )
            order = [mapping[index] for index in range(len(first_coords))]
            distance = kabsch_rmsd(first_coords, second_coords[order])
            best = distance if best is None else min(best, distance)

        metadata = {
            "bond_scale": self.bond_scale,
            "connectivity_checked": True,
            "connectivity_match": best is not None,
            "comparison_complete": True,
            "isomorphisms_evaluated": isomorphisms_evaluated,
        }
        if best is None:
            return ComparisonResult(False, None, None, self.method, self.threshold, metadata)
        equivalent = None if self.threshold is None else best < self.threshold
        return ComparisonResult(True, float(best), equivalent, self.method, self.threshold, metadata)


__all__ = ["GraphRMSDComparator", "infer_molecular_graph"]
