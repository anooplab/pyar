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


def selected_atom_indices(symbols, atom_mode):
    """Select RMSD atoms; hydrogen-only systems always use every atom."""
    heavy = [index for index, symbol in enumerate(symbols) if str(symbol).upper() != "H"]
    return heavy if atom_mode == "heavy" and heavy else list(range(len(symbols)))


def _mapping_graph(graph, atom_mode):
    """Contract terminal hydrogens without losing their connectivity constraints.

    Hydrogen permutations do not affect heavy-atom RMSD. Contracting only leaf
    hydrogens avoids factorial enumeration while retaining bridging hydrogens
    and hydrogen-only components in the mapping graph.
    """
    graph = graph.copy()
    for node in graph:
        graph.nodes[node]["terminal_hydrogens"] = 0
    if atom_mode == "heavy" and any(data["element"] != "H" for _, data in graph.nodes(data=True)):
        leaves = [node for node, data in graph.nodes(data=True)
                  if data["element"] == "H" and graph.degree[node] == 1
                  and graph.nodes[next(iter(graph.neighbors(node)))]["element"] != "H"]
        for node in leaves:
            parent = next(iter(graph.neighbors(node)))
            graph.nodes[parent]["terminal_hydrogens"] += 1
            graph.remove_node(node)
    return graph


class GraphRMSDComparator:
    """Only calculate RMSD when inferred molecular graphs are isomorphic."""

    method = "element-labeled-graph-isomorphism-kabsch-rmsd"

    def __init__(self, threshold=None, bond_scale=1.15, max_isomorphisms=10000,
                 atom_mode="all"):
        if not np.isfinite(bond_scale) or bond_scale <= 0:
            raise ValueError("bond_scale must be finite and positive")
        if threshold is not None and (not np.isfinite(threshold) or threshold < 0):
            raise ValueError("threshold must be finite and nonnegative")
        if not isinstance(max_isomorphisms, int) or max_isomorphisms <= 0:
            raise ValueError("max_isomorphisms must be a positive integer")
        self.threshold = threshold
        self.bond_scale = float(bond_scale)
        self.max_isomorphisms = max_isomorphisms
        if atom_mode not in {"all", "heavy"}:
            raise ValueError("atom_mode must be 'all' or 'heavy'")
        self.atom_mode = atom_mode

    def compare(self, first, second) -> ComparisonResult:
        for attribute in ("charge", "multiplicity"):
            left, right = getattr(first, attribute, None), getattr(second, attribute, None)
            if left is not None and right is not None and left != right:
                return ComparisonResult(
                    False, None, None, self.method, self.threshold,
                    {"electronic_state_match": False, "mismatch_attribute": attribute},
                )
        if (
            not len(first.atoms_list) or len(first.atoms_list) != len(second.atoms_list)
            or Counter(first.atoms_list) != Counter(second.atoms_list)
        ):
            return ComparisonResult(False, None, None, self.method, self.threshold)

        first_graph = _mapping_graph(infer_molecular_graph(first, self.bond_scale), self.atom_mode)
        second_graph = _mapping_graph(infer_molecular_graph(second, self.bond_scale), self.atom_mode)
        matcher = nx.algorithms.isomorphism.GraphMatcher(
            first_graph, second_graph,
            node_match=lambda left, right: (left["element"], left["terminal_hydrogens"])
            == (right["element"], right["terminal_hydrogens"]),
        )
        best = None
        isomorphisms_evaluated = 0
        first_coords = np.asarray(first.coordinates, dtype=float)
        second_coords = np.asarray(second.coordinates, dtype=float)
        indices = selected_atom_indices(first.atoms_list, self.atom_mode)
        for mapping in matcher.isomorphisms_iter():
            isomorphisms_evaluated += 1
            if isomorphisms_evaluated > self.max_isomorphisms:
                return ComparisonResult(
                    True, None, None, self.method, self.threshold,
                    {
                        "bond_scale": self.bond_scale,
                        "atom_mode": self.atom_mode,
                        "connectivity_checked": True,
                        "connectivity_match": True,
                        "comparison_complete": False,
                        "isomorphism_limit": self.max_isomorphisms,
                        "isomorphisms_evaluated": self.max_isomorphisms,
                    },
                )
            order = [mapping[index] for index in indices]
            distance = kabsch_rmsd(first_coords[indices], second_coords[order])
            best = distance if best is None else min(best, distance)

        metadata = {
            "bond_scale": self.bond_scale,
            "atom_mode": self.atom_mode,
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
