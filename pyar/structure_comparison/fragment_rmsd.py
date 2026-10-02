"""Bounded whole-fragment matching under one global proper rotation.

Fragment connectivity comes from coordinates, independently for each geometry.
No original growth-fragment assignment is trusted after chemical reactions.
Distances are valid mapped RMSD upper bounds, not guaranteed global minima.
"""

from itertools import permutations, product
from math import factorial, prod

import networkx as nx
import numpy as np

from pyar.structure_comparison.graph_rmsd import infer_molecular_graph, selected_atom_indices
from pyar.structure_comparison.models import ComparisonResult
from pyar.structure_comparison.rmsd import kabsch_rmsd, kabsch_rotation


class FragmentRMSDComparator:
    """Match isomorphic fragments, preserving their placement and orientation."""

    method = "fragment-matched-global-kabsch-rmsd-upper-bound"

    def __init__(self, *, atom_mode="all", bond_scale=1.15, max_mappings=10000,
                 max_fragment_permutations=720):
        if atom_mode not in {"all", "heavy"}:
            raise ValueError("atom_mode must be 'all' or 'heavy'")
        if not np.isfinite(bond_scale) or bond_scale <= 0:
            raise ValueError("bond_scale must be finite and positive")
        for value in (max_mappings, max_fragment_permutations):
            if isinstance(value, bool) or not isinstance(value, int) or value < 1:
                raise ValueError("Fragment matching limits must be positive integers")
        self.atom_mode = atom_mode
        self.bond_scale = bond_scale
        self.max_mappings = max_mappings
        self.max_fragment_permutations = max_fragment_permutations

    def compare(self, first, second):
        metadata = {"atom_mode": self.atom_mode, "bond_scale": self.bond_scale,
                    "distance_is_upper_bound": True}
        for attribute in ("charge", "multiplicity"):
            left, right = getattr(first, attribute, None), getattr(second, attribute, None)
            if left is not None and right is not None and left != right:
                return ComparisonResult(False, None, None, self.method, metadata=metadata)
        graphs = [infer_molecular_graph(molecule, self.bond_scale) for molecule in (first, second)]
        components = [[sorted(component) for component in nx.connected_components(graph)] for graph in graphs]
        if (len(components[0]) < 2 or len(components[0]) != len(components[1])
                or any(len(component) < 2 for parts in components for component in parts)):
            return ComparisonResult(False, None, None, self.method,
                                    metadata={**metadata, "reason": "requires matching multi-atom fragments"})

        # Classify fragments by exact element-labelled isomorphism, never by formula alone.
        representatives = []
        groups = [[], []]
        for side in range(2):
            for component in components[side]:
                graph = graphs[side].subgraph(component)
                for group_id, representative in enumerate(representatives):
                    if nx.is_isomorphic(graph, representative, node_match=lambda a, b: a["element"] == b["element"]):
                        break
                else:
                    group_id = len(representatives)
                    representatives.append(graph)
                groups[side].append(group_id)
        group_ids = sorted(set(groups[0]))
        if sorted(groups[0]) != sorted(groups[1]):
            return ComparisonResult(False, None, None, self.method,
                                    metadata={**metadata, "reason": "fragment topologies differ"})
        permutation_count = prod(factorial(groups[0].count(group)) for group in group_ids)
        if permutation_count > self.max_fragment_permutations:
            return ComparisonResult(True, None, None, self.method,
                                    metadata={**metadata, "reason": "fragment permutation limit", "comparison_complete": False})

        # Precompute each allowed internal graph mapping with a bounded budget.
        mappings = {}
        total = 0
        for left_index, left in enumerate(components[0]):
            for right_index, right in enumerate(components[1]):
                if groups[0][left_index] != groups[1][right_index]:
                    continue
                matcher = nx.algorithms.isomorphism.GraphMatcher(
                    graphs[0].subgraph(left), graphs[1].subgraph(right),
                    node_match=lambda a, b: a["element"] == b["element"],
                )
                orders = []
                for mapping in matcher.isomorphisms_iter():
                    total += 1
                    if total > self.max_mappings:
                        return ComparisonResult(True, None, None, self.method,
                                                metadata={**metadata, "reason": "internal mapping limit", "comparison_complete": False})
                    orders.append([mapping[index] for index in left])
                mappings[left_index, right_index] = orders

        coordinates = [np.asarray(molecule.coordinates, dtype=float) for molecule in (first, second)]
        centered = [value - value.mean(axis=0) for value in coordinates]
        selected = selected_atom_indices(first.atoms_list, self.atom_mode)
        left_groups = [[i for i, group in enumerate(groups[0]) if group == target] for target in group_ids]
        right_groups = [[i for i, group in enumerate(groups[1]) if group == target] for target in group_ids]
        best = None
        evaluated = 0
        for group_orders in product(*(permutations(group) for group in right_groups)):
            assignment = dict(pair for left, right in zip(left_groups, group_orders) for pair in zip(left, right))
            # Anchor fits cover the otherwise undetermined global orientation of
            # a dimer's centre-to-centre axis. Other fragments are never fitted
            # independently when the final placement is scored.
            for anchor, left in enumerate(components[0]):
                for anchor_order in mappings[anchor, assignment[anchor]]:
                    evaluated += 1
                    if evaluated > self.max_mappings:
                        return ComparisonResult(True, None, None, self.method,
                                                metadata={**metadata, "reason": "rotation hypothesis limit", "comparison_complete": False})
                    first_anchor = centered[0][left]
                    second_anchor = centered[1][anchor_order]
                    rotation = kabsch_rotation(first_anchor - first_anchor.mean(axis=0),
                                               second_anchor - second_anchor.mean(axis=0))
                    order = np.empty(len(coordinates[0]), dtype=int)
                    for fragment, indices in enumerate(components[0]):
                        options = mappings[fragment, assignment[fragment]]
                        scored_indices = [index for index in indices if index in selected]
                        offsets = [indices.index(index) for index in scored_indices]
                        local_order = min(options, key=lambda option: np.sum(
                            (centered[0][scored_indices] - centered[1][np.asarray(option)[offsets]] @ rotation) ** 2))
                        order[indices] = local_order
                    distance = kabsch_rmsd(coordinates[0][selected], coordinates[1][order[selected]])
                    best = distance if best is None else min(best, distance)
        metadata.update(fragment_count=len(components[0]), rotation_hypotheses=evaluated,
                        comparison_complete=True)
        return ComparisonResult(True, float(best), None, self.method, metadata=metadata)


__all__ = ["FragmentRMSDComparator"]
