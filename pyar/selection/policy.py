"""System-aware defaults for exploratory structural clustering.

Classification uses workflow hints when supplied and otherwise makes only
conservative inferences from XYZ-only adjacency. It does not treat formula or
geometric graphs as definitive chemical identity assignments.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass

import networkx as nx
import numpy as np

from pyar.structure_comparison.coordinate_graph import infer_coordinate_graph

SYSTEM_TYPES = (
    "auto",
    "molecular",
    "isomers",
    "atomic-cluster",
    "molecular-aggregate",
    "unknown",
)

# SOAP is the planned first representation for atomic/molecular aggregates.
# MBTR remains the general molecular and unknown-system representation; it
# includes pair and angular terms in PyAR's current implementation.
FEATURE_POLICY = {
    "molecular": "mbtr",
    "isomers": "mbtr",
    "atomic-cluster": "soap",
    "molecular-aggregate": "soap",
    "unknown": "mbtr",
}
FEATURE_FALLBACK_POLICY = {
    "molecular": ("soap", "distance-histogram"),
    "isomers": ("soap", "distance-histogram"),
    "atomic-cluster": ("mbtr", "distance-histogram"),
    "molecular-aggregate": ("mbtr", "distance-histogram"),
    "unknown": ("soap", "distance-histogram"),
}
DEFAULT_CLUSTER_ALGORITHM = "hybrid"
SELECTION_POLICY_VERSION = 2

_METAL_AND_NOBLE_CLUSTER_ELEMENTS = {
    "Li", "Be", "Na", "Mg", "Al", "K", "Ca", "Sc", "Ti", "V", "Cr",
    "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Rb", "Sr", "Y", "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Cs",
    "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy",
    "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir",
    "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "Fr", "Ra", "Ac", "Th",
    "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md",
    "No", "Lr", "He", "Ne", "Ar", "Kr", "Xe", "Rn",
}


@dataclass(frozen=True)
class SystemClassification:
    """Inferred or user-specified system class with an auditable rationale."""

    system_type: str
    confidence: str
    reason: str
    topology_group_ids: tuple[int, ...] = ()

    def to_dict(self):
        return {
            "system_type": self.system_type,
            "confidence": self.confidence,
            "reason": self.reason,
            "topology_group_ids": list(self.topology_group_ids),
        }


def _graphs_isomorphic(left, right):
    matcher = nx.algorithms.isomorphism.GraphMatcher(
        left,
        right,
        node_match=lambda first, second: first["element"] == second["element"],
    )
    return matcher.is_isomorphic()


def _topology_groups(molecules, graphs=None):
    representatives = []
    group_ids = []
    if graphs is None:
        graphs = [infer_coordinate_graph(molecule) for molecule in molecules]
    for graph in graphs:
        for group_id, representative in enumerate(representatives):
            if _graphs_isomorphic(graph, representative):
                group_ids.append(group_id)
                break
        else:
            group_ids.append(len(representatives))
            representatives.append(graph)
    return tuple(group_ids)


def normalize_system_type(system_type="auto"):
    """Normalize the shared API and CLI system-type vocabulary."""
    hint = str(system_type or "auto").strip().lower().replace("_", "-")
    aliases = {
        "molecule": "molecular", "conformer": "molecular", "conformers": "molecular",
        "aggregate": "molecular-aggregate", "noncovalent-aggregate": "molecular-aggregate",
        "cluster": "atomic-cluster",
    }
    hint = aliases.get(hint, hint)
    if hint not in SYSTEM_TYPES:
        raise ValueError(f"Unknown system type {system_type!r}. Choose one of: {', '.join(SYSTEM_TYPES)}")
    return hint


def classify_system_pool(
    molecules, system_type="auto", *, coordinate_model="covalent-radii",
    bond_scale=1.15, bond_cutoff=None,
):
    """Classify an XYZ pool conservatively for feature selection.

    XYZ-only inference recognizes metal/noble-gas elemental clusters,
    disconnected molecular aggregates, and graph-distinct molecular isomers.
    Homonuclear main-group pools (for example carbon) remain unknown because
    they may represent either molecules or covalent atomic clusters.
    """
    molecules = list(molecules)
    hint = normalize_system_type(system_type)
    for molecule in molecules:
        coordinates = np.asarray(molecule.coordinates, dtype=float)
        if (not len(molecule.atoms_list)
                or coordinates.shape != (len(molecule.atoms_list), 3)
                or not np.isfinite(coordinates).all()):
            raise ValueError("System classification requires nonempty structures with finite (natoms, 3) coordinates")

    def make_graphs():
        return [infer_coordinate_graph(
            molecule, model=coordinate_model, scale=bond_scale, cutoff=bond_cutoff
        ) for molecule in molecules]

    if hint != "auto":
        reason = "system type supplied by the caller"
        groups = ()
        if hint in {"isomers", "molecular-aggregate"}:
            try:
                groups = _topology_groups(molecules, make_graphs())
            except KeyError as exc:
                reason += f"; coordinate adjacency unavailable for element {exc.args[0]!r}"
        return SystemClassification(
            hint, "explicit", reason, groups,
        )
    if not molecules:
        return SystemClassification("unknown", "low", "empty structure pool")

    compositions = [Counter(str(atom).capitalize() for atom in molecule.atoms_list) for molecule in molecules]
    symbols = set().union(*(composition.keys() for composition in compositions))
    same_composition = all(composition == compositions[0] for composition in compositions[1:])
    if symbols <= _METAL_AND_NOBLE_CLUSTER_ELEMENTS and all(
        len(molecule.atoms_list) > 1 for molecule in molecules
    ):
        return SystemClassification(
            "atomic-cluster", "medium",
            "pool contains only metal/noble-gas atoms; no molecular bond-order model was applied",
        )
    if len(symbols) == 1 and same_composition:
        return SystemClassification(
            "unknown", "low",
            "homonuclear XYZ structures can be molecules or covalent clusters; specify --system-type",
        )

    try:
        graphs = make_graphs()
    except KeyError as exc:
        # Lack of a covalent radius must not disable a coordinate histogram.
        return SystemClassification(
            "unknown", "low", f"coordinate adjacency unavailable for element {exc.args[0]!r}",
        )
    component_counts = [nx.number_connected_components(graph) for graph in graphs]
    if all(count > 1 for count in component_counts) and all(
        any(len(component) > 1 for component in nx.connected_components(graph))
        for graph in graphs
    ):
        return SystemClassification(
            "molecular-aggregate", "medium",
            f"all structures contain multiple components including a multi-atom fragment under {coordinate_model} adjacency",
            _topology_groups(molecules, graphs),
        )

    if same_composition and all(count == 1 for count in component_counts):
        topology_group_ids = _topology_groups(molecules, graphs)
        if len(set(topology_group_ids)) > 1:
            return SystemClassification(
                "isomers", "medium",
                "same-composition connected structures have different inferred element-labelled graphs",
                topology_group_ids,
            )
        return SystemClassification(
            "molecular", "medium",
            "same-composition connected structures share an inferred element-labelled graph",
            topology_group_ids,
        )

    return SystemClassification(
        "unknown", "low",
        "XYZ-only composition and adjacency do not determine a unique system class",
    )


def resolve_clustering_policy(system_type="unknown", feature="auto", algorithm="auto"):
    """Resolve automatic feature/clusterer choices and explain the policy."""
    normalized_type = normalize_system_type(system_type or "unknown")
    if normalized_type == "auto":
        normalized_type = "unknown"
    if normalized_type not in FEATURE_POLICY:
        raise ValueError(f"Unknown system type {system_type!r}")
    normalized_feature = str(feature or "auto").strip().lower().replace("_", "-")
    normalized_feature = {"histogram": "distance-histogram", "pair-distance": "distance-histogram"}.get(
        normalized_feature, normalized_feature
    )
    if normalized_feature not in {"auto", "mbtr", "soap", "distance-histogram"}:
        raise ValueError(f"Unknown feature {feature!r}")
    if normalized_feature == "auto":
        selected_feature = FEATURE_POLICY[normalized_type]
        feature_reason = f"automatic feature policy for {normalized_type}"
    else:
        selected_feature = normalized_feature
        feature_reason = "feature explicitly requested"
    normalized_algorithm = str(algorithm or "auto").strip().lower()
    if normalized_algorithm not in {"auto", "hybrid", "hdbscan", "agglomerative", "ward", "dbscan", "optics", "maxmin", "max-min", "max_min"}:
        raise ValueError(f"Unknown clustering algorithm {algorithm!r}")
    if normalized_algorithm == "auto":
        selected_algorithm = DEFAULT_CLUSTER_ALGORITHM
        algorithm_reason = "HDBSCAN with agglomerative operational fallback for unknown cluster counts"
    else:
        selected_algorithm = normalized_algorithm
        algorithm_reason = "cluster-label algorithm explicitly requested"
    return {
        "system_type": normalized_type,
        "feature": selected_feature,
        "feature_fallbacks": (
            FEATURE_FALLBACK_POLICY[normalized_type] if normalized_feature == "auto"
            else {"mbtr": ("soap", "distance-histogram"),
                  "soap": ("mbtr", "distance-histogram"),
                  "distance-histogram": ()}[selected_feature]
        ),
        "algorithm": selected_algorithm,
        "reason": f"{feature_reason}; {algorithm_reason}",
    }


__all__ = [
    "DEFAULT_CLUSTER_ALGORITHM",
    "FEATURE_FALLBACK_POLICY",
    "FEATURE_POLICY",
    "SYSTEM_TYPES",
    "SELECTION_POLICY_VERSION",
    "SystemClassification",
    "classify_system_pool",
    "normalize_system_type",
    "resolve_clustering_policy",
]
