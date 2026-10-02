"""Coordinate-only structural descriptions for XYZ geometries.

This module intentionally infers adjacency only. XYZ coordinates do not encode
bond orders, charge, spin, or whether a short contact is chemically bonded.
The result is therefore a geometric graph plus diagnostics, not a chemical
identity assignment.
"""

from __future__ import annotations

from collections import Counter

import networkx as nx
import numpy as np

from pyar.data import new_atomic_data


def infer_coordinate_graph(molecule, *, model="covalent-radii", scale=1.15, cutoff=None):
    """Infer an element-labelled adjacency graph from XYZ coordinates only.

    ``covalent-radii`` connects atom pairs within ``scale`` times the sum of
    their covalent radii. ``distance-cutoff`` uses an explicit Angstrom cutoff.
    ``none`` returns only labelled vertices. Edge distances are retained; no
    bond orders or charges are assigned.
    """
    model = str(model).strip().lower().replace("_", "-")
    aliases = {"covalent": "covalent-radii", "radius": "covalent-radii"}
    model = aliases.get(model, model)
    if model not in {"covalent-radii", "distance-cutoff", "none"}:
        raise ValueError("model must be 'covalent-radii', 'distance-cutoff', or 'none'")
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError("scale must be a finite positive number")
    if model == "distance-cutoff" and (
        cutoff is None or not np.isfinite(cutoff) or cutoff <= 0
    ):
        raise ValueError("distance-cutoff model requires a finite positive cutoff in Angstrom")

    symbols = [str(atom).capitalize() for atom in molecule.atoms_list]
    coordinates = np.asarray(molecule.coordinates, dtype=float)
    if coordinates.shape != (len(symbols), 3) or not np.isfinite(coordinates).all():
        raise ValueError("molecule coordinates must be a finite (natoms, 3) array")
    graph = nx.Graph()
    graph.add_nodes_from((index, {"element": element}) for index, element in enumerate(symbols))
    if model in {"covalent-radii", "distance-cutoff"}:
        radii = (
            [float(new_atomic_data.covalent_radius[element]) for element in symbols]
            if model == "covalent-radii"
            else None
        )
        for left in range(len(symbols)):
            for right in range(left + 1, len(symbols)):
                distance = float(np.linalg.norm(coordinates[left] - coordinates[right]))
                radius_sum = radii[left] + radii[right] if radii is not None else None
                limit = (
                    float(scale) * radius_sum
                    if model == "covalent-radii"
                    else float(cutoff)
                )
                if distance <= limit:
                    edge_data = {"distance": distance, "normalized_distance": distance / limit}
                    if radius_sum is not None:
                        edge_data["radius_normalized_distance"] = distance / radius_sum
                    graph.add_edge(left, right, **edge_data)
    return graph


def analyze_coordinate_structure(molecule, *, model="covalent-radii", scale=1.15, cutoff=None):
    """Return a JSON-ready, charge- and bond-order-free structure summary."""
    graph = infer_coordinate_graph(molecule, model=model, scale=scale, cutoff=cutoff)
    components = [sorted(int(index) for index in component) for component in nx.connected_components(graph)]
    elements = Counter(data["element"] for _, data in graph.nodes(data=True))
    component_count = len(components)
    return {
        "name": str(getattr(molecule, "name", "<unnamed>")),
        "method": "coordinate-only-adjacency",
        "model": str(model),
        "scale": float(scale),
        "cutoff_angstrom": None if cutoff is None else float(cutoff),
        "natoms": graph.number_of_nodes(),
        "composition": dict(sorted(elements.items())),
        "edges": [
            {
                "atom_indices": [int(left), int(right)],
                "elements": [graph.nodes[left]["element"], graph.nodes[right]["element"]],
                "distance": float(data["distance"]),
                "normalized_distance": float(data["normalized_distance"]),
            }
            for left, right, data in graph.edges(data=True)
        ],
        "degrees": [int(graph.degree[index]) for index in graph.nodes],
        "components": components,
        "component_count": component_count,
        "geometric_component_pattern": (
            "single-connected-component" if component_count == 1
            else "multiple-components" if component_count > 1
            else "empty"
        ),
        "declared_fragment_count": len(getattr(molecule, "fragments", ()) or ()),
        "interpretation": "geometric adjacency only; not a bond-order or chemical-identity assignment",
    }


def compare_coordinate_structures(
    before, after, *, model="covalent-radii", scale=1.15, cutoff=None
):
    """Compare coordinate-derived summaries without assuming atom mapping.

    The summary reports composition, component count, and graph edge counts.
    Edge identities are only compared when atom ordering is identical, which
    is the case for PyAR growth geometries assembled from their input fragments.
    """
    before_graph = infer_coordinate_graph(before, model=model, scale=scale, cutoff=cutoff)
    after_graph = infer_coordinate_graph(after, model=model, scale=scale, cutoff=cutoff)
    same_order = [str(atom).capitalize() for atom in before.atoms_list] == [
        str(atom).capitalize() for atom in after.atoms_list
    ]
    before_edges = {tuple(sorted(edge)) for edge in before_graph.edges}
    after_edges = {tuple(sorted(edge)) for edge in after_graph.edges}
    result = {
        "same_atom_count": before_graph.number_of_nodes() == after_graph.number_of_nodes(),
        "same_composition": Counter(data["element"] for _, data in before_graph.nodes(data=True))
        == Counter(data["element"] for _, data in after_graph.nodes(data=True)),
        "same_atom_order": same_order,
        "input_edge_count": len(before_edges),
        "output_edge_count": len(after_edges),
        "input_component_count": nx.number_connected_components(before_graph),
        "output_component_count": nx.number_connected_components(after_graph),
    }
    if same_order:
        result["added_edges"] = [list(edge) for edge in sorted(after_edges - before_edges)]
        result["removed_edges"] = [list(edge) for edge in sorted(before_edges - after_edges)]
    else:
        result["edge_changes"] = "not calculated because atom order differs"
    return result


def analyze_growth_transition(
    input_parts, output, *, model="covalent-radii", scale=1.15, cutoff=None
):
    """Summarize a growth step as input fragment graphs versus product graph."""
    input_parts = list(input_parts)
    input_reports = [
        analyze_coordinate_structure(part, model=model, scale=scale, cutoff=cutoff)
        for part in input_parts
    ]
    output_report = analyze_coordinate_structure(
        output, model=model, scale=scale, cutoff=cutoff
    )
    input_symbols = [
        str(atom).capitalize() for part in input_parts for atom in part.atoms_list
    ]
    output_symbols = [str(atom).capitalize() for atom in output.atoms_list]
    input_edges = set()
    offset = 0
    input_components = 0
    for report in input_reports:
        input_edges.update(
            tuple(sorted((left + offset, right + offset)))
            for left, right in (edge["atom_indices"] for edge in report["edges"])
        )
        input_components += report["component_count"]
        offset += report["natoms"]
    output_edges = {
        tuple(sorted(edge["atom_indices"])) for edge in output_report["edges"]
    }
    same_order = input_symbols == output_symbols
    comparison = {
        "same_atom_count": len(input_symbols) == len(output_symbols),
        "same_composition": Counter(input_symbols) == Counter(output_symbols),
        "same_atom_order": same_order,
        "input_component_count": input_components,
        "output_component_count": output_report["component_count"],
        "input_fragment_count": len(input_parts),
        "input_edge_count": len(input_edges),
        "output_edge_count": len(output_edges),
    }
    if same_order:
        comparison["added_edges"] = [list(edge) for edge in sorted(output_edges - input_edges)]
        comparison["removed_edges"] = [list(edge) for edge in sorted(input_edges - output_edges)]
    else:
        comparison["edge_changes"] = "not calculated because atom order differs"
    return {
        "method": "coordinate-only-adjacency",
        "model": model,
        "scale": float(scale),
        "cutoff_angstrom": None if cutoff is None else float(cutoff),
        "input_parts": input_reports,
        "output": output_report,
        "input_to_output": comparison,
        "limitations": [
            "No bond orders, charges, spin states, or fragment intent are inferred.",
            "An edge change is a coordinate-threshold change, not proof of a chemical reaction.",
        ],
    }


__all__ = [
    "analyze_coordinate_structure", "analyze_growth_transition",
    "compare_coordinate_structures", "infer_coordinate_graph",
]
