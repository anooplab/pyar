"""Conservative graph-first policy for removing duplicate geometries."""

from __future__ import annotations

import json
import os
import subprocess
import sys
from collections import Counter

import numpy as np

from pyar.structure_comparison.graph_rmsd import GraphRMSDComparator
from pyar.structure_comparison.models import ComparisonResult
from pyar.structure_comparison.coordinate_graph import infer_coordinate_graph
from pyar.structure_comparison.rmsd import kabsch_rmsd


def _validated_irmsd_distance(first, second, matched_first, matched_second, distance, bond_scale):
    """Verify the native atom correspondence before trusting its distance."""
    from types import SimpleNamespace

    matched = [SimpleNamespace(atoms_list=list(molecule.symbols),
                               coordinates=np.asarray(molecule.positions, dtype=float))
               for molecule in (matched_first, matched_second)]
    for original, reordered in zip((first, second), matched):
        if Counter(original.atoms_list) != Counter(reordered.atoms_list):
            raise ValueError("iRMSD changed the element composition")
    if matched[0].atoms_list != matched[1].atoms_list:
        raise ValueError("iRMSD correspondence does not preserve element labels")
    graphs = [infer_coordinate_graph(molecule, scale=bond_scale) for molecule in matched]
    edges = [{tuple(sorted(edge)) for edge in graph.edges} for graph in graphs]
    if edges[0] != edges[1]:
        raise ValueError("iRMSD correspondence does not preserve coordinate adjacency")
    verified = kabsch_rmsd(matched[0].coordinates, matched[1].coordinates)
    if (not np.isfinite(distance) or distance < 0
            or not np.isclose(verified, distance, rtol=1e-5, atol=1e-7)):
        raise ValueError("iRMSD distance disagrees with the returned atom correspondence")
    return max(float(distance), verified)


_IRMSD_BIDIRECTIONAL_WORKER = r'''
import json
import sys
import numpy as np
import irmsd
from types import SimpleNamespace
from pyar.structure_comparison.deduplication_policy import _validated_irmsd_distance
from pyar.structure_comparison.graph_rmsd import selected_atom_indices
from pyar.structure_comparison.rmsd import kabsch_rmsd

payload = json.load(sys.stdin)
first = irmsd.Molecule(symbols=payload["first_symbols"], positions=np.asarray(payload["first"], dtype=float))
second = irmsd.Molecule(symbols=payload["second_symbols"], positions=np.asarray(payload["second"], dtype=float))
for left, right in ((first, second), (second, first)):
    distance, matched_left, matched_right = irmsd.get_irmsd_molecule(left, right, iinversion=2)
    distance = _validated_irmsd_distance(
        SimpleNamespace(atoms_list=left.symbols, coordinates=left.positions),
        SimpleNamespace(atoms_list=right.symbols, coordinates=right.positions),
        matched_left, matched_right, float(distance), payload["bond_scale"],
    )
    if payload.get("atom_mode", "all") == "heavy":
        indices = selected_atom_indices(matched_left.symbols, "heavy")
        distance = kabsch_rmsd(np.asarray(matched_left.positions)[indices],
                               np.asarray(matched_right.positions)[indices])
    print("PYAR_IRMSD_DISTANCE:" + json.dumps(float(distance), allow_nan=False), flush=True)
'''


def _run_bidirectional_irmsd(first, second, timeout, bond_scale=1.15, atom_mode="all"):
    """Run native iRMSD in a child so warnings/fallbacks are observable."""
    payload = {
        "first_symbols": list(first.atoms_list),
        "second_symbols": list(second.atoms_list),
        "first": np.asarray(first.coordinates, dtype=float).tolist(),
        "second": np.asarray(second.coordinates, dtype=float).tolist(),
        "bond_scale": bond_scale,
        "atom_mode": atom_mode,
    }
    environment = dict(os.environ)
    environment.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    try:
        completed = subprocess.run(
            [sys.executable, "-c", _IRMSD_BIDIRECTIONAL_WORKER],
            input=json.dumps(payload, allow_nan=False), text=True,
            capture_output=True, timeout=timeout, env=environment, check=False,
        )
    except subprocess.TimeoutExpired:
        return None, "timeout"
    except (OSError, TypeError, ValueError):
        return None, "launch_or_input_error"

    prefix = "PYAR_IRMSD_DISTANCE:"
    lines = completed.stdout.splitlines()
    values = [line[len(prefix):] for line in lines if line.startswith(prefix)]
    diagnostics = [line for line in lines if not line.startswith(prefix)]
    if completed.returncode != 0:
        return None, "backend_error"
    if diagnostics or completed.stderr.strip() or len(values) != 2:
        return None, "backend_diagnostic_or_incomplete_output"
    try:
        distances = tuple(float(value) for value in values)
    except (TypeError, ValueError):
        return None, "invalid_output"
    if any(not np.isfinite(value) or value < 0.0 for value in distances):
        return None, "invalid_distance"
    return distances, "ok"


class GraphFirstDeduplicationComparator:
    """Use graph RMSD first, with a strict bidirectional iRMSD second check.

    iRMSD is consulted only after graph isomorphism was established but the
    graph's mapping enumeration exceeded its limit. Both iRMSD argument orders
    must complete without native output or errors, and both distances must pass
    the threshold. Every other uncertain outcome remains an abstention.
    """

    method = "graph-rmsd-with-bidirectional-irmsd-on-incomplete-mapping"

    def __init__(self, threshold=None, bond_scale=1.15, max_isomorphisms=10000,
                 irmsd_timeout=30.0, atom_mode="all"):
        self.threshold = threshold
        self.irmsd_timeout = float(irmsd_timeout)
        if not np.isfinite(self.irmsd_timeout) or self.irmsd_timeout <= 0:
            raise ValueError("irmsd_timeout must be finite and positive")
        self.bond_scale = float(bond_scale)
        self.atom_mode = atom_mode
        self.graph_comparator = GraphRMSDComparator(
            threshold=threshold, bond_scale=bond_scale,
            max_isomorphisms=max_isomorphisms,
            atom_mode=atom_mode,
        )

    def compare(self, first, second) -> ComparisonResult:
        graph_result = self.graph_comparator.compare(first, second)
        metadata = dict(graph_result.metadata)
        if (graph_result.metadata.get("comparison_complete") is not False
                or graph_result.metadata.get("connectivity_match") is not True):
            return graph_result

        kwargs = {"atom_mode": self.atom_mode} if self.atom_mode != "all" else {}
        distances, status = _run_bidirectional_irmsd(
            first, second, self.irmsd_timeout, self.bond_scale, **kwargs,
        )
        metadata.update({"fallback_method": "irmsd", "fallback_status": status})
        if distances is None:
            return ComparisonResult(
                graph_result.compatible, None, None, self.method, self.threshold, metadata,
            )

        conservative_distance = max(distances)
        equivalent = None if self.threshold is None else conservative_distance < self.threshold
        metadata["irmsd_distances_angstrom"] = list(distances)
        metadata["fallback_mapping_verified"] = True
        return ComparisonResult(
            graph_result.compatible, conservative_distance, equivalent,
            self.method, self.threshold, metadata,
        )


__all__ = ["GraphFirstDeduplicationComparator"]
