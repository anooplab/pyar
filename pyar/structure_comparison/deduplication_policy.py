"""Conservative graph-first policy for removing duplicate geometries."""

from __future__ import annotations

import json
import os
import subprocess
import sys

import numpy as np

from pyar.structure_comparison.graph_rmsd import GraphRMSDComparator
from pyar.structure_comparison.models import ComparisonResult


_IRMSD_BIDIRECTIONAL_WORKER = r'''
import json
import sys
import numpy as np
import irmsd

payload = json.load(sys.stdin)
first = irmsd.Molecule(symbols=payload["first_symbols"], positions=np.asarray(payload["first"], dtype=float))
second = irmsd.Molecule(symbols=payload["second_symbols"], positions=np.asarray(payload["second"], dtype=float))
for left, right in ((first, second), (second, first)):
    distance, _, _ = irmsd.get_irmsd_molecule(left, right, iinversion=2)
    print("PYAR_IRMSD_DISTANCE:" + json.dumps(float(distance), allow_nan=False), flush=True)
'''


def _run_bidirectional_irmsd(first, second, timeout):
    """Run native iRMSD in a child so warnings/fallbacks are observable."""
    payload = {
        "first_symbols": list(first.atoms_list),
        "second_symbols": list(second.atoms_list),
        "first": np.asarray(first.coordinates, dtype=float).tolist(),
        "second": np.asarray(second.coordinates, dtype=float).tolist(),
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
                 irmsd_timeout=30.0):
        self.threshold = threshold
        self.irmsd_timeout = float(irmsd_timeout)
        self.graph_comparator = GraphRMSDComparator(
            threshold=threshold, bond_scale=bond_scale,
            max_isomorphisms=max_isomorphisms,
        )

    def compare(self, first, second) -> ComparisonResult:
        graph_result = self.graph_comparator.compare(first, second)
        metadata = dict(graph_result.metadata)
        if (graph_result.metadata.get("comparison_complete") is not False
                or graph_result.metadata.get("connectivity_match") is not True):
            return graph_result

        distances, status = _run_bidirectional_irmsd(first, second, self.irmsd_timeout)
        metadata.update({"fallback_method": "irmsd", "fallback_status": status})
        if distances is None:
            return ComparisonResult(
                graph_result.compatible, None, None, self.method, self.threshold, metadata,
            )

        conservative_distance = max(distances)
        equivalent = None if self.threshold is None else conservative_distance < self.threshold
        metadata["irmsd_distances_angstrom"] = list(distances)
        return ComparisonResult(
            graph_result.compatible, conservative_distance, equivalent,
            self.method, self.threshold, metadata,
        )


__all__ = ["GraphFirstDeduplicationComparator"]
