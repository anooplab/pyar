"""Adapt generated and refined conformers to the shared comparison policy."""

from collections import Counter
from types import SimpleNamespace

import numpy as np

from pyar.structure_comparison import ComparisonResult, GraphFirstDeduplicationComparator


class ConformerComparisons:
    """Compare complete structures while scoring the requested atom subset.

    Unavailable geometries and comparison failures abstain. Counters describe
    actual comparison outcomes, including diversity comparisons without a
    threshold; an abstention never authorizes removal.
    """

    def __init__(self, atom_mode="heavy", *, charge=None, multiplicity=1):
        if atom_mode not in {"heavy", "all"}:
            raise ValueError("atom_mode must be 'heavy' or 'all'")
        self.atom_mode = atom_mode
        self.charge = charge
        self.multiplicity = multiplicity
        self._comparators = {None: GraphFirstDeduplicationComparator(atom_mode=atom_mode)}
        self._counts = {}

    def _molecule(self, record):
        if record.molecule is not None:
            return record.molecule
        molecule = record.rdkit_molecule
        if molecule is None:
            raise ValueError("Conformer geometry is unavailable")
        conformer_id = (
            record.source_conf_id if record.rdkit_conf_id is None else record.rdkit_conf_id
        )
        conformer = molecule.GetConformer(int(conformer_id))
        coordinates = []
        symbols = []
        for atom in molecule.GetAtoms():
            position = conformer.GetAtomPosition(atom.GetIdx())
            symbols.append(atom.GetSymbol())
            coordinates.append([position.x, position.y, position.z])
        return SimpleNamespace(
            atoms_list=symbols, coordinates=np.asarray(coordinates, dtype=float),
            charge=(self.charge if self.charge is not None else
                    sum(atom.GetFormalCharge() for atom in molecule.GetAtoms())),
            multiplicity=self.multiplicity,
        )

    def compare(self, first, second, threshold=None, *, stage="deduplication"):
        """Return a shared policy result and record a compact outcome summary."""
        if threshold not in self._comparators:
            self._comparators[threshold] = GraphFirstDeduplicationComparator(
                threshold=threshold, atom_mode=self.atom_mode,
            )
        try:
            result = self._comparators[threshold].compare(
                self._molecule(first), self._molecule(second),
            )
        except Exception as exc:
            result = ComparisonResult(
                False, None, None, GraphFirstDeduplicationComparator.method, threshold,
                {"comparison_error": type(exc).__name__, "atom_mode": self.atom_mode},
            )
        counts = self._counts.setdefault(stage, Counter())
        counts["comparisons"] += 1
        if result.metadata.get("comparison_error"):
            counts["uncertain"] += 1
            counts["error_" + result.metadata["comparison_error"]] += 1
        elif not result.compatible:
            counts["incompatible"] += 1
        elif result.distance is None:
            counts["uncertain"] += 1
        elif result.equivalent is True:
            counts["equivalent"] += 1
        elif result.equivalent is False:
            counts["distinct"] += 1
        else:
            counts["distance_only"] += 1
        if result.metadata.get("fallback_method"):
            counts["irmsd_attempts"] += 1
            counts["irmsd_status_" + str(result.metadata.get("fallback_status", "unknown"))] += 1
            if result.metadata.get("fallback_mapping_verified"):
                counts["irmsd_verified"] += 1
        return result

    def summary(self):
        """Return JSON-compatible provenance for the comparisons performed."""
        comparator = self._comparators[None]
        return {
            "method": GraphFirstDeduplicationComparator.method,
            "policy_version": 1,
            "atom_mode": self.atom_mode,
            "bond_scale": comparator.bond_scale,
            "max_isomorphisms": comparator.graph_comparator.max_isomorphisms,
            "irmsd_timeout_seconds": comparator.irmsd_timeout,
            "uncertainty_policy": "keep",
            "stages": {stage: dict(counts) for stage, counts in self._counts.items()},
        }
