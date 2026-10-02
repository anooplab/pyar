"""Reproduce optimized-conformer retention at current production cutoffs.

Run from the repository root with the project Python environment. Saved basin
labels are operational references, not independent proofs of distinct minima.
"""
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from pyar.structure_comparison import GraphFirstDeduplicationComparator
from pyar.structure_comparison.fragment_rmsd import FragmentRMSDComparator
from pyar.selection.basin_memory import _persist_basin_registry, _load_basin_registry

ROOT = Path(__file__).resolve().parents[3]


def read_geometry(path):
    lines = path.read_text().splitlines()
    rows = [line.split() for line in lines[2:2 + int(lines[0])]]
    return SimpleNamespace(atoms_list=[row[0] for row in rows],
                           coordinates=np.asarray([row[1:4] for row in rows], float),
                           charge=0, multiplicity=1)


def main():
    output = {"reference_caveat": "Operational graph RMSD labels at 0.35 A, not ground-truth minima", "pools": []}
    water = np.array([[0., 0., 0.], [.9572, 0., 0.], [-.239987, .927297, 0.]])
    seven_waters = SimpleNamespace(
        atoms_list=["O", "H", "H"] * 7,
        coordinates=np.concatenate([water + [3.2 * i, 0., 0.] for i in range(7)]),
    )
    fragment_result = FragmentRMSDComparator().compare(seven_waters, seven_waters)
    output["seven_water_identity_fragment_comparison"] = {
        "distance": fragment_result.distance, "metadata": dict(fragment_result.metadata),
        "purpose": "Mapping-budget probe only; generated chain is not an optimized physical benchmark",
    }
    base = SimpleNamespace(atoms_list=["O", "H", "H"], coordinates=water,
                           charge=0, multiplicity=1, name="original_water", energy=-1.)
    distinct = SimpleNamespace(atoms_list=["O", "H", "H"], coordinates=water * 1.1,
                               charge=0, multiplicity=1, name="distorted_water", energy=0.)
    translations = [SimpleNamespace(atoms_list=base.atoms_list, coordinates=water + [i, 0., 0.],
                                    charge=0, multiplicity=1, name=f"translated_water_{i}", energy=-1.)
                    for i in (1, 2, 3)]
    registry = str(Path(__file__).with_name("comparison_archive_fixture.json"))
    _persist_basin_registry(registry, [distinct, base], existing_entries=[], max_entries=3)
    for molecule in translations:
        _persist_basin_registry(registry, [molecule], max_entries=3)
    entries = _load_basin_registry(registry)
    output["basin_archive_rigid_duplicate_eviction"] = {
        "initial_names": [distinct.name, base.name],
        "final_names": [entry["name"] for entry in entries],
        "added_rmsds_to_original": [GraphFirstDeduplicationComparator(threshold=.001).compare(base, mol).distance for mol in translations],
        "distinct_initial_geometry_evicted": distinct.name not in [entry["name"] for entry in entries],
    }
    comparator = GraphFirstDeduplicationComparator(atom_mode="heavy")
    for name in ("CAMVES_I", "FGG55", "WG01"):
        directory = ROOT / "benchmarks/clustering_scientific/runs" / name
        reference = json.loads((directory / "reference_basins.json").read_text())
        records = sorted(reference["structures"], key=lambda row: (row["xtb_energy_hartree"], row["name"]))
        geometries = [read_geometry(directory / row["optimized_xyz"]) for row in records]
        cache = {}
        def distance(i, j):
            key = tuple(sorted((i, j)))
            if key not in cache:
                cache[key] = comparator.compare(geometries[i], geometries[j]).distance
            return cache[key]
        conditions = []
        for cutoff in (.35, .5, .75):
            kept = []
            cross_label_removals = []
            for i, record in enumerate(records):
                for j in kept:
                    rmsd = distance(i, j)
                    if rmsd is not None and rmsd < cutoff:
                        if record["basin_id"] != records[j]["basin_id"]:
                            cross_label_removals.append({"removed": record["name"], "kept": records[j]["name"], "distance": rmsd})
                        break
                else:
                    kept.append(i)
            conditions.append({"cutoff": cutoff, "retained_count": len(kept),
                               "retained_reference_labels": len({records[i]["basin_id"] for i in kept}),
                               "cross_label_removals": cross_label_removals})
        output["pools"].append({"name": name, "input_count": len(records),
                                "reference_labels": len({r["basin_id"] for r in records}), "conditions": conditions})
        print(name, [(c["cutoff"], c["retained_count"], c["retained_reference_labels"], len(c["cross_label_removals"])) for c in conditions], flush=True)
    target = Path(__file__).with_suffix(".json")
    target.write_text(json.dumps(output, indent=2) + "\n")


if __name__ == "__main__":
    main()
