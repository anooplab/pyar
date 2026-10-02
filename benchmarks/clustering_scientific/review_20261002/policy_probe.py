"""Independent metric checks and focused clustering/archive review probes.
Run from repository root with .venv/bin/python on this file.
Production sources and original benchmark artifacts are read-only.
"""
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from pyar.selection.basin_memory import _persist_basin_registry, _load_basin_registry
from pyar.selection.clusterers import cluster_molecules
from pyar.structure_comparison import GraphRMSDComparator

HERE = Path(__file__).resolve().parent
BENCH = HERE.parent


def manual_auc(truth, distances):
    positive = distances[truth]
    negative = distances[~truth]
    return float(np.mean((positive[:, None] < negative[None, :]).astype(float)
                         + .5 * (positive[:, None] == negative[None, :])))


def main():
    evidence = {}
    archived = []
    for name, length, shift in (("distinct_B", 1.4, 0), ("A", .74, 0),
                                ("translated_A_1", .74, 10), ("translated_A_2", .74, 20)):
        archived.append(SimpleNamespace(name=name, atoms_list=["H", "H"],
            coordinates=np.array([[shift, 0., 0.], [shift + length, 0., 0.]]),
            energy=-1., charge=0, multiplicity=1))
    archive = HERE / "archive_probe.json"
    _persist_basin_registry(archive, archived, existing_entries=[], max_entries=3)
    stored = _load_basin_registry(archive)
    comparator = GraphRMSDComparator(threshold=1e-8)
    evidence["archive_capacity"] = {
        "configured_capacity": 3,
        "stored_names": [row["name"] for row in stored],
        "translated_copies_verified_equivalent": all(
            comparator.compare(archived[1], mol).equivalent is True for mol in archived[2:]),
        "distinct_geometry_evicted": not any(row["name"] == "distinct_B" for row in stored),
        "fixture_note": "Synthetic H2 coordinates test storage invariance; no physical basin labels claimed",
    }

    pool=[]
    for i, smiles in enumerate(("CCO", "COC")):
        molecule=Chem.AddHs(Chem.MolFromSmiles(smiles))
        assert AllChem.EmbedMolecule(molecule, randomSeed=42) == 0
        pool.append(SimpleNamespace(name=smiles, atoms_list=[a.GetSymbol() for a in molecule.GetAtoms()],
            coordinates=np.array(molecule.GetConformer().GetPositions()), energy=float(i)))
    result=cluster_molecules(pool, feature="distance-histogram", algorithm="agglomerative",
                             distance_metric="graph-rmsd", maximum_number_of_clusters=2)
    evidence["mixed_topology_graph_rmsd"] = result.to_dict()

    pool=[SimpleNamespace(name=str(i), atoms_list=["H", "H"],
          coordinates=np.array([[0., 0., 0.], [.7+.3*i, 0., 0.]]), energy=float(i))
          for i in range(4)]
    evidence["explicit_dbscan_fallback"] = cluster_molecules(
        pool, feature="distance-histogram", algorithm="dbscan",
        algorithm_options={"eps":1e-12, "min_samples":2},
        maximum_number_of_clusters=2).to_dict()

    summary=json.loads((BENCH / "meta_analysis/meta_analysis.json").read_text())
    checks=[]
    for name in ("CAMVES_I", "FGG55", "WG01"):
        report=json.loads((BENCH / "runs" / name / "results/comparison.json").read_text())
        refs=json.loads((BENCH / "runs" / name / "reference_basins.json").read_text())
        labels={r["name"]:r["basin_id"] for r in refs["structures"]}
        basin=np.array([labels[f"frame_{i:04d}"] for i in range(len(labels))])
        left,right=np.triu_indices(len(basin),1)
        truth=basin[left] == basin[right]
        for row in report["conditions"]:
            cluster=np.array(row["diagnostics"]["labels"])
            pred=(cluster[left] == cluster[right]) & (cluster[left]>=0) & (cluster[right]>=0)
            counts={"tp":int(sum(pred & truth)), "fp":int(sum(pred & ~truth)),
                    "fn":int(sum(~pred & truth)), "tn":int(sum(~pred & ~truth))}
            algorithm="auto" if row["algorithm_requested"]=="hybrid" else row["algorithm_requested"]
            saved=next(r for r in summary["conformer_pair_confusion"]
                       if r["dataset"]==name and r["algorithm_requested"]==algorithm
                       and r["feature_requested"]==row["feature_requested"])
            assert all(saved[key] == value for key,value in counts.items())
            checks.append({"dataset":name,"algorithm":algorithm,"feature":row["feature_requested"],**counts})
    evidence["conformer_confusion_independent_check"]={"conditions_verified":len(checks),"results":checks}
    water_root=BENCH / "runs/water_similarity"
    reference=np.loadtxt(water_root / "reference_fragment_rmsd_upper_bound_angstrom.csv",delimiter=",")
    pair=np.triu_indices(len(reference),1)
    water=[]
    for row in summary["water_proxy_roc_auc_confusion"]:
        distances=np.loadtxt(water_root / f"{row['feature']}_euclidean_feature_distances.csv",delimiter=",")[pair]
        truth=reference[pair] <= row["reference_rmsd_proxy_cutoff_angstrom"]
        auc=manual_auc(truth,distances)
        assert np.isclose(auc,row["proxy_pair_auc"])
        if row["specificity_target"] is None:
            continue
        pred=np.zeros_like(truth) if row["feature_distance_threshold"] is None else distances<=row["feature_distance_threshold"]
        tp,fp,fn,tn=(int(sum(pred & truth)),int(sum(pred & ~truth)),int(sum(~pred & truth)),int(sum(~pred & ~truth)))
        assert (tp,fp,fn,tn)==tuple(row[key] for key in ("tp","fp","fn","tn"))
        water.append({"feature":row["feature"],"specificity_target":row["specificity_target"],
                      "tp":tp,"fp":fp,"fn":fn,"tn":tn,"auc":auc,
                      "specificity":tn/(tn+fp),"precision":tp/(tp+fp) if tp+fp else None,
                      "false_positive_rate":fp/(fp+tn),
                      "false_discovery_fraction":fp/(fp+tp) if fp+tp else None})
    evidence["water_independent_metrics"]=water
    Path(__file__).with_suffix(".json").write_text(json.dumps(evidence,indent=2)+"\n")
    print(json.dumps({key: value for key,value in evidence.items() if key!="conformer_confusion_independent_check"},indent=2))


if __name__=="__main__":
    main()
