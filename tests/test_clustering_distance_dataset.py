"""Independent audit and score checks for constructed distance benchmarks."""

import json
from pathlib import Path

import numpy as np
import pytest

from pyar.scripts.benchmark_distances import load_dataset
from pyar.selection.clusterers import cluster_molecules


DATASET = Path(__file__).parents[1] / "benchmarks" / "clustering_distances" / "dataset.json"


def test_dataset_has_verified_rigid_copies_and_fragment_geometries():
    manifest, frames, audit = load_dataset(DATASET)
    assert audit["geometry_count"] == len(frames) == 84
    assert audit["pool_count"] == 7
    assert audit["witness_count"] == 42
    assert len(audit["manifest_sha256"]) == len(audit["structures_sha256"]) == 64
    assert "not optimized PES basins" in manifest["label_interpretation"]


def test_audit_rejects_a_corrupted_rigid_transform(tmp_path):
    manifest = json.loads(DATASET.read_text())
    manifest["pools"][0]["rigid_transform_witnesses"][0]["translation"][0] += 1
    (tmp_path / "dataset.json").write_text(json.dumps(manifest))
    (tmp_path / "structures.xyz").write_bytes(DATASET.with_name("structures.xyz").read_bytes())
    with pytest.raises(ValueError, match="transformation witness"):
        load_dataset(tmp_path / "dataset.json")


@pytest.mark.parametrize("pool_name,metric", [("water_dimer_orientation", "fragment-rmsd"),
                                              ("water_dimer_beyond_cutoff", "fragment-rmsd"),
                                              ("water_ammonia_packing", "soap-rematch"),
                                              ("butane_torsions", "graph-rmsd")])
def test_structural_clustering_recovers_constructed_families(pool_name, metric):
    from sklearn.metrics import adjusted_rand_score

    manifest, frames, _ = load_dataset(DATASET)
    pool = next(pool for pool in manifest["pools"] if pool["id"] == pool_name)
    molecules = [frames[record["id"]] for record in pool["records"]]
    result = cluster_molecules(molecules, distance_metric=metric, algorithm="agglomerative",
                              maximum_number_of_clusters=2, system_type=pool["system_type"])
    assert result.distance_used == metric
    assert adjusted_rand_score([record["family"] for record in pool["records"]], result.labels) == 1
    assert np.isfinite(result.distance_matrix).all()
