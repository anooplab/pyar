"""Regressions from the September 30 clustering and deduplication review."""

from types import SimpleNamespace
from unittest import mock

import numpy as np
import pytest

from pyar.selection import clustering
from pyar.selection.clusterers import _run_algorithm, cluster_molecules
from pyar.selection.policy import classify_system_pool


def _pool():
    return [SimpleNamespace(name=f"m{i}", atoms_list=["H", "H"], energy=float(i),
                            coordinates=np.array([[0., 0., 0.], [distance, 0., 0.]]))
            for i, distance in enumerate((0.4, 0.8, 1.4, 2.0))]


@pytest.mark.parametrize("metric", ["euclidean", "manhattan", "cosine"])
def test_real_hdbscan_supports_every_advertised_metric(metric):
    pytest.importorskip("hdbscan")
    values = np.array([[1, 0], [1, .01], [1, -.01], [-1, 0], [-1, .01], [-1, -.01]])
    labels = _run_algorithm(values, "hdbscan", 2, {"metric": metric, "min_samples": 1})
    assert labels[0] == labels[1] == labels[2] >= 0
    assert labels[3] == labels[4] == labels[5] >= 0
    assert labels[0] != labels[3]


@pytest.mark.parametrize("size", [1, 2])
def test_strict_connectivity_is_enforced_even_when_pool_fits_budget(size):
    with mock.patch("pyar.sampling.trial_generator.broken", return_value=True):
        assert clustering.choose_geometries(
            _pool()[:size], maximum_number_of_seeds=4, connectivity_policy="strict",
        ) == []


@pytest.mark.parametrize("pool", [[], _pool()[:1], _pool()])
def test_zero_seed_budget_is_rejected_for_every_pool_size(pool):
    with pytest.raises(ValueError, match="maximum_number_of_seeds"):
        clustering.choose_geometries(pool, maximum_number_of_seeds=0)


@pytest.mark.parametrize("options", [{"eps": -1}, {"xi": float("nan")}, {"min_samples": 0}])
def test_invalid_clustering_options_cannot_trigger_a_silent_fallback(options):
    with pytest.raises(ValueError):
        clustering.choose_geometries(_pool()[:1], algorithm_options=options)


def test_collapsed_histogram_preserves_every_candidate_with_provenance():
    pool = _pool()[:2]
    pool[1].coordinates[1, 0] = .41  # Both distances lie in the same histogram bin.
    with mock.patch("pyar.selection.clusterers._run_algorithm") as runner:
        result = cluster_molecules(pool, feature="distance-histogram")
    runner.assert_not_called()
    assert result.algorithm_used == "singleton-preservation"
    assert result.labels.tolist() == [0, 1]
    assert "no variation" in result.algorithm_fallbacks[0]["reason"]


def test_benchmark_reports_the_actual_pruned_selection_run():
    from pyar.scripts.benchmark_clustering import _benchmark_algorithm

    pool = _pool()
    with mock.patch("pyar.selection.clustering.remove_similar", return_value=pool[1:]), mock.patch(
        "pyar.selection.clusterers._run_algorithm", return_value=np.array([0, 0, 1])
    ) as runner:
        row = _benchmark_algorithm(pool, "agglomerative", 1, "distance-histogram")
    runner.assert_called_once()
    assert row["algorithm_used"] == "agglomerative"
    assert row["selected_names"] == ["m1"]


def test_mixed_metal_cluster_is_not_classified_as_a_molecule():
    molecule = SimpleNamespace(atoms_list=["Au", "Ag"], coordinates=[[0, 0, 0], [2.5, 0, 0]])
    assert classify_system_pool([molecule]).system_type == "atomic-cluster"


def test_disconnected_atoms_are_not_molecular_aggregate_evidence():
    molecule = SimpleNamespace(atoms_list=["C", "H"], coordinates=[[0, 0, 0], [5, 0, 0]])
    assert classify_system_pool([molecule]).system_type == "unknown"


def test_aggregate_pools_keep_different_constituent_topologies_separate():
    atoms = ["O", "O", "H", "H"]
    water_and_oxygen = SimpleNamespace(atoms_list=atoms,
                                      coordinates=[[0, 0, 0], [5, 0, 0], [.95, 0, 0], [-.95, 0, 0]])
    two_hydroxyls = SimpleNamespace(atoms_list=atoms,
                                   coordinates=[[0, 0, 0], [5, 0, 0], [.95, 0, 0], [5.95, 0, 0]])
    result = cluster_molecules([water_and_oxygen, two_hydroxyls], feature="distance-histogram")
    assert result.system_type == "molecular-aggregate"
    assert result.topology_group_ids == (0, 1)
    assert result.labels.tolist() == [0, 1]


def test_classification_uses_the_requested_coordinate_model():
    molecule = SimpleNamespace(atoms_list=["C", "H"], coordinates=[[0, 0, 0], [.9, 0, 0]])
    assert classify_system_pool([molecule]).system_type == "molecular"
    assert classify_system_pool([molecule], coordinate_model="none").system_type == "unknown"


def test_comparison_failure_retains_candidates():
    with mock.patch("pyar.structure_comparison.GraphFirstDeduplicationComparator") as comparator:
        comparator.return_value.compare.side_effect = RuntimeError("comparison unavailable")
        result = clustering.remove_similar(_pool()[:2])
    assert [molecule.name for molecule in result] == ["m0", "m1"]


def test_missing_covalent_radius_does_not_disable_the_histogram_fallback():
    from pyar.selection.features import compute_feature_matrix

    molecule = SimpleNamespace(name="unknown-radius", atoms_list=["Og", "H"],
                               coordinates=np.array([[0, 0, 0], [2, 0, 0]]))
    with mock.patch("pyar.selection.policy.infer_coordinate_graph", side_effect=KeyError("Og")):
        result = compute_feature_matrix([molecule], "distance-histogram")
    assert result.system_type == "unknown"
    assert result.values.shape[0] == 1


def test_cluster_cli_report_tracks_the_selected_geometries(tmp_path, monkeypatch, capsys):
    import json
    from pyar.scripts import clustering as script

    paths = []
    for i, distance in enumerate((.4, .8, 1.4, 2.0)):
        path = tmp_path / f"m{i}.xyz"
        path.write_text(f"2\nm{i}: energy={i}\nH 0 0 0\nH {distance} 0 0\n")
        paths.append(str(path))
    output = tmp_path / "report.json"
    monkeypatch.setattr("sys.argv", ["pyar-clustering", *paths, "-n", "2",
                                    "-a", "agglomerative", "--feature", "distance-histogram",
                                    "--report-output", str(output)])
    script.main()
    report = json.loads(output.read_text())
    assert report["input_files"] == paths
    assert len(report["labels"]) == len(report["candidate_names"])
    assert report["selected_count"] == len(report["selected_names"]) <= 2
    assert report["algorithm_used"] == "agglomerative"


def test_workflow_positional_distance_argument_retains_its_original_slot():
    import inspect
    from pyar.workflows._growth import add_one

    bound = inspect.signature(add_one).bind(
        "aggregate", [], None, 8, {}, 4, None, "off", "soap", "agglomerative", "cosine",
    )
    assert bound.arguments["selection_distance"] == "cosine"
