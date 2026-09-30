from types import SimpleNamespace
from unittest import mock

import numpy as np
import pytest

from pyar.selection.clusterers import (
    _cluster_agglomerative,
    cluster_molecules,
    determine_dbscan_params,
)


def _molecules(n=4):
    return [
        SimpleNamespace(
            name=f"m{i}", atoms_list=["H", "H"],
            coordinates=np.array([[0.0, 0.0, 0.0], [0.5 + i, 0.0, 0.0]]),
        )
        for i in range(n)
    ]


def test_agglomerative_returns_requested_number_of_groups_for_tied_distances():
    values = np.array([[0.0], [1.0], [10.0], [11.0]])
    labels = _cluster_agglomerative(values, 2)
    assert len(set(labels)) == 2


def test_clusterer_failure_falls_back_with_provenance():
    with mock.patch("pyar.selection.clusterers._run_algorithm", side_effect=[RuntimeError("missing HDBSCAN"), np.array([0, 0, 1, 1])]):
        result = cluster_molecules(_molecules(), feature="distance-histogram", algorithm="hybrid")
    assert result.algorithm_used == "agglomerative"
    assert result.algorithm_fallbacks[0]["algorithm"] == "hdbscan"
    assert "missing HDBSCAN" in result.algorithm_fallbacks[0]["reason"]
    assert result.number_of_clusters == 2


def test_all_noise_result_uses_fallback_clusterer():
    with mock.patch("pyar.selection.clusterers._run_algorithm", side_effect=[np.full(4, -1), np.array([0, 0, 1, 1])]):
        result = cluster_molecules(_molecules(), feature="distance-histogram", algorithm="dbscan")
    assert result.algorithm_used == "agglomerative"
    assert result.number_of_noise_points == 0


def test_bad_label_shape_is_rejected_and_fallback_is_used():
    with mock.patch("pyar.selection.clusterers._run_algorithm", side_effect=[np.array([0]), np.array([0, 0, 1, 1])]):
        result = cluster_molecules(_molecules(), feature="distance-histogram", algorithm="optics")
    assert result.algorithm_used == "agglomerative"
    assert "shape" in result.algorithm_fallbacks[0]["reason"]


def test_dbscan_epsilon_is_euclidean_not_squared_distance():
    epsilon, min_samples = determine_dbscan_params(np.array([[0.0], [3.0], [7.0]]))
    assert epsilon == pytest.approx(3.5)
    assert min_samples == 2


def test_unknown_algorithm_is_not_silently_changed():
    with pytest.raises(ValueError, match="Unknown clustering algorithm"):
        cluster_molecules(_molecules(), algorithm="not-a-clusterer")


def test_total_clusterer_failure_preserves_each_structure_as_a_cluster():
    with mock.patch("pyar.selection.clusterers._run_algorithm", side_effect=RuntimeError("backend failure")):
        result = cluster_molecules(
            _molecules(), feature="distance-histogram", algorithm="dbscan",
        )
    assert result.algorithm_used == "singleton-preservation"
    assert result.labels.tolist() == [0, 1, 2, 3]
    assert len(result.algorithm_fallbacks) == 3
