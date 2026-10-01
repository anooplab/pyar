"""Structural distance invariants, fallback scales, and selection integration."""

import json
from types import SimpleNamespace
from unittest import mock

import numpy as np
import pytest

from pyar.selection.clusterers import ClusteringResult, _run_algorithm, cluster_molecules, determine_dbscan_params
from pyar.selection.clustering import choose_geometries
from pyar.selection.distances import pairwise_distances, validate_distance_matrix
from pyar.selection.structural_distances import (
    _rematch_similarity, compute_distance_matrix, validate_distance_options,
)
from pyar.structure_comparison.fragment_rmsd import FragmentRMSDComparator


def aggregate(separation=3.5, *, name="aggregate", energy=0, orientation=0):
    water = np.array([[0., 0., 0.], [.9572, 0., 0.], [-.239, .927, 0.]])
    angle = np.deg2rad(orientation)
    rotation = np.array([[np.cos(angle), -np.sin(angle), 0],
                         [np.sin(angle), np.cos(angle), 0], [0, 0, 1]])
    return SimpleNamespace(name=name, energy=energy, atoms_list=["O", "H", "H"] * 2,
                           coordinates=np.vstack((water, water @ rotation + [separation, 0, 0])))


def transformed(molecule):
    rng = np.random.default_rng(29)
    order = rng.permutation(len(molecule.atoms_list))
    rotation, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    rotation[:, -1] *= np.linalg.det(rotation)
    return SimpleNamespace(name="transformed", energy=1,
                           atoms_list=np.asarray(molecule.atoms_list)[order].tolist(),
                           coordinates=molecule.coordinates[order] @ rotation + [6, -1, 3])


@pytest.mark.parametrize("metric", ["graph-rmsd", "fragment-rmsd", "soap-rematch"])
def test_aggregate_matrix_is_invariant_and_distinguishes_packing(metric):
    if metric == "soap-rematch":
        pytest.importorskip("dscribe")
    original = aggregate()
    molecules = [original, transformed(original), aggregate(4.5, name="expanded")]
    snapshots = [molecule.coordinates.copy() for molecule in molecules]
    result = compute_distance_matrix(molecules, metric, allow_fallbacks=False)
    assert result.used == metric
    validate_distance_matrix(result.values, 3)
    assert result.values[0, 1] < 1e-6
    assert result.values[0, 2] > .01
    np.testing.assert_allclose(result.values[0, 2], result.values[1, 2], atol=1e-7)
    json.dumps(result.to_dict(), allow_nan=False)
    for molecule, snapshot in zip(molecules, snapshots):
        np.testing.assert_array_equal(molecule.coordinates, snapshot)


def test_fragment_matching_retains_relative_orientation_and_separation():
    result = compute_distance_matrix([aggregate(), aggregate(orientation=100)],
                                     "fragment-rmsd", allow_fallbacks=False)
    assert result.values[0, 1] > .1
    assert result.units == "angstrom"
    assert result.parameters["distance_is_upper_bound"] is True
    result = FragmentRMSDComparator().compare(aggregate(), aggregate(4.5))
    assert result.distance == pytest.approx(.5)


def test_fragment_matching_does_not_trust_stale_workflow_fragments():
    first = aggregate()
    second = aggregate()
    second.coordinates[3:] -= [2.9, 0, 0]  # Fragments now overlap/connect geometrically.
    first.fragments = second.fragments = [{"atoms": list(range(3))}, {"atoms": list(range(3, 6))}]
    result = FragmentRMSDComparator().compare(first, second)
    assert not result.compatible
    assert result.distance is None


def test_fragment_mapping_limit_abstains_explicitly():
    result = FragmentRMSDComparator(max_mappings=1).compare(aggregate(), aggregate())
    assert result.distance is None
    assert result.metadata["comparison_complete"] is False


def test_graph_matrix_requires_compatible_topology_without_fallbacks():
    first = SimpleNamespace(atoms_list=["C", "O"], coordinates=np.array([[0, 0, 0], [1.4, 0, 0]]))
    second = SimpleNamespace(atoms_list=["C", "O"], coordinates=np.array([[0, 0, 0], [4., 0, 0]]))
    with pytest.raises(RuntimeError, match="no graph-rmsd distance"):
        compute_distance_matrix([first, second], "graph-rmsd", allow_fallbacks=False)


def test_fallback_rebuilds_the_entire_matrix_in_one_scale():
    pool = [aggregate(), aggregate(4), aggregate(5)]
    features = np.array([[0], [2], [7]])
    with mock.patch("pyar.selection.structural_distances._mapped_matrix", side_effect=RuntimeError("mapping unavailable")), \
         mock.patch("pyar.selection.structural_distances._soap_rematch_matrix", side_effect=ImportError("DScribe missing")):
        result = compute_distance_matrix(pool, "fragment-rmsd", feature_values=features)
    assert result.used == "euclidean"
    assert result.units == "standardized-feature-units"
    assert [failure["distance"] for failure in result.fallbacks] == ["fragment-rmsd", "graph-rmsd", "soap-rematch"]
    np.testing.assert_array_equal(result.values, pairwise_distances(features))


def test_invalid_input_cannot_be_hidden_by_a_distance_fallback():
    molecule = aggregate()
    molecule.coordinates[0, 0] = float("nan")
    with pytest.raises(ValueError, match="finite"):
        compute_distance_matrix([molecule, aggregate()], "fragment-rmsd", feature_values=[[1], [2]])


@pytest.mark.parametrize("options", [{"atom_mode": "oxygen"}, {"max_mappings": 0},
                                     {"rematch_alpha": 0}, {"soap_cutoff": float("nan")},
                                     {"rematch_iterations": True}, {"unknown": 4}])
def test_invalid_scientific_options_fail_before_any_fallback(options):
    with pytest.raises(ValueError):
        validate_distance_options(options)


@pytest.mark.parametrize("matrix", [[[0, 1], [2, 0]], [[1]], [[0, float("nan")], [float("nan"), 0]],
                                    [[0, -1], [-1, 0]], [[0, 1, 2]]])
def test_invalid_precomputed_matrix_is_rejected(matrix):
    with pytest.raises(ValueError):
        validate_distance_matrix(matrix)


@pytest.mark.parametrize("algorithm", ["agglomerative", "hdbscan", "dbscan", "optics"])
def test_real_clusterers_consume_structural_matrix(algorithm):
    pytest.importorskip("sklearn")
    if algorithm == "hdbscan":
        pytest.importorskip("hdbscan")
    pool = [aggregate(s, name=f"m{i}") for i, s in enumerate((3.4, 3.5, 3.6, 6.4, 6.5, 6.6))]
    result = cluster_molecules(pool, distance_metric="fragment-rmsd", algorithm=algorithm,
                              maximum_number_of_clusters=2,
                              algorithm_options={"min_samples": 2, "min_cluster_size": 2, "eps": .2})
    assert result.distance_used == "fragment-rmsd"
    assert result.algorithm_used == algorithm
    assert result.labels[0] == result.labels[1] == result.labels[2] >= 0
    assert result.labels[3] == result.labels[4] == result.labels[5] >= 0
    assert result.labels[0] != result.labels[3]


def test_dbscan_estimates_epsilon_from_actual_precomputed_distances():
    distances = np.array([[0., 3., 7.], [3., 0., 4.], [7., 4., 0.]])
    assert determine_dbscan_params(distances, metric="precomputed")[0] == pytest.approx(3.5)


def test_graph_distances_do_not_require_descriptor_packages():
    with mock.patch("pyar.selection.features.compute_feature_matrix", side_effect=AssertionError("descriptor used")):
        result = cluster_molecules([aggregate(), aggregate(4)], distance_metric="graph-rmsd",
                                  algorithm="agglomerative", maximum_number_of_clusters=2)
    assert result.feature_used == "coordinates"
    assert result.feature_values.shape == (2, 0)


def test_descriptor_fallback_discards_explicit_radius_in_different_units():
    with mock.patch("pyar.selection.structural_distances._mapped_matrix", side_effect=RuntimeError("no mapping")), \
         mock.patch("pyar.selection.structural_distances._soap_rematch_matrix", side_effect=ImportError("no DScribe")), \
         mock.patch("pyar.selection.clusterers._run_algorithm", return_value=np.array([0, 0])) as runner:
        result = cluster_molecules([aggregate(), aggregate(4)], feature="distance-histogram",
                                  distance_metric="fragment-rmsd", algorithm="dbscan",
                                  algorithm_options={"eps": .1})
    assert result.distance_used == "euclidean"
    assert "eps" not in runner.call_args.args[3]
    assert result.distance_parameters["requested_eps_discarded"] == .1


def test_cluster_minima_trimming_reuses_structural_distances():
    pool = [aggregate(name=f"m{i}", energy=i) for i in range(4)]
    positions = np.array([[0.], [.1], [5.], [.2]])
    distance = pairwise_distances(positions)
    clustered = ClusteringResult(np.arange(4), "auto", "coordinates", "auto", "agglomerative",
                                 "fragment-rmsd", np.empty((4, 0)), distance_matrix=distance,
                                 distance_used="fragment-rmsd", distance_units="angstrom")
    with mock.patch("pyar.selection.clustering.remove_similar", return_value=pool), \
         mock.patch("pyar.selection.clustering.cluster_molecules", return_value=clustered):
        selected = choose_geometries(pool, maximum_number_of_seeds=2, distance_metric="fragment-rmsd",
                                    persist_basin_memory=False, apply_basin_memory=False)
    assert [molecule.name for molecule in selected] == ["m0", "m2"]


def test_graph_zero_distances_are_valid_even_without_descriptor_variation():
    pool = [aggregate(), transformed(aggregate())]
    result = cluster_molecules(pool, distance_metric="graph-rmsd", algorithm="agglomerative")
    assert result.algorithm_used == "agglomerative"
    assert result.labels.tolist() == [0, 0]


def test_local_soap_cutoff_collapse_preserves_candidates():
    pytest.importorskip("dscribe")
    result = cluster_molecules([aggregate(15), aggregate(20)], distance_metric="soap-rematch",
                              algorithm="agglomerative")
    assert result.distance_used == "soap-rematch"
    assert result.algorithm_used == "singleton-preservation"
    assert result.labels.tolist() == [0, 1]
    np.testing.assert_array_equal(result.distance_matrix, np.zeros((2, 2)))


def test_mean_soap_cutoff_collapse_does_not_turn_roundoff_into_distances():
    pytest.importorskip("dscribe")
    pool = [aggregate(15), transformed(aggregate(15)), aggregate(20)]
    result = cluster_molecules(pool, feature="soap", algorithm="agglomerative")
    assert result.algorithm_used == "singleton-preservation"


def test_explicit_aggregate_type_still_separates_constituent_topologies():
    first = SimpleNamespace(name="water-oxygen", atoms_list=["O", "O", "H", "H"],
                            coordinates=np.array([[0, 0, 0], [5, 0, 0], [.95, 0, 0], [-.95, 0, 0]]))
    second = SimpleNamespace(name="hydroxyls", atoms_list=first.atoms_list,
                             coordinates=np.array([[0, 0, 0], [5, 0, 0], [.95, 0, 0], [5.95, 0, 0]]))
    result = cluster_molecules([first, second], system_type="molecular-aggregate",
                              distance_metric="soap-rematch", algorithm="agglomerative")
    assert result.labels.tolist() == [0, 1]


@pytest.mark.parametrize("metric", ["graph-rmsd", "fragment-rmsd", "soap-rematch"])
def test_cli_labels_reports_structural_distance_and_parameters(tmp_path, monkeypatch, metric):
    from pyar.scripts import clustering as command

    pool = [aggregate(name="compact"), aggregate(4.5, name="expanded")]
    paths = []
    for molecule in pool:
        path = tmp_path / (molecule.name + ".xyz")
        lines = [str(len(molecule.atoms_list)), f"{molecule.name}: energy=0"]
        lines.extend(f"{atom} {x} {y} {z}" for atom, (x, y, z) in zip(molecule.atoms_list, molecule.coordinates))
        path.write_text("\n".join(lines) + "\n")
        paths.append(str(path))
    report = tmp_path / "report.json"
    monkeypatch.setattr("sys.argv", ["pyar-clustering", *paths, "--mode", "labels", "--distance", metric,
                                    "--distance-atom-mode", "all", "-a", "agglomerative",
                                    "--report-output", str(report)])
    command.main()
    payload = json.loads(report.read_text())
    assert payload["distance_requested"] == payload["distance_used"] == metric
    assert payload["distance_parameters"]
    assert len(payload["labels"]) == 2


def test_bounded_rematch_matches_dscribe_reference():
    kernels = pytest.importorskip("dscribe.kernels")
    local = np.array([[.9, .3], [.2, .8], [.4, .5]])
    expected = kernels.REMatchKernel(alpha=.5, threshold=1e-12).get_global_similarity(local)
    assert _rematch_similarity(local, alpha=.5) == pytest.approx(expected, abs=1e-7)


def test_rematch_nonconvergence_has_an_iteration_limit():
    with pytest.raises(RuntimeError, match="iteration limit"):
        _rematch_similarity([[1, .3], [.4, .8], [.1, .2]], alpha=.01, max_iterations=1)


def test_small_alpha_rematch_remains_finite_in_log_domain():
    assert np.isfinite(_rematch_similarity(np.eye(3), alpha=.0001))


def test_grouped_precomputed_clustering_extracts_square_submatrices():
    from pyar.selection.clusterers import _run_grouped_algorithm
    distances = pairwise_distances([[0], [1], [10], [11]])
    with mock.patch("pyar.selection.clusterers._run_algorithm", return_value=np.array([0, 0])) as run:
        labels = _run_grouped_algorithm(distances, "agglomerative", 1, {"metric": "precomputed"}, (0, 0, 1, 1))
    assert labels.tolist() == [0, 0, 1, 1]
    assert run.call_args_list[0].args[0].shape == (2, 2)
    np.testing.assert_array_equal(run.call_args_list[0].args[0], distances[:2, :2])
