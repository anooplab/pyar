"""Behavioral checks for graph-first, conservative iRMSD fallback."""

from types import SimpleNamespace
from unittest import mock

import numpy as np

from pyar.structure_comparison import ComparisonResult, GraphFirstDeduplicationComparator


def _molecule(name, symbols, coordinates):
    return SimpleNamespace(name=name, atoms_list=symbols,
                           coordinates=np.asarray(coordinates, dtype=float))


def _incomplete_graph_result():
    return ComparisonResult(
        True, None, None, "element-labeled-graph-rmsd", 0.1,
        {"connectivity_match": True, "comparison_complete": False,
         "isomorphisms_evaluated": 10000},
    )


def test_irmsd_fallback_only_runs_after_graph_match_with_incomplete_search():
    first = _molecule("a", ["C", "C"], [[0, 0, 0], [1.5, 0, 0]])
    second = _molecule("b", ["C", "C"], [[0, 0, 0], [1.5, 0, 0]])
    graph = mock.Mock(compare=mock.Mock(return_value=_incomplete_graph_result()))

    with mock.patch("pyar.structure_comparison.deduplication_policy.GraphRMSDComparator",
                    return_value=graph):
        with mock.patch("pyar.structure_comparison.deduplication_policy._run_bidirectional_irmsd",
                        return_value=((0.02, 0.04), "ok")) as run_irmsd:
            result = GraphFirstDeduplicationComparator(threshold=0.1).compare(first, second)

    run_irmsd.assert_called_once()
    assert result.equivalent is True
    assert result.distance == 0.04
    assert result.metadata["connectivity_match"] is True
    assert result.metadata["irmsd_distances_angstrom"] == [0.02, 0.04]


def test_irmsd_fallback_uses_conservative_maximum_of_both_directions():
    first = _molecule("a", ["C", "C"], [[0, 0, 0], [1.5, 0, 0]])
    second = _molecule("b", ["C", "C"], [[0, 0, 0], [1.5, 0, 0]])
    graph = mock.Mock(compare=mock.Mock(return_value=_incomplete_graph_result()))

    with mock.patch("pyar.structure_comparison.deduplication_policy.GraphRMSDComparator",
                    return_value=graph):
        with mock.patch("pyar.structure_comparison.deduplication_policy._run_bidirectional_irmsd",
                        return_value=((0.02, 0.11), "ok")):
            result = GraphFirstDeduplicationComparator(threshold=0.1).compare(first, second)

    assert result.distance == 0.11
    assert result.equivalent is False


def test_irmsd_diagnostics_or_failure_leave_the_comparison_uncertain():
    first = _molecule("a", ["C", "C"], [[0, 0, 0], [1.5, 0, 0]])
    second = _molecule("b", ["C", "C"], [[0, 0, 0], [1.5, 0, 0]])
    graph = mock.Mock(compare=mock.Mock(return_value=_incomplete_graph_result()))

    with mock.patch("pyar.structure_comparison.deduplication_policy.GraphRMSDComparator",
                    return_value=graph):
        with mock.patch("pyar.structure_comparison.deduplication_policy._run_bidirectional_irmsd",
                        return_value=(None, "backend_diagnostic_or_incomplete_output")):
            result = GraphFirstDeduplicationComparator(threshold=0.1).compare(first, second)

    assert result.compatible
    assert result.distance is None
    assert result.equivalent is None
    assert result.metadata["fallback_status"] == "backend_diagnostic_or_incomplete_output"


def test_irmsd_does_not_run_when_graphs_are_not_isomorphic():
    first = _molecule("a", ["C", "C"], [[0, 0, 0], [1.5, 0, 0]])
    second = _molecule("b", ["C", "C"], [[0, 0, 0], [2.0, 0, 0]])
    graph_result = ComparisonResult(
        False, None, None, "element-labeled-graph-rmsd", 0.1,
        {"connectivity_match": False, "comparison_complete": True},
    )
    graph = mock.Mock(compare=mock.Mock(return_value=graph_result))

    with mock.patch("pyar.structure_comparison.deduplication_policy.GraphRMSDComparator",
                    return_value=graph):
        with mock.patch("pyar.structure_comparison.deduplication_policy._run_bidirectional_irmsd") as run_irmsd:
            result = GraphFirstDeduplicationComparator(threshold=0.1).compare(first, second)

    run_irmsd.assert_not_called()
    assert result is graph_result


def test_real_irmsd_resolves_graph_mapping_cap_for_au13_duplicate():
    import pytest

    pytest.importorskip("irmsd")
    from test_structure_comparison_dataset import _load_xyz_frames

    frames = _load_xyz_frames()
    result = GraphFirstDeduplicationComparator(threshold=0.1).compare(
        frames["au13_ico"], frames["au13_ico_rotperm"],
    )

    assert result.equivalent is True
    assert result.metadata["comparison_complete"] is False
    assert result.metadata["fallback_method"] == "irmsd"
    assert result.metadata["fallback_status"] == "ok"
    assert result.distance < 1e-12
