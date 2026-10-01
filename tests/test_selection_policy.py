from types import SimpleNamespace
from unittest import mock

import numpy as np

from pyar.selection.features import compute_feature_matrix
from pyar.selection.policy import classify_system_pool, resolve_clustering_policy


def _molecule(atoms, coordinates):
    return SimpleNamespace(atoms_list=atoms, coordinates=np.asarray(coordinates, dtype=float))


def test_homonuclear_carbon_classification_stays_unknown():
    pool = [_molecule(["C", "C"], [[0, 0, 0], [1.4, 0, 0]])]
    result = classify_system_pool(pool)
    assert result.system_type == "unknown"
    assert result.confidence == "low"
    assert "specify --system-type" in result.reason


def test_explicit_system_class_controls_feature_policy():
    assert resolve_clustering_policy("atomic-cluster", "auto", "auto")["feature"] == "soap"
    assert resolve_clustering_policy("molecular", "auto", "auto")["feature"] == "mbtr"


def test_automatic_atomic_cluster_feature_falls_back_and_records_policy():
    pool = [_molecule(["Au", "Au"], [[0, 0, 0], [2.8, 0, 0]])]
    calls = []

    def compute(molecules, feature, species):
        calls.append(feature)
        if feature == "soap":
            raise RuntimeError("SOAP unavailable")
        return np.ones((len(molecules), 3))

    with mock.patch("pyar.selection.features._compute_feature_matrix", side_effect=compute):
        result = compute_feature_matrix(pool)
    assert calls == ["soap", "mbtr"]
    assert result.name == "mbtr"
    assert result.system_type == "atomic-cluster"
    assert result.fallbacks[0]["feature"] == "soap"


def test_auto_fallback_does_not_replace_explicit_feature_order():
    pool = [_molecule(["H", "H"], [[0, 0, 0], [0.74, 0, 0]])]
    calls = []

    def compute(molecules, feature, species):
        calls.append(feature)
        if feature == "mbtr":
            raise RuntimeError("MBTR unavailable")
        return np.ones((len(molecules), 3))

    with mock.patch("pyar.selection.features._compute_feature_matrix", side_effect=compute):
        result = compute_feature_matrix(pool, "mbtr")
    assert calls == ["mbtr", "soap"]
    assert result.name == "soap"
