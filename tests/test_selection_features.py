from types import SimpleNamespace
from unittest import mock

import numpy as np
import pytest

from pyar.selection.features import compute_feature_matrix, standardize_features


def _molecule(name, coordinates, atoms=("H", "H")):
    return SimpleNamespace(name=name, atoms_list=list(atoms), coordinates=np.asarray(coordinates, dtype=float))


def test_distance_histogram_is_rigid_motion_and_atom_order_invariant():
    atoms = ("H", "H", "H")
    original = _molecule("a", [[0, 0, 0], [1, 0, 0], [0, 2, 0]], atoms)
    transformed = _molecule("b", [[4, 5, 6], [4, 7, 6], [5, 5, 6]], atoms)
    reordered = _molecule("c", [[0, 2, 0], [0, 0, 0], [1, 0, 0]], atoms)

    values = compute_feature_matrix([original, transformed, reordered], "distance-histogram").values

    np.testing.assert_allclose(values[0], values[1])
    np.testing.assert_allclose(values[0], values[2])


def test_feature_fallback_records_each_failed_attempt_and_keeps_pool_rows():
    molecules = [_molecule("a", [[0, 0, 0], [1, 0, 0]]), _molecule("b", [[0, 0, 0], [2, 0, 0]])]
    with mock.patch("pyar.representations.mbtr_descriptor", side_effect=RuntimeError("MBTR unavailable")):
        with mock.patch("pyar.representations.soap_structure_descriptor", side_effect=RuntimeError("SOAP unavailable")):
            result = compute_feature_matrix(molecules, "mbtr")

    assert result.name == "distance-histogram"
    assert result.values.shape[0] == len(molecules)
    assert [item["feature"] for item in result.fallbacks] == ["mbtr", "soap"]
    assert all("unavailable" in item["reason"] for item in result.fallbacks)


def test_mixed_composition_pool_uses_one_consistent_species_vocabulary():
    molecules = [_molecule("water", [[0, 0, 0], [1, 0, 0], [0, 1, 0]], ("H", "H", "O")),
                 _molecule("hydrogen", [[0, 0, 0], [1, 0, 0]])]
    result = compute_feature_matrix(molecules, "distance-histogram")
    assert result.species == ("H", "O")
    assert result.values.shape[0] == 2
    assert np.isfinite(result.values).all()


def test_constant_features_do_not_invent_diversity():
    np.testing.assert_array_equal(standardize_features(np.ones((3, 4))), np.zeros((3, 1)))


def test_near_constant_large_features_do_not_amplify_roundoff():
    values = np.full((4, 2), 1e6) + np.array([[0], [1e-8], [-1e-8], [2e-8]])
    np.testing.assert_array_equal(standardize_features(values), np.zeros((4, 1)))


def test_small_resolved_feature_variation_is_still_scaled():
    result = standardize_features([[1e-6], [2e-6], [3e-6]])
    assert result.shape == (3, 1)
    assert result.std() == pytest.approx(1.)


def test_unknown_feature_fails_with_supported_names():
    with pytest.raises(ValueError, match="mbtr, soap, distance-histogram"):
        compute_feature_matrix([_molecule("a", [[0, 0, 0], [1, 0, 0]])], "ani")
