import numpy as np
import pytest

from pyar.selection.distances import pairwise_distances


@pytest.mark.parametrize("metric", ["euclidean", "manhattan", "cosine"])
def test_distance_matrix_is_symmetric_finite_and_has_zero_diagonal(metric):
    result = pairwise_distances([[0, 0], [1, 0], [1, 1]], metric)
    np.testing.assert_allclose(result, result.T)
    np.testing.assert_allclose(np.diag(result), 0.0)
    assert np.isfinite(result).all()


def test_cosine_distance_defines_zero_vector_pairs_without_nan():
    result = pairwise_distances([[0, 0], [0, 0], [1, 0]], "cosine")
    assert result[0, 1] == pytest.approx(0.0)
    assert result[0, 2] == pytest.approx(1.0)


def test_unknown_distance_metric_fails_explicitly():
    with pytest.raises(ValueError, match="Unknown distance metric"):
        pairwise_distances([[0.0], [1.0]], "rmsd")
