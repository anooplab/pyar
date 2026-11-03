"""Analytical A-B-C angles and their existing hydrogen-bond consumer."""

import numpy as np
import pytest

from pyar.property import calculate_angle, hydrogen_bond_analysis

pytestmark = [pytest.mark.scientific, pytest.mark.tier1]


@pytest.mark.parametrize("expected", [0.0, 30.0, 60.0, 90.0, 120.0, 180.0])
@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("transform", [False, True])
def test_bond_angle_is_at_the_middle_atom_in_degrees(expected, reverse, transform):
    theta = np.deg2rad(expected)
    # Unequal ray lengths keep even the zero-angle endpoints distinct.
    points = np.array(
        [
            [1.4, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [2.1 * np.cos(theta), 2.1 * np.sin(theta), 0.0],
        ]
    )
    if reverse:
        points = points[::-1].copy()
    if transform:
        points = 3 * points @ np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]]) + [3, -2, 5]
    original = points.copy()

    actual = calculate_angle(*points)

    # Analytical degrees; tolerance covers floating-point coordinate transformations.
    assert actual == pytest.approx(expected, abs=1e-12)
    np.testing.assert_array_equal(points, original)


def test_parallel_rays_do_not_produce_nan_from_roundoff():
    with np.errstate(invalid="raise"):
        angle = calculate_angle(np.ones(3), np.zeros(3), 2 * np.ones(3))
    assert angle == pytest.approx(0.0, abs=1e-12)


@pytest.mark.parametrize("expected_angle", [90.0, 159.0, 161.0, 180.0])
@pytest.mark.parametrize("separation", [2.0, 2.25])
def test_hydrogen_bond_analysis_uses_donor_hydrogen_acceptor_angle(
    capsys, expected_angle, separation
):
    direction = np.deg2rad(180.0 - expected_angle)
    coordinates = np.array(
        [
            [-1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [separation * np.cos(direction), separation * np.sin(direction), 0.0],
        ]
    )
    qualifies = expected_angle > 160.0
    original = coordinates.copy()

    result = hydrogen_bond_analysis(
        coordinates, np.array([0.66, 0.31, 0.66]), [8, 1, 8], ["O", "H", "O"]
    )

    assert result.shape == (3, 3)
    assert result.dtype.kind == "f"
    assert result[1, 2] == pytest.approx(separation if qualifies else 0.0, abs=1e-14)
    output = capsys.readouterr().out
    assert ("H(2) - O(3) =" in output) is qualifies
    np.testing.assert_array_equal(coordinates, original)
