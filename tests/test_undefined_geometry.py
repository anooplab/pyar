"""Owner-approved exact degeneracy policy; no near-degeneracy tolerance."""

from types import SimpleNamespace

import numpy as np
import pytest

from pyar.property import calculate_angle, calculate_dihedral, hydrogen_bond_analysis
from pyar.representations import make_internal_coordinates

pytestmark = [pytest.mark.scientific, pytest.mark.tier1]


@pytest.mark.parametrize(
    "coordinates",
    [
        [[0, 0, 0], [0, 0, 0], [1, 0, 0]],
        [[1, 0, 0], [0, 0, 0], [0, 0, 0]],
        [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
    ],
)
def test_zero_angle_ray_raises_value_error_without_mutating_input(coordinates):
    points = np.array(coordinates, dtype=float)
    before = points.copy()
    with pytest.raises(ValueError, match="zero-length"):
        calculate_angle(*points)
    np.testing.assert_array_equal(points, before)


@pytest.mark.parametrize(
    "coordinates",
    [
        [[1, 0, 0], [0, 0, 0], [0, 0, 0], [0, 1, 0]],
        [[0, 0, -1], [0, 0, 0], [0, 0, 1], [1, 0, 1]],
        [[1, 0, 0], [0, 0, 0], [0, 0, 1], [0, 0, 2]],
        [[0, 0, -1], [0, 0, 0], [0, 0, 1], [0, 0, 2]],
        [[0, 0, 0], [0, 0, 0], [0, 0, 1], [1, 0, 1]],
        [[1, 0, 0], [0, 0, 0], [0, 0, 1], [0, 0, 1]],
        # Raw cross products are exactly zero; projection roundoff must not define a plane.
        [[-1, -1, -1], [0, 0, 0], [1, 1, 1], [1, 2, 1]],
    ],
)
def test_undefined_torsion_raises_value_error_without_mutating_input(coordinates):
    points = np.array(coordinates, dtype=float)
    before = points.copy()
    with pytest.raises(ValueError, match="torsion"):
        calculate_dihedral(*points)
    np.testing.assert_array_equal(points, before)


@pytest.mark.parametrize("expected", [-150.0, -90.0, -30.0, 0.0, 30.0, 90.0, 150.0])
@pytest.mark.parametrize("transform", [False, True])
def test_defined_signed_torsions_preserve_existing_convention(expected, transform):
    theta = np.deg2rad(expected)
    points = np.array(
        [
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
            [np.cos(theta), -np.sin(theta), 1.0],
        ]
    )
    if transform:
        points = 3 * points @ np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]]) + [3, -2, 5]
    before = points.copy()
    assert calculate_dihedral(*points) == pytest.approx(expected, abs=1e-12)
    np.testing.assert_array_equal(points, before)


def test_defined_planar_half_turn_remains_180_degrees():
    points = np.array(
        [[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [-1.0, 0.0, 1.0]]
    )
    assert abs(calculate_dihedral(*points)) == 180.0


@pytest.mark.parametrize("nonzero_component", [1.0e-12, 1.0e-200])
def test_nearly_collinear_torsion_is_not_rejected_by_an_added_tolerance(
    nonzero_component,
):
    points = np.array(
        [
            [nonzero_component, 0.0, -1.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 1.0, 1.0],
        ]
    )
    assert calculate_dihedral(*points) == pytest.approx(-90.0, abs=1e-12)


@pytest.mark.parametrize(
    "expected,endpoint", [(0.0, [1.0, 0.0, 0.0]), (180.0, [-1.0, 0.0, 0.0])]
)
def test_short_nonzero_angle_ray_is_defined(expected, endpoint):
    assert (
        calculate_angle(np.array([1.0e-100, 0.0, 0.0]), np.zeros(3), np.array(endpoint))
        == expected
    )


def test_internal_coordinate_consumer_propagates_undefined_torsion():
    coordinates = np.array(
        [[0.0, 0.0, 0.0], [1.4, 0.0, 0.0], [2.8, 0.0, 0.0], [4.2, 0.0, 0.0]]
    )
    molecule = SimpleNamespace(
        coordinates=coordinates, covalent_radius=np.full(4, 0.76)
    )
    before = coordinates.copy()
    with pytest.raises(ValueError, match="torsion"):
        make_internal_coordinates(molecule)
    np.testing.assert_array_equal(coordinates, before)


def test_hydrogen_bond_consumer_does_not_publish_an_undefined_angle(capsys):
    coordinates = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
    with pytest.raises(ValueError, match="zero-length"):
        hydrogen_bond_analysis(
            coordinates, np.array([0.66, 0.31, 0.66]), [8, 1, 8], ["O", "H", "O"]
        )
    assert capsys.readouterr().out == ""
