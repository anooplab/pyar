"""Analytical internal-coordinate contracts, independent of external backends."""

from types import SimpleNamespace

import numpy as np
import pytest

from pyar.representations import make_internal_coordinates


def molecule(coordinates):
    return SimpleNamespace(
        coordinates=np.asarray(coordinates, dtype=float),
        covalent_radius=np.full(len(coordinates), 0.76),
    )


@pytest.mark.scientific
@pytest.mark.tier1
@pytest.mark.parametrize("handedness", [-1, 1])
@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("transform", [False, True])
def test_internal_torsion_uses_coordinates_and_preserves_signed_degrees(
    handedness, reverse, transform
):
    # Three orthogonal 1.4-Angstrom bonds form a chain: the nonbonded
    # diagonals exceed the existing carbon-radius connectivity cutoff.
    xyz = np.array([[0, 1.4, 0], [0, 0, 0], [1.4, 0, 0], [1.4, 0, 1.4 * handedness]])
    if reverse:
        xyz = xyz[::-1].copy()
    if transform:
        xyz = xyz @ np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]]) + [3, -2, 5]
    mol = molecule(xyz)
    original = mol.coordinates.copy()

    bonds, angles, torsions = make_internal_coordinates(mol)

    assert [entry[0] for entry in bonds] == [[0, 1], [1, 2], [2, 3]]
    np.testing.assert_allclose([entry[1] for entry in bonds], 1.4, atol=1e-14)
    assert [entry[0] for entry in angles] == [[0, 1, 2], [1, 2, 3]]
    np.testing.assert_allclose([entry[1] for entry in angles], 90.0, atol=1e-12)
    assert len(torsions) == 1
    assert torsions[0][0] == [0, 1, 2, 3]
    # Exact analytic quarter turn; tolerance only accommodates floating-point rotation.
    assert torsions[0][1] == pytest.approx(-90 * handedness, abs=1e-12)
    np.testing.assert_array_equal(mol.coordinates, original)


def test_ring_torsions_keep_atom_order_and_reverse_deduplication():
    mol = molecule([[0, 0, 0], [1.4, 0, 0], [1.4, 1.4, 0], [0, 1.4, 0]])

    bonds, angles, torsions = make_internal_coordinates(mol)

    assert len(bonds) == 4 and len(angles) == 4 and len(torsions) == 4
    np.testing.assert_allclose([entry[1] for entry in angles], 90.0, atol=1e-12)
    indices = [tuple(entry[0]) for entry in torsions]
    assert len(set(indices)) == 4
    assert not any(tuple(reversed(entry)) in indices for entry in indices)
    np.testing.assert_allclose([entry[1] for entry in torsions], 0.0, atol=1e-12)


@pytest.mark.parametrize(
    "coordinates, expected_counts",
    [
        ([[0, 0, 0]], (0, 0, 0)),
        ([[0, 0, 0], [1.4, 0, 0]], (1, 0, 0)),
        ([[0, 0, 0], [1.4, 0, 0], [0, 1.4, 0], [0, 0, 1.4]], (3, 3, 0)),
    ],
)
def test_small_or_branched_inputs_keep_the_three_list_schema(
    coordinates, expected_counts
):
    result = make_internal_coordinates(molecule(coordinates))

    assert isinstance(result, list) and len(result) == 3
    assert tuple(map(len, result)) == expected_counts
