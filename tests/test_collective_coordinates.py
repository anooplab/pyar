"""Tests for the interfragment collective-coordinate API."""

import unittest

import numpy as np

from pyar.biases.collective_coordinates import evaluate_contact_coordinate
from pyar.data.units import angstrom2bohr


class CollectiveCoordinateTests(unittest.TestCase):
    def setUp(self):
        self.symbols = ["H", "C", "H"]
        self.coordinates = angstrom2bohr(np.asarray([
            [0.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
            [3.0, 1.0, 0.0],
        ]))

    def test_softmin_gradient_matches_finite_difference_and_preserves_order(self):
        # Deliberately non-contiguous fragments exercise global gradient placement.
        fragments = [[2], [0, 1]]
        q, gradient, diagnostics = evaluate_contact_coordinate(
            fragments, self.symbols, self.coordinates, kind="softmin", beta=1.7
        )
        displacement = 1.0e-5
        forward = self.coordinates.copy()
        backward = self.coordinates.copy()
        forward[2, 0] += displacement
        backward[2, 0] -= displacement
        forward_q, _, _ = evaluate_contact_coordinate(
            fragments, self.symbols, forward, kind="softmin", beta=1.7
        )
        backward_q, _, _ = evaluate_contact_coordinate(
            fragments, self.symbols, backward, kind="softmin", beta=1.7
        )

        self.assertAlmostEqual(gradient[2, 0], (forward_q - backward_q) / (2.0 * displacement), places=7)
        np.testing.assert_allclose(gradient.sum(axis=0), np.zeros(3), atol=1.0e-12)
        self.assertEqual(diagnostics.closest_pair, (2, 1))
        self.assertEqual(diagnostics.pair_count, 2)
        self.assertIsNotNone(diagnostics.effective_pair_count)

    def test_afir_gradient_and_contact_diagnostics(self):
        q, gradient, diagnostics = evaluate_contact_coordinate(
            [[0], [1, 2]], self.symbols, self.coordinates, kind="afir", contact_factor=1.2
        )
        self.assertGreater(q, 0.0)
        self.assertEqual(gradient.shape, (3, 3))
        np.testing.assert_allclose(gradient.sum(axis=0), np.zeros(3), atol=1.0e-12)
        self.assertEqual(diagnostics.contact_count, sum(pair.is_contact for pair in diagnostics.contacts))
        self.assertIn("contacts", diagnostics.as_dict())

    def test_rejects_overlapping_fragments(self):
        with self.assertRaisesRegex(ValueError, "must not overlap"):
            evaluate_contact_coordinate([[0, 1], [1, 2]], self.symbols, self.coordinates)


if __name__ == "__main__":
    unittest.main()
