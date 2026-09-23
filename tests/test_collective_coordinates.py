"""Tests for the interfragment collective-coordinate API."""

import unittest

from autograd import grad
import autograd.numpy as anp
import numpy as np

from pyar.biases.afir import get_covalent_radius
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

    def _autograd_gradient(self, fragments, symbols, coordinates, kind, *, beta=1.0, distance_power=6.0):
        """Independent derivative reference retained outside the production path."""
        left, right = (tuple(fragment) for fragment in fragments)
        radii = anp.asarray([get_covalent_radius(symbol) for symbol in symbols])

        def coordinate(all_coordinates):
            differences = all_coordinates[list(left), None, :] - all_coordinates[None, list(right), :]
            distances = anp.sqrt(anp.sum(differences ** 2, axis=2)).reshape(-1)
            radii_sums = (radii[list(left), None] + radii[None, list(right)]).reshape(-1)
            if kind == "softmin":
                scaled = -beta * (distances - radii_sums)
                maximum = anp.max(scaled)
                return -(maximum + anp.log(anp.mean(anp.exp(scaled - maximum)))) / beta
            weights = (radii_sums / distances) ** distance_power
            return anp.sum(weights * distances) / anp.sum(weights)

        return np.asarray(grad(coordinate)(coordinates), dtype=float)

    def _finite_difference_gradient(self, fragments, symbols, coordinates, kind, *, beta=1.0, distance_power=6.0):
        displacement = 1.0e-5
        reference = np.zeros_like(coordinates)
        for atom_index in range(len(coordinates)):
            for component in range(3):
                forward = coordinates.copy()
                backward = coordinates.copy()
                forward[atom_index, component] += displacement
                backward[atom_index, component] -= displacement
                forward_q, _, _ = evaluate_contact_coordinate(
                    fragments, symbols, forward, kind=kind, beta=beta, distance_power=distance_power
                )
                backward_q, _, _ = evaluate_contact_coordinate(
                    fragments, symbols, backward, kind=kind, beta=beta, distance_power=distance_power
                )
                reference[atom_index, component] = (forward_q - backward_q) / (2.0 * displacement)
        return reference

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

    def test_analytical_gradients_match_autograd_and_finite_differences(self):
        symmetric_symbols = ["C", "H", "H", "C"]
        symmetric_coordinates = angstrom2bohr(np.asarray([
            [-1.5, 0.0, 0.0],
            [-1.5, 1.0, 0.0],
            [1.5, 1.0, 0.0],
            [1.5, 0.0, 0.0],
        ]))
        cases = (
            (self.symbols, self.coordinates, [[2], [0, 1]]),
            (symmetric_symbols, symmetric_coordinates, [[3, 1], [0, 2]]),
        )
        parameters_by_kind = {
            "softmin": ({"beta": 0.4}, {"beta": 1.7}, {"beta": 8.0}),
            "afir": ({"distance_power": 2.0}, {"distance_power": 6.0}, {"distance_power": 10.0}),
        }
        for symbols, coordinates, fragments in cases:
            for kind, parameters_list in parameters_by_kind.items():
                for parameters in parameters_list:
                    with self.subTest(kind=kind, parameters=parameters, fragments=fragments):
                        _, analytical, _ = evaluate_contact_coordinate(
                            fragments, symbols, coordinates, kind=kind, **parameters
                        )
                        autograd_reference = self._autograd_gradient(
                            fragments, symbols, coordinates, kind, **parameters
                        )
                        finite_difference_reference = self._finite_difference_gradient(
                            fragments, symbols, coordinates, kind, **parameters
                        )
                        np.testing.assert_allclose(
                            analytical, autograd_reference, rtol=1.0e-10, atol=1.0e-11
                        )
                        np.testing.assert_allclose(
                            analytical, finite_difference_reference, rtol=1.0e-6, atol=1.0e-7
                        )

    def test_seeded_random_geometries_match_independent_gradients(self):
        rng = np.random.default_rng(20260923)
        symbols = ["C", "H", "N", "O", "C", "H"]
        fragments = [[0, 2, 4], [1, 3, 5]]
        for case_index in range(5):
            coordinates = rng.normal(size=(len(symbols), 3)) * 1.8
            # Keep fragments separated while retaining varied, non-planar geometries.
            coordinates[[1, 3, 5], 0] += 3.5
            for kind, parameters in (
                ("afir", {"distance_power": (2.0, 6.0, 10.0)[case_index % 3]}),
                ("softmin", {"beta": (0.4, 1.7, 8.0)[case_index % 3]}),
            ):
                with self.subTest(case=case_index, kind=kind):
                    _, analytical, _ = evaluate_contact_coordinate(
                        fragments, symbols, coordinates, kind=kind, **parameters
                    )
                    autograd_reference = self._autograd_gradient(
                        fragments, symbols, coordinates, kind, **parameters
                    )
                    finite_difference_reference = self._finite_difference_gradient(
                        fragments, symbols, coordinates, kind, **parameters
                    )
                    np.testing.assert_allclose(analytical, autograd_reference, rtol=1e-9, atol=1e-10)
                    np.testing.assert_allclose(analytical, finite_difference_reference, rtol=2e-6, atol=2e-7)

    def test_coordinates_and_gradients_are_rigid_motion_covariant(self):
        angle = 0.73
        rotation = np.asarray([
            [np.cos(angle), -np.sin(angle), 0.0],
            [np.sin(angle), np.cos(angle), 0.0],
            [0.0, 0.0, 1.0],
        ])
        translation = np.asarray([1.2, -0.7, 2.3])
        for kind, parameters in (("afir", {"distance_power": 4.0}),
                                 ("softmin", {"beta": 2.1})):
            with self.subTest(kind=kind):
                q, gradient, _ = evaluate_contact_coordinate(
                    [[2], [0, 1]], self.symbols, self.coordinates, kind=kind, **parameters
                )
                translated_q, translated_gradient, _ = evaluate_contact_coordinate(
                    [[2], [0, 1]], self.symbols, self.coordinates + translation,
                    kind=kind, **parameters
                )
                rotated_coordinates = self.coordinates @ rotation.T
                rotated_q, rotated_gradient, _ = evaluate_contact_coordinate(
                    [[2], [0, 1]], self.symbols, rotated_coordinates, kind=kind, **parameters
                )
                self.assertAlmostEqual(q, translated_q, places=12)
                self.assertAlmostEqual(q, rotated_q, places=12)
                np.testing.assert_allclose(gradient, translated_gradient, atol=1e-12)
                np.testing.assert_allclose(rotated_gradient, gradient @ rotation.T, atol=1e-11)

    def test_rejects_overlapping_fragments(self):
        with self.assertRaisesRegex(ValueError, "must not overlap"):
            evaluate_contact_coordinate([[0, 1], [1, 2]], self.symbols, self.coordinates)


if __name__ == "__main__":
    unittest.main()
