"""Tests for the soft-minimum intermolecular bias potential."""

import unittest

import numpy as np

from pyar.biases.softmin import softmin
from pyar.data.units import angstrom2bohr


class SoftminBiasTests(unittest.TestCase):
    def setUp(self):
        self.fragment_indices = [[0], [1, 2]]
        self.atoms = ["C", "H", "H"]
        self.coordinates = angstrom2bohr(
            np.asarray(
                [
                    [0.0, 0.0, 0.0],
                    [3.0, 0.0, 0.0],
                    [3.0, 1.0, 0.0],
                ]
            )
        )

    def test_zero_force_has_zero_energy_and_forces(self):
        energy, forces = softmin(
            self.fragment_indices, self.atoms, self.coordinates, gamma=0.0
        )

        self.assertEqual(energy, 0.0)
        np.testing.assert_allclose(forces, np.zeros((3, 3)))

    def test_bias_attracts_fragments_and_conserves_net_force(self):
        _, forces = softmin(
            self.fragment_indices, self.atoms, self.coordinates, gamma=100.0
        )

        self.assertGreater(forces[0, 0], 0.0)
        self.assertLess(forces[1:, 0].sum(), 0.0)
        np.testing.assert_allclose(forces.sum(axis=0), np.zeros(3), atol=1.0e-12)

    def test_reported_forces_match_energy_finite_difference(self):
        displacement = 1.0e-5
        energy, forces = softmin(
            self.fragment_indices, self.atoms, self.coordinates, gamma=100.0, beta=1.7
        )
        shifted = self.coordinates.copy()
        shifted[0, 0] += displacement
        shifted_energy, _ = softmin(
            self.fragment_indices, self.atoms, shifted, gamma=100.0, beta=1.7
        )

        numerical_force = -(shifted_energy - energy) / displacement
        self.assertAlmostEqual(forces[0, 0], numerical_force, places=7)

    def test_invalid_beta_is_rejected(self):
        for beta in (0.0, -1.0, "nan", "not-a-number"):
            with self.subTest(beta=beta):
                with self.assertRaisesRegex(ValueError, "(?i)soft-min beta"):
                    softmin(
                        self.fragment_indices, self.atoms, self.coordinates, gamma=100.0, beta=beta
                    )


if __name__ == "__main__":
    unittest.main()
