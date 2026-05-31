import unittest

import numpy as np

from pyar.core.molecule import Molecule
from pyar.sampling import trial_generator as trial_generation


class TrialGenerationTests(unittest.TestCase):
    def test_broken_remains_available_from_canonical_module(self):
        separated = Molecule(
            ["H", "H"],
            np.array([[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]]),
        )

        self.assertTrue(trial_generation.broken(separated))

    def test_placement_finishes_when_radial_approach_misses_seed_atoms(self):
        seed = Molecule(
            ["C", "C"],
            np.array([[0.0, 3.0, 0.0], [0.0, -3.0, 0.0]]),
        )
        monomer = Molecule(["C"], np.array([[0.0, 0.0, 0.0]]))
        vector = np.array([1.0, 0.0, 0.0, 0.0, 0.0, 0.0])

        merged = trial_generation.merge_two_molecules(
            vector,
            seed,
            monomer,
            distance_scaling=1.2,
        )

        placed = merged.coordinates[-1]
        minimum_distance = min(
            np.linalg.norm(placed - coordinate)
            for coordinate in merged.coordinates[:-1]
        )
        threshold = (seed.covalent_radius[0] + monomer.covalent_radius[0]) * 1.2
        self.assertEqual(merged.number_of_atoms, 3)
        self.assertGreaterEqual(minimum_distance, threshold)
        self.assertLess(minimum_distance, threshold + 0.06)


if __name__ == "__main__":
    unittest.main()
