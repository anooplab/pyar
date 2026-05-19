import importlib
import sys
import unittest
from pathlib import Path

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


class FormulaGeneratorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.aggregator = importlib.import_module("pyar.aggregator")

    def test_generate_molecule_from_formula_builds_valid_geometry(self):
        rng = np.random.default_rng(0)
        molecule = self.aggregator.generate_molecule_from_formula("H4", rng=rng)

        self.assertEqual(molecule.name, "H4")
        self.assertEqual(molecule.title, "H4")
        self.assertEqual(molecule.atoms_list, ["H", "H", "H", "H"])
        self.assertEqual(molecule.coordinates.shape, (4, 3))
        self.assertTrue(np.isfinite(molecule.coordinates).all())

        min_distance = min(
            np.linalg.norm(molecule.coordinates[i] - molecule.coordinates[j])
            for i in range(molecule.number_of_atoms)
            for j in range(i + 1, molecule.number_of_atoms)
        )
        self.assertGreater(min_distance, 0.5)

    def test_generate_molecule_from_formula_accepts_lowercase_input(self):
        molecule = self.aggregator.generate_molecule_from_formula("h4")

        self.assertEqual(molecule.atoms_list, ["H", "H", "H", "H"])

    def test_generate_molecule_from_formula_accepts_mixed_formula(self):
        molecule = self.aggregator.generate_molecule_from_formula("c5h4")

        self.assertEqual(molecule.atoms_list, ["C"] * 5 + ["H"] * 4)

    def test_generate_molecule_from_formula_rejects_invalid_formula(self):
        with self.assertRaises(ValueError):
            self.aggregator.generate_molecule_from_formula("H4-")

    def test_generate_molecule_from_formula_rejects_unknown_symbol(self):
        with self.assertRaisesRegex(ValueError, "Unknown element symbol: Xx"):
            self.aggregator.generate_molecule_from_formula("Xx2")


if __name__ == "__main__":
    unittest.main()
