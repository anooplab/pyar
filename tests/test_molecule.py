import tempfile
import unittest
from pathlib import Path

import numpy as np

from pyar.Molecule import Molecule, XYZParseError, parse_xyz, read_xyz


class MoleculeTests(unittest.TestCase):
    def test_parse_xyz_reads_energy_when_present(self):
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "a.xyz"
            path.write_text(
                "2\n"
                "job: -1.23\n"
                "H 0 0 0\n"
                "H 0 0 0.74\n"
            )
            atoms, coords, name, title, energy = parse_xyz(str(path))
            self.assertEqual(atoms, ["H", "H"])
            self.assertEqual(coords.shape, (2, 3))
            self.assertEqual(title, "job: -1.23")
            self.assertAlmostEqual(float(energy), -1.23)

    def test_parse_xyz_raises_on_bad_header(self):
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "bad.xyz"
            path.write_text("notanint\nx\n")
            with self.assertRaises(XYZParseError):
                parse_xyz(str(path))

    def test_read_xyz_exits_on_parse_failure(self):
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "bad.xyz"
            path.write_text("notanint\nx\n")
            with self.assertRaises(SystemExit):
                read_xyz(str(path))

    def test_coordinates_normalization_and_centroid_update(self):
        mol = Molecule(["H"], [[0, 0, 0]])
        self.assertEqual(mol.coordinates.shape, (1, 3))
        self.assertTrue(np.allclose(mol.centroid, [0.0, 0.0, 0.0]))
        mol.translate(np.array([1.0, 0.0, 0.0]))
        self.assertTrue(np.allclose(mol.coordinates[0], [-1.0, 0.0, 0.0]))

    def test_non_mutating_transforms_leave_original_unchanged(self):
        mol = Molecule(["H"], np.array([[0.0, 0.0, 0.0]]))
        moved = mol.translated(np.array([1.0, 2.0, 3.0]))
        self.assertTrue(np.allclose(mol.coordinates, [[0.0, 0.0, 0.0]]))
        self.assertTrue(np.allclose(moved.coordinates, [[-1.0, -2.0, -3.0]]))

    def test_to_ase_atoms_roundtrip_when_available(self):
        try:
            import ase  # noqa: F401
        except Exception:
            self.skipTest("ASE is not installed")

        mol = Molecule(["H", "H"], np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.74]]))
        atoms = mol.to_ase_atoms()
        mol2 = Molecule.from_ase_atoms(atoms)
        self.assertEqual(mol2.atoms_list, ["H", "H"])
        self.assertTrue(np.allclose(mol2.coordinates, mol.coordinates))


if __name__ == "__main__":
    unittest.main()

