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

    def test_parse_xyz_reads_energy_from_whitespace_title(self):
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "a.xyz"
            path.write_text(
                "1\n"
                "job -1.25\n"
                "H 0 0 0\n"
            )
            self.assertAlmostEqual(parse_xyz(str(path))[-1], -1.25)

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

    def test_constructor_rejects_mismatched_atoms_and_coordinates(self):
        with self.assertRaisesRegex(ValueError, "same number of atoms"):
            Molecule(["H", "H"], [[0, 0, 0]])

    def test_constructor_rejects_invalid_fragment_indices(self):
        with self.assertRaisesRegex(ValueError, "fragment atom index"):
            Molecule(["H"], [[0, 0, 0]], fragments=[[1]])

    def test_coordinates_can_be_cleared_after_failed_job(self):
        mol = Molecule(["H"], [[0, 0, 0]])
        mol.coordinates = None
        self.assertIsNone(mol.coordinates)
        self.assertIsNone(mol.centroid)

    def test_copy_preserves_failed_job_without_coordinates(self):
        mol = Molecule(["H"], [[0, 0, 0]], fragments=[[0]])
        mol.coordinates = None

        copied = mol.copy()

        self.assertIsNone(copied.coordinates)
        self.assertEqual(copied.number_of_atoms, 1)
        self.assertEqual(copied.fragments, [[0]])

    def test_centroid_reflects_in_place_coordinate_changes(self):
        mol = Molecule(["H", "H"], [[0, 0, 0], [2, 0, 0]])
        mol.coordinates[:] += np.array([4.0, 0.0, 0.0])
        self.assertTrue(np.allclose(mol.centroid, [5.0, 0.0, 0.0]))

    def test_moments_of_inertia_tensor_does_not_translate_molecule(self):
        mol = Molecule(["H", "H"], [[0, 0, 0], [2, 0, 0]])
        before = np.array(mol.coordinates, copy=True)
        _ = mol.moments_of_inertia_tensor
        self.assertTrue(np.allclose(mol.coordinates, before))

    def test_repr_is_human_readable(self):
        mol = Molecule(["H"], [[0, 0, 0]], name="example", energy=-1.23)
        self.assertEqual(
            repr(mol),
            "Molecule(name='example', atoms=1, energy=-1.23)",
        )

    def test_non_mutating_transforms_leave_original_unchanged(self):
        mol = Molecule(["H"], np.array([[0.0, 0.0, 0.0]]))
        moved = mol.translated(np.array([1.0, 2.0, 3.0]))
        self.assertTrue(np.allclose(mol.coordinates, [[0.0, 0.0, 0.0]]))
        self.assertTrue(np.allclose(moved.coordinates, [[-1.0, -2.0, -3.0]]))

    def test_rotate_3d_preserves_centroid(self):
        mol = Molecule(["H", "H"], [[2.0, 1.0, 0.0], [4.0, 1.0, 0.0]])
        centroid = mol.centroid.copy()
        mol.rotate_3d((0.0, 0.0, np.pi / 2))
        self.assertTrue(np.allclose(mol.centroid, centroid))

    def test_align_centers_molecule_at_centre_of_mass(self):
        mol = Molecule(["C", "H"], [[3.0, 0.0, 0.0], [4.0, 1.0, 0.0]])
        mol.align()
        centre_of_mass = np.average(mol.coordinates, axis=0, weights=mol.atomic_mass)
        self.assertTrue(np.allclose(centre_of_mass, [0.0, 0.0, 0.0]))

    def test_merged_with_is_explicit_addition_api(self):
        first = Molecule(["C"], [[0.0, 0.0, 0.0]])
        second = Molecule(["H"], [[1.0, 0.0, 0.0]])
        merged = first.merged_with(second)
        self.assertEqual(merged.fragments, [[0], [1]])
        self.assertEqual((first + second).atoms_list, merged.atoms_list)

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
