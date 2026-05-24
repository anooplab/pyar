import tempfile
import unittest
from pathlib import Path

import numpy as np

from pyar.Molecule import Molecule
from pyar.molecule_geometry import moments_of_inertia_tensor, rotate_3d
from pyar.molecule_io import (
    as_coordinates_array,
    parse_xyz,
    validate_fragments,
    write_turbomole_coord,
    write_xyz,
)
from pyar.molecule_merge import combine_multiplicity, merged_with


class MoleculeHelperTests(unittest.TestCase):
    def test_write_xyz_roundtrip(self):
        mol = Molecule(["H", "H"], [[0.0, 0.0, 0.0], [0.0, 0.0, 0.74]], name="h2", title="h2")
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "a.xyz"
            write_xyz(mol, path)
            atoms, coords, name, title, energy = parse_xyz(str(path))
            self.assertEqual(atoms, ["H", "H"])
            self.assertEqual(name, str(path)[:-4])
            self.assertEqual(title, "h2: None")
            self.assertIsNone(energy)
            self.assertTrue(np.allclose(coords, mol.coordinates))

    def test_as_coordinates_array_normalizes_shape(self):
        coords = as_coordinates_array([[0.0, 0.0, 0.0]])
        self.assertEqual(coords.shape, (1, 3))
        with self.assertRaisesRegex(ValueError, "coordinates must have shape"):
            as_coordinates_array([0.0, 0.0, 0.0])

    def test_validate_fragments_rejects_out_of_range_indices(self):
        self.assertEqual(validate_fragments([[0], [1, 2]], 3), [[0], [1, 2]])
        with self.assertRaisesRegex(ValueError, "fragment atom index"):
            validate_fragments([[3]], 3)

    def test_write_turbomole_coord_creates_expected_file(self):
        mol = Molecule(["H"], [[0.0, 0.0, 0.0]])
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "coord"
            write_turbomole_coord(mol, path)
            self.assertIn("$coord", path.read_text())
            self.assertIn("$end", path.read_text())

    def test_moments_of_inertia_tensor_does_not_translate(self):
        mol = Molecule(["H", "H"], [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
        before = np.array(mol.coordinates, copy=True)
        _ = moments_of_inertia_tensor(mol)
        self.assertTrue(np.allclose(mol.coordinates, before))

    def test_rotate_3d_preserves_centroid(self):
        mol = Molecule(["H", "H"], [[2.0, 1.0, 0.0], [4.0, 1.0, 0.0]])
        before = np.array(mol.centroid, copy=True)
        rotate_3d(mol, (0.0, 0.0, np.pi / 2))
        self.assertTrue(np.allclose(mol.centroid, before))

    def test_combine_multiplicity_uses_existing_heuristic(self):
        self.assertEqual(combine_multiplicity(1, 2), 2)
        self.assertEqual(combine_multiplicity(2, 2), 1)

    def test_merged_with_preserves_fragment_structure(self):
        first = Molecule(["C"], [[0.0, 0.0, 0.0]], name="a")
        second = Molecule(["H"], [[1.0, 0.0, 0.0]], name="b")
        merged = merged_with(first, second)
        self.assertEqual(merged.fragments, [[0], [1]])
        self.assertEqual(merged.name, "a + b")
        self.assertEqual(merged.scftype, "rhf")


if __name__ == "__main__":
    unittest.main()
