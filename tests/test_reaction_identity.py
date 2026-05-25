import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np

from pyar.Molecule import Molecule
from pyar import reaction_identity


class ReactionIdentityTests(unittest.TestCase):
    def test_write_disconnected_reference_separates_fragments(self):
        molecule = Molecule(
            ["C", "H"],
            np.array([[0.0, 0.0, 0.0], [1.1, 0.0, 0.0]]),
            fragments=[[0], [1]],
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            output = Path(tmpdir) / "reference.xyz"
            reaction_identity.write_disconnected_reference(molecule, output)
            reference = Molecule.from_xyz(output)

        self.assertAlmostEqual(reference.coordinates[1, 0] - reference.coordinates[0, 0], 101.1)
        self.assertAlmostEqual(molecule.coordinates[1, 0] - molecule.coordinates[0, 0], 1.1)

    def test_molecule_identity_uses_openbabel_helpers(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            reference = Path(tmpdir) / "reference.xyz"
            reference.write_text("1\nref\nH 0 0 0\n")
            with mock.patch.object(reaction_identity.babel, "make_inchi_string_from_xyz", return_value="ref-inchi"), \
                mock.patch.object(reaction_identity.babel, "make_smile_string_from_xyz", return_value="ref-smile"):
                identity = reaction_identity.molecule_identity_from_xyz(reference)

        self.assertEqual(identity, {"inchi": "ref-inchi", "smiles": "ref-smile"})

    def test_molecule_identity_rejects_incomplete_identifiers(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            xyz_path = Path(tmpdir) / "candidate.xyz"
            xyz_path.write_text("1\ncandidate\nH 0 0 0\n")
            with mock.patch.object(reaction_identity.babel, "make_inchi_string_from_xyz", return_value="inchi"), \
                mock.patch.object(reaction_identity.babel, "make_smile_string_from_xyz", return_value=""):
                with self.assertRaisesRegex(ValueError, "complete product identity"):
                    reaction_identity.molecule_identity_from_xyz(xyz_path)

    def test_same_molecular_identity_uses_inchi_as_canonical_key(self):
        first = {"inchi": "same-inchi", "smiles": "serialization-a"}
        same_product = {"inchi": "same-inchi", "smiles": "serialization-b"}
        different_product = {"inchi": "other-inchi", "smiles": "serialization-a"}

        self.assertTrue(reaction_identity.same_molecular_identity(first, same_product))
        self.assertFalse(reaction_identity.same_molecular_identity(first, different_product))


if __name__ == "__main__":
    unittest.main()
