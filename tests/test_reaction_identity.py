import os
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np

from pyar.core.molecule import Molecule
from pyar import reaction_identity
from pyar.workflows import reaction


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

    def test_reaction_product_changed_uses_either_canonical_identifier(self):
        start = {"inchi": "same-inchi", "smiles": "same-smiles"}
        changed_smiles = {"inchi": "same-inchi", "smiles": "changed-smiles"}
        changed_inchi = {"inchi": "changed-inchi", "smiles": "same-smiles"}
        unchanged = {"inchi": "same-inchi", "smiles": "same-smiles"}

        self.assertTrue(reaction_identity.reaction_product_changed(start, changed_smiles))
        self.assertTrue(reaction_identity.reaction_product_changed(start, changed_inchi))
        self.assertFalse(reaction_identity.reaction_product_changed(start, unchanged))

    def test_separated_reactant_identity_uses_disconnected_reference(self):
        reactant_a = Molecule(["H"], np.array([[0.0, 0.0, 0.0]]), name="a")
        reactant_b = Molecule(["H"], np.array([[0.0, 0.0, 1.0]]), name="b")
        expected = {"inchi": "reactants-inchi", "smiles": "reactants-smiles"}

        with mock.patch.object(
            reaction_identity,
            "molecule_identity_from_xyz",
            return_value=expected,
        ) as identity_from_xyz:
            identity = reaction_identity.separated_reactant_identity(
                reactant_a,
                reactant_b,
            )

        self.assertEqual(identity, expected)
        identity_from_xyz.assert_called_once()

    def test_orientation_final_coordinate_path_prefers_relaxed_snapshot(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            gamma_directory = Path(tmpdir)
            orientation_directory = gamma_directory / "orientation_abcdefgh"
            optimized_directory = orientation_directory / "job_0001_abcdefgh"
            optimized_directory.mkdir(parents=True)
            (orientation_directory / "result_relax.xyz").write_text("1\nrelaxed\nH 0 0 0\n")
            (optimized_directory / "result_0001_abcdefgh.xyz").write_text(
                "1\noptimized\nH 0 0 0\n"
            )

            final_coordinate = reaction._orientation_final_coordinate_path(
                str(gamma_directory),
                "orientation_abcdefgh",
                "0001_abcdefgh",
            )

        self.assertEqual(final_coordinate, "orientation_abcdefgh/result_relax.xyz")

    def test_orientation_final_coordinate_path_uses_relax_job_when_snapshot_is_missing(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            gamma_directory = Path(tmpdir)
            orientation_directory = gamma_directory / "orientation_abcdefgh"
            relax_directory = orientation_directory / "job_relax"
            biased_directory = orientation_directory / "job_0001_abcdefgh"
            relax_directory.mkdir(parents=True)
            biased_directory.mkdir()
            (relax_directory / "result_relax.xyz").write_text("1\nrelax job\nH 0 0 0\n")
            (biased_directory / "result_0001_abcdefgh.xyz").write_text(
                "1\nbiased job\nH 0 0 0\n"
            )

            final_coordinate = reaction._orientation_final_coordinate_path(
                str(gamma_directory),
                "orientation_abcdefgh",
                "0001_abcdefgh",
            )

        self.assertEqual(final_coordinate, "orientation_abcdefgh/job_relax/result_relax.xyz")

    def test_reaction_optimize_all_logs_existing_biased_result_for_survivor(self):
        class DummyReactionMolecule:
            def __init__(self):
                self.name = "trial_abcdefgh"
                self.energy = 0.0
                self.fragments = [[0], [1]]
                self.coordinates = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])

            def mol_to_xyz(self, filename):
                Path(filename).write_text("1\norientation\nH 0 0 0\n")

            def copy(self):
                copied = DummyReactionMolecule()
                copied.name = self.name
                copied.energy = self.energy
                copied.coordinates = self.coordinates.copy()
                return copied

            def is_bonded(self):
                return False

        class DummyRunState:
            def current_survivor_molecules(self):
                return []

            def record_job(self, *args, **kwargs):
                return None

        def fake_optimise(molecule, qc_params):
            result_directory = Path(f"job_{molecule.name}")
            result_directory.mkdir()
            (result_directory / f"result_{molecule.name}.xyz").write_text(
                "1\noptimized\nH 0 0 0\n"
            )
            molecule.energy = -1.0
            return True

        with tempfile.TemporaryDirectory() as tmpdir:
            cwd = Path.cwd()
            try:
                os.chdir(tmpdir)
                product_dir = Path(tmpdir) / "products"
                product_dir.mkdir()
                with mock.patch.object(
                    reaction,
                    "molecule_identity_from_xyz",
                    return_value={"inchi": "reactants-inchi", "smiles": "reactants-smiles"},
                ), \
                    mock.patch.object(reaction, "optimise", side_effect=fake_optimise):
                    with self.assertLogs("pyar.workflows.reaction", level="INFO") as logs:
                        reaction.optimize_all(
                            gamma_id="0001",
                            orientations=[DummyReactionMolecule()],
                            run_state=DummyRunState(),
                            product_dir=str(product_dir),
                            qc_param={"gamma": 1.0},
                        )
            finally:
                os.chdir(cwd)

        logged_messages = "\n".join(logs.output)
        self.assertIn(
            "final_coordinate=orientation_abcdefgh/job_0001_abcdefgh/"
            "result_0001_abcdefgh.xyz",
            logged_messages,
        )
        self.assertNotIn("final_coordinate=orientation_abcdefgh/result_relax.xyz", logged_messages)

    def test_reaction_optimize_all_does_not_use_connectivity_filter(self):
        class DummyReactionMolecule:
            def __init__(self):
                self.name = "orientation"
                self.energy = 0.0
                self.fragments = [[0]]

            def mol_to_xyz(self, filename):
                Path(filename).write_text("1\norientation\nH 0 0 0\n")

            def copy(self):
                return DummyReactionMolecule()

            def is_bonded(self):
                return True

        class DummyRunState:
            def current_survivor_molecules(self):
                return []

            def record_job(self, *args, **kwargs):
                return None

            def record_product(self, *args, **kwargs):
                return None

        molecule = DummyReactionMolecule()
        identity = {"inchi": "same-inchi", "smiles": "same-smiles"}

        with tempfile.TemporaryDirectory() as tmpdir:
            cwd = Path.cwd()
            try:
                with mock.patch.object(reaction, "optimise", return_value=True):
                    with mock.patch.object(reaction, "relax_without_afir_bias", return_value=True):
                        with mock.patch.object(reaction, "molecule_identity_from_xyz", return_value=identity):
                            with mock.patch.object(reaction, "reaction_product_changed", return_value=False) as product_changed:
                                with mock.patch("pyar.selection.deduplication._prefer_connected_structures") as prefer_connected:
                                    os.chdir(tmpdir)
                                    reaction.optimize_all(
                                        gamma_id="0001",
                                        orientations=[molecule],
                                        run_state=DummyRunState(),
                                        product_dir=str(Path(tmpdir) / "products"),
                                        qc_param={"gamma": 1.0},
                                    )
            finally:
                os.chdir(cwd)

        prefer_connected.assert_not_called()
        product_changed.assert_called_once()

    def test_reaction_optimize_all_promotes_smiles_only_identity_changes(self):
        class DummyReactionMolecule:
            def __init__(self):
                self.name = "orientation"
                self.energy = 0.0
                self.fragments = [[0]]

            def mol_to_xyz(self, filename):
                Path(filename).write_text("1\norientation\nH 0 0 0\n")

            def copy(self):
                return DummyReactionMolecule()

            def is_bonded(self):
                return True

        class DummyRunState:
            def current_survivor_molecules(self):
                return []

            def record_job(self, *args, **kwargs):
                return None

            def record_product(self, *args, **kwargs):
                self.recorded_product = args

        molecule = DummyReactionMolecule()
        start_identity = {"inchi": "same-inchi", "smiles": "start-smiles"}
        current_identity = {"inchi": "same-inchi", "smiles": "current-smiles"}

        with tempfile.TemporaryDirectory() as tmpdir:
            cwd = Path.cwd()
            try:
                with mock.patch.object(reaction, "optimise", return_value=True):
                    with mock.patch.object(reaction, "relax_without_afir_bias", return_value=True):
                        with mock.patch.object(
                            reaction,
                            "molecule_identity_from_xyz",
                            side_effect=[start_identity, current_identity],
                        ):
                            with mock.patch.object(reaction, "reaction_product_changed", return_value=True):
                                with mock.patch.object(reaction.shutil, "copy", return_value=None) as copy_product:
                                    with mock.patch.object(
                                        reaction.reaction_analysis,
                                        "analyse_reaction_trace",
                                        return_value={"candidate_ts_directory": "candidate_ts"},
                                    ):
                                        os.chdir(tmpdir)
                                        product_dir = Path(tmpdir) / "products"
                                        product_dir.mkdir()
                                        reaction.optimize_all(
                                            gamma_id="0001",
                                            orientations=[molecule],
                                            run_state=DummyRunState(),
                                            product_dir=str(product_dir),
                                            qc_param={"gamma": 1.0},
                                        )
                                        copy_product.assert_called_once()
            finally:
                os.chdir(cwd)

    def test_reaction_optimize_all_compares_against_original_reactants(self):
        class DummyReactionMolecule:
            def __init__(self):
                self.name = "orientation"
                self.energy = 0.0
                self.fragments = [[0]]

            def mol_to_xyz(self, filename):
                Path(filename).write_text("1\norientation\nH 0 0 0\n")

            def copy(self):
                return DummyReactionMolecule()

            def is_bonded(self):
                return True

        class DummyRunState:
            def __init__(self):
                self.data = {
                    "reactant_identity": {
                        "inchi": "reactants-inchi",
                        "smiles": "reactants-smiles",
                    }
                }
                self.recorded_jobs = []

            def current_survivor_molecules(self):
                return []

            def record_job(self, *args, **kwargs):
                self.recorded_jobs.append(args)

            def record_product(self, *args, **kwargs):
                raise AssertionError("false product should not be recorded")

        molecule = DummyReactionMolecule()
        start_identity = {"inchi": "survivor-inchi", "smiles": "survivor-smiles"}
        current_identity = {
            "inchi": "reactants-inchi",
            "smiles": "reactants-smiles",
        }
        run_state = DummyRunState()

        with tempfile.TemporaryDirectory() as tmpdir:
            cwd = Path.cwd()
            try:
                with mock.patch.object(reaction, "optimise", return_value=True):
                    with mock.patch.object(reaction, "relax_without_afir_bias", return_value=True):
                        with mock.patch.object(
                            reaction,
                            "molecule_identity_from_xyz",
                            side_effect=[start_identity, current_identity],
                        ):
                            os.chdir(tmpdir)
                            product_dir = Path(tmpdir) / "products"
                            product_dir.mkdir()
                            survivors = reaction.optimize_all(
                                gamma_id="0001",
                                orientations=[molecule],
                                run_state=run_state,
                                product_dir=str(product_dir),
                                qc_param={"gamma": 1.0},
                            )
            finally:
                os.chdir(cwd)

        self.assertEqual(len(survivors), 1)
        self.assertEqual(run_state.recorded_jobs[0][1], 1.0)
        self.assertEqual(run_state.recorded_jobs[0][2], True)


if __name__ == "__main__":
    unittest.main()
