"""Tests for structured reaction restart state."""

import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from pyar.Molecule import Molecule
from pyar.reaction_state import ReactionRunState, ReactionStateError


class ReactionRunStateTests(unittest.TestCase):
    def setUp(self):
        self.reactant_a = Molecule(["C"], np.array([[0.0, 0.0, 0.0]]), name="a")
        self.reactant_b = Molecule(["H"], np.array([[0.0, 0.0, 1.0]]), name="b")
        self.orientation = Molecule(
            ["C", "H"],
            np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]]),
            name="orientation_0",
            fragments=[[0], [1]],
            charge=0,
            multiplicity=2,
        )
        self.request = {
            "gamma_schedule": [100.0, 200.0],
            "orientations": 1,
            "backend_parameters": {"software": "xtb_turbo"},
            "site": None,
            "proximity_factor": 2.3,
            "reactants": [],
        }

    def test_create_writes_json_and_roundtrips_pending_geometry_metadata(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            Path(tmpdir, "reaction").mkdir()
            state = ReactionRunState.create(
                tmpdir,
                self.request,
                [self.orientation],
                (self.reactant_a, self.reactant_b),
            )
            loaded_state = ReactionRunState.load(tmpdir, self.request)
            molecule = loaded_state.pending_molecules()[0]

            with state.state_file.open() as fp:
                raw_state = json.load(fp)

        self.assertEqual(raw_state["version"], 2)
        self.assertEqual(raw_state["gamma_schedule"], [100.0, 200.0])
        self.assertEqual(molecule.fragments, [[0], [1]])
        self.assertEqual(molecule.multiplicity, 2)
        np.testing.assert_allclose(molecule.coordinates, self.orientation.coordinates)

    def test_load_rejects_changed_request(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            Path(tmpdir, "reaction").mkdir()
            ReactionRunState.create(
                tmpdir,
                self.request,
                [self.orientation],
                (self.reactant_a, self.reactant_b),
            )
            changed = {**self.request, "gamma_schedule": [100.0, 300.0]}

            with self.assertRaisesRegex(ReactionStateError, "does not match"):
                ReactionRunState.load(tmpdir, changed)

    def test_record_job_and_complete_cycle_are_persisted(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            Path(tmpdir, "reaction").mkdir()
            state = ReactionRunState.create(
                tmpdir,
                self.request,
                [self.orientation],
                (self.reactant_a, self.reactant_b),
            )
            state.record_job("job_0", 100.0, "success", [], [self.orientation])
            state.complete_cycle(100.0, [self.orientation])
            loaded = ReactionRunState.load(tmpdir, self.request)

        self.assertEqual(loaded.data["current_cycle"], 1)
        self.assertEqual(loaded.data["completed_jobs"][0]["job_name"], "job_0")
        self.assertEqual(loaded.remaining_gamma_schedule(), [200.0])

    def test_current_cycle_survivors_are_available_after_restart(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            Path(tmpdir, "reaction").mkdir()
            state = ReactionRunState.create(
                tmpdir,
                self.request,
                [self.orientation],
                (self.reactant_a, self.reactant_b),
            )
            state.record_job("job_0", 100.0, "success", [], [self.orientation])
            resumed = ReactionRunState.load(tmpdir, self.request)
            survivors = resumed.current_survivor_molecules()

        self.assertEqual(len(survivors), 1)
        self.assertEqual(survivors[0].fragments, [[0], [1]])

    def test_product_identity_is_persisted_immediately_for_restart(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            products = Path(tmpdir, "reaction", "products")
            products.mkdir(parents=True)
            state = ReactionRunState.create(
                tmpdir,
                self.request,
                [self.orientation],
                (self.reactant_a, self.reactant_b),
            )
            product_path = products / "product_0.xyz"
            self.orientation.mol_to_xyz(str(product_path))
            state.record_product("product_0", 100.0, "inchi-0", "smiles-0", product_path)
            resumed = ReactionRunState.load(tmpdir, self.request)

        self.assertEqual(
            resumed.saved_product_identities(),
            {"product_0": ("inchi-0", "smiles-0")},
        )

    def test_missing_geometry_snapshot_reports_restart_state_error(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            Path(tmpdir, "reaction").mkdir()
            state = ReactionRunState.create(
                tmpdir,
                self.request,
                [self.orientation],
                (self.reactant_a, self.reactant_b),
            )
            pending_path = state.reaction_directory / state.data["pending_orientations"][0]["path"]
            pending_path.unlink()

            with self.assertLogs("pyar.molecule", level="ERROR"):
                with self.assertRaisesRegex(ReactionStateError, "Could not restore reaction geometry snapshot"):
                    state.pending_molecules()

    def test_legacy_import_maps_unambiguous_remaining_schedule(self):
        legacy = {"0200": [self.orientation]}
        with tempfile.TemporaryDirectory() as tmpdir:
            Path(tmpdir, "reaction").mkdir()
            state = ReactionRunState.migrate_legacy(tmpdir, legacy, self.request)

        self.assertEqual(state.data["current_cycle"], 1)
        self.assertEqual(state.remaining_gamma_schedule(), [200.0])

    def test_legacy_import_rejects_collapsed_fractional_schedule(self):
        request = {**self.request, "gamma_schedule": [0.1, 0.2]}
        with tempfile.TemporaryDirectory() as tmpdir:
            Path(tmpdir, "reaction").mkdir()
            with self.assertRaisesRegex(ReactionStateError, "cannot be safely resumed"):
                ReactionRunState.migrate_legacy(tmpdir, {"0000": [self.orientation]}, request)


if __name__ == "__main__":
    unittest.main()
