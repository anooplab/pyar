"""Tests for structured reaction restart state."""

import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from pyar.core.molecule import Molecule
from pyar.state.reaction_v2 import ReactionRunStateV2 as ReactionRunState
from pyar.state.reaction import ReactionStateError
from pyar.workflows import reaction as reaction_workflow


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
        self.sampling = {
            "direction_method": "fibonacci",
            "rotation_method": "halton",
            "number_of_orientations": 1,
            "sequence_offset": 0,
            "oversample_factor": 8,
            "use_angles": True,
        }

    def test_create_writes_json_and_roundtrips_pending_geometry_metadata(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            Path(tmpdir, "reaction").mkdir()
            state = ReactionRunState.create(
                tmpdir,
                self.request,
                [self.orientation],
                (self.reactant_a, self.reactant_b),
                sampling=self.sampling,
            )
            loaded_state = ReactionRunState.load(tmpdir, self.request)
            molecule = loaded_state.pending_molecules()[0]

            with state.state_file.open() as fp:
                raw_state = json.load(fp)

        self.assertEqual(raw_state["workflow_state_schema_version"], 3)
        self.assertEqual(raw_state["gamma_schedule"], [100.0, 200.0])
        self.assertEqual(raw_state["sampling"], self.sampling)
        self.assertEqual(molecule.fragments, [[0], [1]])
        self.assertEqual(molecule.multiplicity, 2)
        np.testing.assert_allclose(molecule.coordinates, self.orientation.coordinates)

    def test_backend_failure_preserves_pending_work_and_reason(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            state = ReactionRunState.create(
                tmpdir, self.request, [self.orientation],
                (self.reactant_a, self.reactant_b),
            )
            state.fail("xTB rejected the input")
            saved = json.loads(state.state_file.read_text())
            self.assertEqual(saved["status"], "failed_backend")
            self.assertEqual(saved["failure"], "xTB rejected the input")
            self.assertEqual(len(saved["pending_orientations"]), 1)

    def test_load_rejects_changed_request(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            Path(tmpdir, "reaction").mkdir()
            ReactionRunState.create(
                tmpdir,
                self.request,
                [self.orientation],
                (self.reactant_a, self.reactant_b),
                sampling=self.sampling,
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
                sampling=self.sampling,
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
                sampling=self.sampling,
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
                sampling=self.sampling,
            )
            product_path = products / "product_0.xyz"
            self.orientation.mol_to_xyz(str(product_path))
            state.record_product("product_0", 100.0, "inchi-0", "smiles-0", product_path)
            resumed = ReactionRunState.load(tmpdir, self.request)

        self.assertEqual(
            resumed.saved_product_identities(),
            {"product_0": ("inchi-0", "smiles-0")},
        )

    def test_record_product_persists_trace_summary_metadata(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            products = Path(tmpdir, "reaction", "products")
            products.mkdir(parents=True)
            state = ReactionRunState.create(
                tmpdir,
                self.request,
                [self.orientation],
                (self.reactant_a, self.reactant_b),
                sampling=self.sampling,
            )
            product_path = products / "product_0.xyz"
            self.orientation.mol_to_xyz(str(product_path))
            state.record_product(
                "product_0",
                100.0,
                "inchi-0",
                "smiles-0",
                product_path,
                trace_summary={
                    "path_summary_csv": "reaction/path_summary.csv",
                    "candidate_ts_directory": "reaction/candidate_ts",
                },
            )
            loaded = ReactionRunState.load(tmpdir, self.request)

        self.assertEqual(
            loaded.data["products"][0]["trace_summary"]["candidate_ts_directory"],
            "discovery/products/product_001/sources/reaction/candidate_ts",
        )

    def test_missing_geometry_snapshot_reports_restart_state_error(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            Path(tmpdir, "reaction").mkdir()
            state = ReactionRunState.create(
                tmpdir,
                self.request,
                [self.orientation],
                (self.reactant_a, self.reactant_b),
                sampling=self.sampling,
            )
            pending_path = state.reaction_directory / state.data["pending_orientations"][0]["path"]
            pending_path.unlink()

            with self.assertRaisesRegex(ReactionStateError, "Invalid reaction snapshot"):
                state.pending_molecules()






def test_reaction_product_decisions_use_explicit_run_registry():
    from unittest.mock import patch

    first_run = {"first": {"inchi": "inchi-1", "smiles": "C"}}
    second_run = {"second": {"inchi": "inchi-2", "smiles": "N"}}
    incoming = {"inchi": "inchi-1", "smiles": "C"}
    with patch.object(
        reaction_workflow, "same_molecular_identity",
        side_effect=lambda left, right: left["inchi"] == right["inchi"],
    ):
        assert reaction_workflow._is_known_product(incoming, first_run)
        assert not reaction_workflow._is_known_product(incoming, second_run)
    assert reaction_workflow._should_continue_after_product([], first_run)
    assert reaction_workflow._should_continue_after_product([], second_run)
    assert not reaction_workflow._should_continue_after_product([object()], first_run)


if __name__ == "__main__":
    unittest.main()
