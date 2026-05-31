"""Tests for structured solvation restart state."""

import tempfile
import unittest

from pyar.core.molecule import Molecule
from pyar.solvation_state import SolvationRunState, SolvationStateError


class SolvationRunStateTests(unittest.TestCase):
    def setUp(self):
        self.seed = Molecule(["C"], [[0.0, 0.0, 0.0]], name="seed")
        self.monomer = Molecule(["H"], [[0.0, 0.0, 1.0]], name="solvent")
        self.request = {
            "aggregate_size": 3,
            "orientations": 8,
            "backend_parameters": {"software": "xtb"},
            "maximum_number_of_seeds": 4,
            "site": None,
            "seeds": [
                {
                    "name": "seed",
                    "atoms": ["C"],
                    "coordinates": [[0.0, 0.0, 0.0]],
                    "charge": 0,
                    "multiplicity": 1,
                    "scftype": "rhf",
                    "fragment_definition": [[0]],
                }
            ],
            "monomer": {
                "name": "solvent",
                "atoms": ["H"],
                "coordinates": [[0.0, 0.0, 1.0]],
                "charge": 0,
                "multiplicity": 1,
                "scftype": "rhf",
                "fragment_definition": [[0]],
            },
        }
        self.sampling = {
            "direction_method": "fibonacci",
            "rotation_method": "halton",
            "number_of_orientations": 8,
            "sequence_offset": 0,
            "oversample_factor": 8,
            "use_angles": False,
        }

    def test_create_and_load_roundtrips_current_seeds(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            state = SolvationRunState.create(
                tmpdir,
                self.request,
                [self.seed],
                sampling=self.sampling,
            )
            loaded = SolvationRunState.load(tmpdir, self.request)
            self.assertEqual(state.data["next_cycle"], 2)
            self.assertEqual([mol.name for mol in loaded.current_seed_molecules()], ["seed"])
            self.assertEqual(loaded.data["sampling"], self.sampling)

    def test_load_rejects_changed_request(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            SolvationRunState.create(tmpdir, self.request, [self.seed], sampling=self.sampling)
            changed = {**self.request, "aggregate_size": 4}

            with self.assertRaisesRegex(SolvationStateError, "does not match"):
                SolvationRunState.load(tmpdir, changed)

    def test_complete_cycle_updates_next_cycle_and_current_seeds(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            state = SolvationRunState.create(
                tmpdir,
                self.request,
                [self.seed],
                sampling=self.sampling,
            )
            next_seed = Molecule(
                ["C", "H"],
                [[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
                name="seed2",
                fragments=[[0], [1]],
            )
            state.complete_cycle(2, [next_seed])
            loaded = SolvationRunState.load(tmpdir, self.request)
            self.assertEqual(state.data["next_cycle"], 3)
            self.assertEqual([mol.name for mol in loaded.current_seed_molecules()], ["seed2"])

    def test_finish_rejects_pending_cycles(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            state = SolvationRunState.create(
                tmpdir,
                self.request,
                [self.seed],
                sampling=self.sampling,
            )
            with self.assertRaisesRegex(SolvationStateError, "remain unfinished"):
                state.finish()

    def test_load_rejects_completed_resume(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            request = {**self.request, "aggregate_size": 1}
            state = SolvationRunState.create(tmpdir, request, [self.seed], sampling=self.sampling)
            state.complete_cycle(2, [self.seed])
            state.finish()

            with self.assertRaisesRegex(SolvationStateError, "already 'completed'"):
                SolvationRunState.load(tmpdir, request)


if __name__ == "__main__":
    unittest.main()
