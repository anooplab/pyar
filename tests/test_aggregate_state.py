"""Tests for structured aggregation restart state."""

import tempfile
import unittest

from pyar.aggregate_state import AggregateRunState, AggregateStateError


class AggregateRunStateTests(unittest.TestCase):
    def setUp(self):
        self.request = {
            "aggregate_sizes": [1, 2],
            "orientations": 8,
            "backend_parameters": {"software": "xtb"},
            "maximum_number_of_seeds": 4,
            "first_pathway": 0,
            "number_of_pathways": 2,
            "site": None,
            "fragments": [{"atoms": ["C"], "coordinates": [[0.0, 0.0, 0.0]]}],
        }

    def test_create_and_load_roundtrips_remaining_pathways(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            state = AggregateRunState.create(tmpdir, self.request, ["abb", "bab"])
            state.complete_pathway("abb", ["ag_a_001_b_001_000/selected/result_one.xyz"])
            loaded = AggregateRunState.load(tmpdir, self.request)

        self.assertEqual(loaded.completed_pathway_count(), 1)
        self.assertEqual(loaded.remaining_pathway_labels(), ["bab"])

    def test_load_rejects_changed_request(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            AggregateRunState.create(tmpdir, self.request, ["abb"])
            changed = {**self.request, "maximum_number_of_seeds": 5}

            with self.assertRaisesRegex(AggregateStateError, "does not match"):
                AggregateRunState.load(tmpdir, changed)

    def test_load_rejects_out_of_sequence_completed_pathways(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            state = AggregateRunState.create(tmpdir, self.request, ["abb", "bab"])
            state.data["completed_pathways"] = [
                {"index": 0, "label": "bab", "status": "completed", "selected_results": []}
            ]
            state.save()

            with self.assertRaisesRegex(AggregateStateError, "out of sequence"):
                AggregateRunState.load(tmpdir, self.request)

    def test_finish_rejects_pending_pathways(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            state = AggregateRunState.create(tmpdir, self.request, ["abb"])

            with self.assertRaisesRegex(AggregateStateError, "remain unfinished"):
                state.finish([])

    def test_finish_retains_results_and_prevents_completed_resume(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            state = AggregateRunState.create(tmpdir, self.request, ["abb"])
            state.complete_pathway("abb", [])
            state.finish(["selected/stoichiometry_CH2/result_final.xyz"])

            self.assertEqual(
                state.data["final_selected_results"],
                ["selected/stoichiometry_CH2/result_final.xyz"],
            )
            with self.assertRaisesRegex(AggregateStateError, "already 'completed'"):
                AggregateRunState.load(tmpdir, self.request)


if __name__ == "__main__":
    unittest.main()
