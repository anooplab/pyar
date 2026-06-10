"""Tests for the solvation workflow."""

import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from pyar.solvation_state import SolvationStateError
from pyar.workflows import solvation
from pyar.workflow_results import SolvationResult


class DummyMolecule:
    def __init__(self, name, n_atoms=1):
        self.name = name
        self.title = "Title"
        self.atoms_list = ["H"] * n_atoms
        self.coordinates = [[0.0, 0.0, 0.0] for _ in range(n_atoms)]
        self.number_of_atoms = n_atoms
        self.charge = 0
        self.multiplicity = 1
        self.scftype = "rhf"
        self.fragments = [[index] for index in range(n_atoms)]
        self.energy = None

    def __len__(self):
        return len(self.atoms_list)

    def mol_to_xyz(self, filename):
        Path(filename).write_text(
            f"{len(self.atoms_list)}\n{self.name}\n"
            + "\n".join("H 0.0 0.0 0.0" for _ in self.atoms_list)
            + "\n"
        )


class SolvationWorkflowTests(unittest.TestCase):
    def test_resolve_solvation_connectivity_policy_accepts_auto_and_off(self):
        self.assertEqual(solvation._resolve_solvation_connectivity_policy("auto"), "off")
        self.assertEqual(solvation._resolve_solvation_connectivity_policy("off"), "off")

    def test_resolve_solvation_connectivity_policy_rejects_prefer_and_strict(self):
        with self.assertRaisesRegex(ValueError, "require connectivity policy 'off'"):
            solvation._resolve_solvation_connectivity_policy("prefer")
        with self.assertRaisesRegex(ValueError, "require connectivity policy 'off'"):
            solvation._resolve_solvation_connectivity_policy("strict")

    def test_solvation_writes_state_and_completes(self):
        seed = DummyMolecule("seed", n_atoms=1)
        solvent = DummyMolecule("solvent", n_atoms=1)
        step_one = DummyMolecule("step_one", n_atoms=1)
        step_two = DummyMolecule("step_two", n_atoms=1)

        with tempfile.TemporaryDirectory() as tmpdir:
            cwd = os.getcwd()
            os.chdir(tmpdir)
            try:
                with mock.patch.object(solvation, "add_one", side_effect=[ [step_one], [step_two] ]):
                    result = solvation.solvate(
                        seeds=[seed],
                        monomer=solvent,
                        aggregate_size=2,
                        hm_orientations=1,
                        qc_params={"software": None},
                        maximum_number_of_seeds=1,
                        site=None,
                    )
                state_exists = Path("solvation", "state.json").is_file()
            finally:
                os.chdir(cwd)

        self.assertIsInstance(result, SolvationResult)
        self.assertEqual(result.workflow, "solvation")
        self.assertEqual(result.status, "completed")
        self.assertTrue(result.state_path.endswith("solvation/state.json"))
        self.assertTrue(state_exists)

    def test_solvation_uses_off_connectivity_policy(self):
        seed = DummyMolecule("seed", n_atoms=1)
        solvent = DummyMolecule("solvent", n_atoms=1)
        step = DummyMolecule("step", n_atoms=1)

        with tempfile.TemporaryDirectory() as tmpdir:
            cwd = os.getcwd()
            os.chdir(tmpdir)
            try:
                with mock.patch.object(solvation, "add_one", return_value=[step]) as add_one:
                    result = solvation.solvate(
                        seeds=[seed],
                        monomer=solvent,
                        aggregate_size=1,
                        hm_orientations=1,
                        qc_params={"software": None},
                        maximum_number_of_seeds=1,
                        site=None,
                    )
                with Path("solvation", "state.json").open() as fp:
                    state = json.load(fp)
            finally:
                os.chdir(cwd)

        self.assertEqual(add_one.call_args.kwargs["connectivity_policy"], "off")
        self.assertEqual(result.metadata["connectivity_policy"], "off")
        self.assertEqual(state["request"]["connectivity_policy"], "off")

    def test_solvation_logs_relative_cycle_path(self):
        seed = DummyMolecule("seed", n_atoms=1)
        solvent = DummyMolecule("solvent", n_atoms=1)

        with tempfile.TemporaryDirectory() as tmpdir:
            cwd = os.getcwd()
            os.chdir(tmpdir)
            try:
                with mock.patch.object(solvation, "add_one", return_value=[]):
                    with self.assertLogs("pyar.workflows.aggregate", level="INFO") as captured:
                        solvation.solvate(
                            seeds=[seed],
                            monomer=solvent,
                            aggregate_size=1,
                            hm_orientations=1,
                            qc_params={"software": None},
                            maximum_number_of_seeds=1,
                            site=None,
                        )
            finally:
                os.chdir(cwd)

        self.assertTrue(any("Solvation cycle path:" in message for message in captured.output))

    def test_solvation_rejects_explicit_prefer_and_strict_connectivity_policy(self):
        seed = DummyMolecule("seed", n_atoms=1)
        solvent = DummyMolecule("solvent", n_atoms=1)

        with tempfile.TemporaryDirectory() as tmpdir:
            cwd = os.getcwd()
            os.chdir(tmpdir)
            try:
                with self.assertRaisesRegex(ValueError, "require connectivity policy 'off'"):
                    solvation.solvate(
                        seeds=[seed],
                        monomer=solvent,
                        aggregate_size=1,
                        hm_orientations=1,
                        qc_params={"software": None},
                        maximum_number_of_seeds=1,
                        site=None,
                        connectivity_policy="prefer",
                    )
                with self.assertRaisesRegex(ValueError, "require connectivity policy 'off'"):
                    solvation.solvate(
                        seeds=[seed],
                        monomer=solvent,
                        aggregate_size=1,
                        hm_orientations=1,
                        qc_params={"software": None},
                        maximum_number_of_seeds=1,
                        site=None,
                        connectivity_policy="strict",
                    )
            finally:
                os.chdir(cwd)

    def test_solvation_resumes_from_saved_state(self):
        seed = DummyMolecule("seed", n_atoms=1)
        solvent = DummyMolecule("solvent", n_atoms=1)
        step_one = DummyMolecule("step_one", n_atoms=1)
        step_two = DummyMolecule("step_two", n_atoms=1)
        first_calls = []

        def fail_on_second_cycle(*args, **kwargs):
            first_calls.append(args[0])
            if len(first_calls) == 2:
                raise RuntimeError("interrupted solvation")
            return [step_one]

        with tempfile.TemporaryDirectory() as tmpdir:
            cwd = os.getcwd()
            os.chdir(tmpdir)
            try:
                with mock.patch.object(solvation, "add_one", side_effect=fail_on_second_cycle):
                    with self.assertRaisesRegex(RuntimeError, "interrupted solvation"):
                        solvation.solvate(
                            seeds=[seed],
                            monomer=solvent,
                            aggregate_size=3,
                            hm_orientations=1,
                            qc_params={"software": None},
                            maximum_number_of_seeds=1,
                            site=None,
                        )

                with mock.patch.object(solvation, "add_one", return_value=[step_two]) as resumed_add:
                    result = solvation.solvate(
                        seeds=[seed],
                        monomer=solvent,
                        aggregate_size=3,
                        hm_orientations=1,
                        qc_params={"software": None},
                        maximum_number_of_seeds=1,
                        site=None,
                    )
            finally:
                os.chdir(cwd)

        self.assertEqual(first_calls, ["002", "003"])
        self.assertEqual(resumed_add.call_count, 2)
        self.assertEqual([call.args[0] for call in resumed_add.call_args_list], ["003", "004"])
        self.assertIsInstance(result, SolvationResult)
        self.assertEqual(result.status, "completed")
        self.assertGreaterEqual(len(result.selected_paths), 1)

    def test_solvation_refuses_changed_request(self):
        seed = DummyMolecule("seed", n_atoms=1)
        solvent = DummyMolecule("solvent", n_atoms=1)

        with tempfile.TemporaryDirectory() as tmpdir:
            cwd = os.getcwd()
            os.chdir(tmpdir)
            try:
                with mock.patch.object(solvation, "add_one", return_value=[seed]):
                    solvation.solvate(
                        seeds=[seed],
                        monomer=solvent,
                        aggregate_size=1,
                        hm_orientations=1,
                        qc_params={"software": None},
                        maximum_number_of_seeds=1,
                        site=None,
                    )

                with self.assertRaisesRegex(SolvationStateError, "does not match"):
                    solvation.solvate(
                        seeds=[seed],
                        monomer=solvent,
                        aggregate_size=2,
                        hm_orientations=1,
                        qc_params={"software": None},
                        maximum_number_of_seeds=1,
                        site=None,
                    )
            finally:
                os.chdir(cwd)


if __name__ == "__main__":
    unittest.main()
