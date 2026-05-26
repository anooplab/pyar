"""Tests for the solvation workflow."""

import os
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from pyar.solvation_state import SolvationStateError
from pyar.workflows import solvation


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

        self.assertIsNone(result)
        self.assertTrue(state_exists)

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
        self.assertIsNone(result)

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
