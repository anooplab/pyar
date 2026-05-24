import os
import tempfile
import unittest
from types import SimpleNamespace
from unittest import mock

import numpy as np

from pyar.interface.xtb_utils import build_xtb_command, xtb_parallel_args


class XtbInterfaceTests(unittest.TestCase):
    def setUp(self):
        self.molecule = SimpleNamespace(
            name="job",
            title="job",
            atoms_list=["C", "H", "H", "H", "H"],
            number_of_atoms=5,
            charge=0,
            multiplicity=1,
            scftype="rhf",
            coordinates=[
                [0.0, 0.0, 0.0],
                [1.0, 1.0, 1.0],
                [1.0, -1.0, -1.0],
                [-1.0, 1.0, -1.0],
                [-1.0, -1.0, 1.0],
            ],
            fragments=[[0], [1, 2, 3, 4]],
        )

    def test_xtb_parallel_args_omits_invalid_values(self):
        self.assertEqual(xtb_parallel_args({"nprocs": 4}), ["--parallel", "4"])
        self.assertEqual(xtb_parallel_args({"nprocs": 0}), [])
        self.assertEqual(xtb_parallel_args({"nprocs": "not-an-int"}), [])

    def test_build_xtb_command_includes_common_qc_settings(self):
        command = build_xtb_command(
            "xtb",
            "input.xyz",
            {"nprocs": 8, "charge": -1, "multiplicity": 2, "scftype": "uhf"},
            opt_threshold="tight",
        )

        self.assertEqual(
            command,
            ["xtb", "input.xyz", "--parallel", "8", "-opt", "tight", "-chrg", "-1", "-uhf", "2"],
        )

    def test_xtb_wrapper_uses_parallel_threads(self):
        from pyar.interface import xtb

        with tempfile.TemporaryDirectory():
            with mock.patch.object(xtb, "require_executable", return_value="xtb"):
                runner = xtb.Xtb(self.molecule, {"opt_threshold": "normal", "nprocs": 16})

        self.assertIn("--parallel", runner.cmd)
        self.assertIn("16", runner.cmd)
        self.assertIn("-opt", runner.cmd)

    def test_xtb_turbo_wrapper_uses_parallel_threads(self):
        from pyar.interface import xtb_turbo

        with tempfile.TemporaryDirectory():
            with mock.patch.object(xtb_turbo, "require_executable", side_effect=["define", "xtb"]):
                runner = xtb_turbo.XtbTurbo(self.molecule, {"nprocs": 12})

        self.assertIn("--parallel", runner.egrad_program)
        self.assertIn("12", runner.egrad_program)
        self.assertIn("-grad", runner.egrad_program)

    def test_xtb_turbo_uses_requested_gamma(self):
        from pyar.interface import xtb_turbo

        molecule = SimpleNamespace(**self.molecule.__dict__)
        molecule.coordinates = np.asarray(self.molecule.coordinates, dtype=float)

        with tempfile.TemporaryDirectory() as tmpdir:
            cwd = os.getcwd()
            os.chdir(tmpdir)
            try:
                with mock.patch.object(xtb_turbo, "require_executable", side_effect=["define", "xtb"]):
                    runner = xtb_turbo.XtbTurbo(molecule, {"nprocs": 12, "gamma": "37.5"})

                with mock.patch.object(xtb_turbo.turbomole, "make_coord"), \
                    mock.patch.object(xtb_turbo.turbomole, "prepare_control"), \
                    mock.patch.object(xtb_turbo.turbomole, "get_coords", return_value=np.zeros((5, 3))), \
                    mock.patch.object(xtb_turbo.turbomole, "rewrite_turbomole_energy_and_gradient_files"), \
                    mock.patch.object(xtb_turbo.turbomole, "update_coord", return_value=False), \
                    mock.patch.object(runner, "calculate_energy_gradient", return_value=(True, [], 1.0, np.zeros((5, 3)))), \
                    mock.patch.object(xtb_turbo.restraints, "isotropic", return_value=(0.0, np.zeros((5, 3)))) as isotropic:
                    status = runner.optimize()
            finally:
                os.chdir(cwd)

        self.assertEqual(status, "UpdateFailed")
        self.assertEqual(runner.gamma, 37.5)
        isotropic.assert_called_once()
        self.assertEqual(isotropic.call_args.args[-1], 37.5)

    def test_xtb_aiqm1_wrapper_uses_parallel_threads(self):
        from pyar.interface import xtb_aiqm1

        with tempfile.TemporaryDirectory():
            with mock.patch.object(xtb_aiqm1, "require_executable", return_value="xtb"):
                runner = xtb_aiqm1.XtbAIQM1(self.molecule, {"opt_threshold": "normal", "nprocs": 4})

        self.assertIn("--parallel", runner.xtb_cmd)
        self.assertIn("4", runner.xtb_cmd)
        self.assertIn("-opt", runner.xtb_cmd)

    def test_xtb_aimnet2_wrapper_uses_parallel_threads(self):
        from pyar.interface import xtb_aimnet2

        with tempfile.TemporaryDirectory():
            with mock.patch.object(xtb_aimnet2, "require_executable", return_value="xtb"):
                runner = xtb_aimnet2.XtbAimnet2(self.molecule, {"opt_threshold": "normal", "nprocs": 6})

        self.assertIn("--parallel", runner.xtb_cmd)
        self.assertIn("6", runner.xtb_cmd)
        self.assertIn("-opt", runner.xtb_cmd)

    def test_turbomole_wrapper_stores_gamma_from_parameters(self):
        from pyar.interface import turbomole

        molecule = SimpleNamespace(**self.molecule.__dict__)
        molecule.coordinates = np.asarray(self.molecule.coordinates, dtype=float)

        with mock.patch.object(turbomole, "require_executable", return_value="define"):
            runner = turbomole.Turbomole(
                molecule,
                {"basis": "def2-SVP", "method": "bp86", "gamma": "22.0"},
            )

        self.assertEqual(runner.gamma, 22.0)


if __name__ == "__main__":
    unittest.main()
