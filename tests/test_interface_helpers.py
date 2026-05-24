import os
import sys
import tempfile
import unittest
from types import SimpleNamespace
from unittest import mock
from pathlib import Path

from pyar import optimiser
from pyar.interface import babel
from pyar.interface import require_executable
from pyar.interface.subprocess_utils import run_command, run_output


class InterfaceHelperTests(unittest.TestCase):
    def test_require_executable_reports_missing_program_cleanly(self):
        with self.assertRaises(FileNotFoundError) as ctx:
            require_executable("__pyar_missing_program__", "MissingTool")

        self.assertEqual(
            str(ctx.exception),
            "MissingTool executable '__pyar_missing_program__' was not found on PATH",
        )

    def test_run_command_accepts_non_string_arguments(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "stdout.txt"
            exit_status = run_command(
                [sys.executable, "-c", "import sys; print(sys.argv[1])", 123],
                stdout_path=output_path,
            )

            self.assertEqual(exit_status, 0)
            self.assertEqual(output_path.read_text().strip(), "123")

    def test_run_output_accepts_pathlike_arguments(self):
        output = run_output(
            [
                Path(sys.executable),
                "-c",
                "import os, sys; print(os.path.basename(sys.argv[0]))",
            ]
        )

        self.assertTrue(output.decode("utf-8").strip())

    def test_optimiser_restores_cwd_on_missing_backend(self):
        molecule = SimpleNamespace(
            name="carbon",
            title="carbon",
            atoms_list=["C"],
            number_of_atoms=1,
            charge=0,
            multiplicity=1,
            scftype="rhf",
            coordinates=[[0.0, 0.0, 0.0]],
        )
        qc_params = {"software": "gaussian"}

        with tempfile.TemporaryDirectory() as tmpdir:
            original_cwd = os.getcwd()
            try:
                os.chdir(tmpdir)
                run_cwd = os.getcwd()
                with mock.patch.dict(os.environ, {"PATH": ""}):
                    status = optimiser.optimise(molecule, qc_params)

                self.assertIsNone(status)
                self.assertEqual(os.getcwd(), run_cwd)
            finally:
                os.chdir(original_cwd)

    def test_product_identifiers_use_modern_obabel_executable(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            xyz_path = Path(tmpdir) / "molecule.xyz"
            xyz_path.write_text("1\nH\nH 0.0 0.0 0.0\n")
            original_cwd = os.getcwd()
            try:
                os.chdir(tmpdir)
                with mock.patch.object(babel, "_obabel_executable", return_value="/usr/bin/obabel"), \
                    mock.patch.object(babel, "run_output", return_value=b"identifier") as run_output_mock:
                    self.assertEqual(babel.make_inchi_string_from_xyz(xyz_path), "identifier")
                    self.assertEqual(babel.make_smile_string_from_xyz(xyz_path), "identifier")
            finally:
                os.chdir(original_cwd)

        commands = [call.args[0] for call in run_output_mock.call_args_list]
        self.assertEqual(commands[0][0], "/usr/bin/obabel")
        self.assertEqual(commands[1][0], "/usr/bin/obabel")


if __name__ == "__main__":
    unittest.main()
