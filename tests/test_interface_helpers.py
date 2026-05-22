import os
import sys
import tempfile
import unittest
from types import SimpleNamespace
from unittest import mock
from pathlib import Path

from pyar import old_optimiser
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

    def test_old_optimiser_reports_missing_backend_and_restores_cwd(self):
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
                    with self.assertRaises(FileNotFoundError) as ctx:
                        old_optimiser.optimise(molecule, qc_params)

                self.assertEqual(
                    str(ctx.exception),
                    "Gaussian executable 'g16' was not found on PATH",
                )
                self.assertEqual(os.getcwd(), run_cwd)
            finally:
                os.chdir(original_cwd)


if __name__ == "__main__":
    unittest.main()
