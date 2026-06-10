import io
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from unittest import mock


class EnergyTableScriptTests(unittest.TestCase):
    def test_energy_table_script_reports_relative_energies(self):
        from pyar.scripts import energy_table as energy_table_script

        with tempfile.TemporaryDirectory() as tmpdir:
            path_a = Path(tmpdir, "a.xyz")
            path_b = Path(tmpdir, "b.xyz")
            path_a.write_text("1\na: -10.0\nH 0 0 0\n")
            path_b.write_text("1\nb: -9.5\nH 0 0 1\n")

            with mock.patch(
                "sys.argv",
                ["pyar-energy-table", str(path_a), str(path_b)],
            ):
                stdout = io.StringIO()
                with redirect_stdout(stdout):
                    energy_table_script.main()

        output = stdout.getvalue()
        self.assertIn("Relative energy table:", output)
        self.assertIn("Global minimum: a (", output)
        self.assertIn("b", output)
        self.assertIn("313.75", output)


if __name__ == "__main__":
    unittest.main()
