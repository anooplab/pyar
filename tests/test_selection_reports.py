"""Tests for selection reporting helpers."""

import io
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

from pyar.selection import reports


class SelectionReportsTests(unittest.TestCase):
    def test_energy_table_output_is_stable(self):
        molecules = [
            SimpleNamespace(name="low", energy=0.0),
            SimpleNamespace(name="high", energy=1.0),
        ]

        stream = io.StringIO()
        reports.print_energy_table(molecules, stream=stream, title="Selection energies:")

        output = stream.getvalue()
        self.assertIn("Selection energies:", output)
        self.assertIn("R. E. (kcal/mol)", output)
        self.assertIn("Global minimum: low", output)

    def test_read_energy_from_xyz_file_uses_trailing_value(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            xyz_file = Path(tmpdir) / "example.xyz"
            xyz_file.write_text("1\ntrial orientation 002_002_ag_a_001_b_004: -40.5243989267341\nH 0 0 0\n")

            energy = reports.read_energy_from_xyz_file(str(xyz_file))

        self.assertAlmostEqual(energy, -40.5243989267341)


if __name__ == "__main__":
    unittest.main()
