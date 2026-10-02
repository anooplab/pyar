import io
import csv
import json
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from types import SimpleNamespace
from unittest import mock


class ClusteringScriptTests(unittest.TestCase):
    def test_cluster_mode_passes_algorithm_and_seed_limit(self):
        from pyar.scripts import clustering as clustering_script

        with tempfile.TemporaryDirectory() as tmpdir:
            path_a = Path(tmpdir, "a.xyz")
            path_b = Path(tmpdir, "b.xyz")
            path_a.write_text("1\na: 0.0\nH 0 0 0\n")
            path_b.write_text("1\nb: 1.0\nH 0 0 1\n")

            molecules = [
                SimpleNamespace(name="a", atoms_list=["H"], coordinates=[[0.0, 0.0, 0.0]], energy=0.0),
                SimpleNamespace(name="b", atoms_list=["H"], coordinates=[[0.0, 0.0, 1.0]], energy=1.0),
            ]

            with mock.patch.object(clustering_script.Molecule, "from_xyz", side_effect=molecules):
                with mock.patch.object(
                    clustering_script.clustering,
                    "choose_geometries",
                    return_value=[molecules[0]],
                ) as chooser:
                    with mock.patch(
                        "sys.argv",
                        ["pyar-clustering", str(path_a), str(path_b), "-a", "maxmin", "-n", "1"],
                    ):
                        stdout = io.StringIO()
                        with redirect_stdout(stdout):
                            clustering_script.main()

        output = stdout.getvalue()
        self.assertIn("Input pool energies:", output)
        self.assertIn("Selected pool energies:", output)
        self.assertIn("Global minimum: a (", output)
        self.assertTrue(output.strip().endswith("a.xyz"))
        self.assertEqual(chooser.call_args.kwargs["algorithm"], "maxmin")
        self.assertEqual(chooser.call_args.kwargs["maximum_number_of_seeds"], 1)

    def test_labels_mode_writes_one_row_per_structure_and_fallback_report(self):
        from pyar.scripts import clustering as clustering_script

        with tempfile.TemporaryDirectory() as tmpdir:
            xyz_files = []
            for name, distance in (("short", 0.7), ("long", 1.4)):
                path = Path(tmpdir, f"{name}.xyz")
                path.write_text(f"2\n{name}: energy={distance}\nH 0 0 0\nH {distance} 0 0\n")
                xyz_files.append(str(path))
            labels_path = Path(tmpdir, "labels.csv")
            report_path = Path(tmpdir, "report.json")
            with mock.patch(
                "sys.argv",
                ["pyar-clustering", *xyz_files, "--mode", "labels", "--feature", "distance-histogram",
                 "-a", "agglomerative", "--labels-output", str(labels_path), "--report-output", str(report_path)],
            ):
                with redirect_stdout(io.StringIO()):
                    clustering_script.main()

            with labels_path.open(newline="", encoding="utf-8") as stream:
                rows = list(csv.DictReader(stream))
            report = json.loads(report_path.read_text())

        self.assertEqual(len(rows), 2)
        self.assertEqual(report["feature_used"], "distance-histogram")
        self.assertEqual(report["algorithm_used"], "agglomerative")
        self.assertEqual(len(report["labels"]), 2)


if __name__ == "__main__":
    unittest.main()
