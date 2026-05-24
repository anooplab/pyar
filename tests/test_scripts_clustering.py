import io
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
        self.assertIn("Global minimum: a", output)
        self.assertTrue(output.strip().endswith("a.xyz"))
        self.assertEqual(chooser.call_args.kwargs["algorithm"], "maxmin")
        self.assertEqual(chooser.call_args.kwargs["maximum_number_of_seeds"], 1)


if __name__ == "__main__":
    unittest.main()
