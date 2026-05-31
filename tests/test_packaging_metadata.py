import builtins
import importlib
import subprocess
import sys
import unittest
from importlib.metadata import version
from unittest import mock


class PackagingMetadataTests(unittest.TestCase):
    def test_distribution_version_and_import_namespace(self):
        import pyar

        self.assertEqual(pyar.__name__, "pyar")
        self.assertEqual(version("pyar-chem"), "1.1.0")

    def test_cli_help_starts(self):
        result = subprocess.run(
            [sys.executable, "-m", "pyar.cli", "--help"],
            check=False,
            capture_output=True,
            text=True,
        )

        self.assertEqual(result.returncode, 0, msg=result.stderr)
        self.assertIn("usage", result.stdout.lower())

    def test_representations_import_without_selection_extra(self):
        original_import = builtins.__import__
        original_module = sys.modules.get("pyar.representations")
        import pyar

        def block_dscribe(name, *args, **kwargs):
            if name == "dscribe" or name.startswith("dscribe."):
                raise ImportError("blocked dscribe")
            return original_import(name, *args, **kwargs)

        sys.modules.pop("pyar.representations", None)
        with mock.patch("builtins.__import__", side_effect=block_dscribe):
            representations = importlib.import_module("pyar.representations")
            with self.assertRaises(ImportError) as ctx:
                representations.mbtr_descriptor(["H"], [[0.0, 0.0, 0.0]])

        self.assertIn("selection", str(ctx.exception))
        self.assertIn("pyar-chem[selection]", str(ctx.exception))
        if original_module is not None:
            sys.modules["pyar.representations"] = original_module
            pyar.representations = original_module

    def test_clustering_import_without_selection_extra(self):
        original_import = builtins.__import__
        original_module = sys.modules.get("pyar.data_analysis.clustering")
        import pyar.data_analysis

        def block_selection_stack(name, *args, **kwargs):
            if name == "hdbscan" or name.startswith("sklearn") or name == "pandas":
                raise ImportError(f"blocked {name}")
            return original_import(name, *args, **kwargs)

        sys.modules.pop("pyar.data_analysis.clustering", None)
        with mock.patch("builtins.__import__", side_effect=block_selection_stack):
            clustering = importlib.import_module("pyar.data_analysis.clustering")
            with self.assertRaises(ImportError) as ctx:
                clustering.hdbscan_clustering([[0.0], [1.0]])

        self.assertIn("selection", str(ctx.exception))
        self.assertIn("pyar-chem[selection]", str(ctx.exception))
        if original_module is not None:
            sys.modules["pyar.data_analysis.clustering"] = original_module
            pyar.data_analysis.clustering = original_module


if __name__ == "__main__":
    unittest.main()
