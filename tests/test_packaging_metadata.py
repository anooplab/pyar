import builtins
import importlib
import subprocess
import sys
import unittest
import zipfile
from importlib.metadata import version
from pathlib import Path
from unittest import mock


class PackagingMetadataTests(unittest.TestCase):
    def test_distribution_version_and_import_namespace(self):
        import pyar

        self.assertEqual(pyar.__name__, "pyar")
        self.assertEqual(version("pyar-chem"), "1.1.1")

    def test_cli_help_starts(self):
        result = subprocess.run(
            [sys.executable, "-m", "pyar.cli", "--help"],
            check=False,
            capture_output=True,
            text=True,
        )

        self.assertEqual(result.returncode, 0, msg=result.stderr)
        self.assertIn("usage", result.stdout.lower())

    def test_built_wheel_excludes_large_model_assets(self):
        wheels = sorted(Path("dist").glob("pyar_chem-*.whl"))
        if not wheels:
            self.skipTest("built wheel not available")

        with zipfile.ZipFile(wheels[-1]) as archive:
            names = set(archive.namelist())

        blocked_exact = {
            "pyar/mlatom/MLatomF",
            "pyar/mlatom/cs.so",
            "pyar/mlatom/ref.json",
        }
        blocked_patterns = (
            ("pyar/AIMNet2/models/", ".jpt"),
            ("pyar/mlatom/aiqm1_model/", ".pt"),
        )
        offenders = sorted(
            name
            for name in names
            if name in blocked_exact
            or any(
                name.startswith(prefix) and name.endswith(suffix)
                for prefix, suffix in blocked_patterns
            )
        )
        self.assertEqual(offenders, [])

    def test_aimnet2_missing_asset_error_describes_external_model_policy(self):
        from pyar.backends import aimnet2_assets

        with mock.patch.object(
            aimnet2_assets,
            "missing_aimnet2_assets",
            return_value=["missing-model.jpt"],
        ):
            with self.assertRaisesRegex(FileNotFoundError, "not bundled"):
                aimnet2_assets.validate_aimnet2_runtime_assets()

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
        original_module = sys.modules.get("pyar.selection.clustering")
        import pyar.selection

        def block_selection_stack(name, *args, **kwargs):
            if name == "hdbscan" or name.startswith("sklearn") or name == "pandas":
                raise ImportError(f"blocked {name}")
            return original_import(name, *args, **kwargs)

        sys.modules.pop("pyar.selection.clustering", None)
        with mock.patch("builtins.__import__", side_effect=block_selection_stack):
            clustering = importlib.import_module("pyar.selection.clustering")
            with self.assertRaises(ImportError) as ctx:
                clustering.hdbscan_clustering([[0.0], [1.0]])

        self.assertIn("selection", str(ctx.exception))
        self.assertIn("pyar-chem[selection]", str(ctx.exception))
        if original_module is not None:
            sys.modules["pyar.selection.clustering"] = original_module
            pyar.selection.clustering = original_module


if __name__ == "__main__":
    unittest.main()
