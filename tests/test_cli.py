import importlib
import os
import sys
import tempfile
import types
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


def make_stub_module(name, **attrs):
    module = types.ModuleType(name)
    for key, value in attrs.items():
        setattr(module, key, value)
    return module


class CliSmokeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._original_cwd = os.getcwd()
        cls._tmpdir = tempfile.TemporaryDirectory()
        os.chdir(cls._tmpdir.name)
        cls.cli = importlib.import_module("pyar.cli")

    @classmethod
    def tearDownClass(cls):
        os.chdir(cls._original_cwd)
        cls._tmpdir.cleanup()

    def setUp(self):
        self._original_argv = sys.argv[:]
        self._original_modules = {
            name: sys.modules.get(name)
            for name in ("pyar.Molecule", "pyar.aggregator", "pyar.reactor", "pyar.scan")
        }
        self._install_stub_modules()

    def tearDown(self):
        sys.argv = self._original_argv
        for name, module in self._original_modules.items():
            if module is None:
                sys.modules.pop(name, None)
            else:
                sys.modules[name] = module

        pyar_pkg = sys.modules.get("pyar")
        if pyar_pkg is not None:
            for attr in ("Molecule", "aggregator", "reactor", "scan"):
                if hasattr(pyar_pkg, attr):
                    delattr(pyar_pkg, attr)

    def _install_stub_modules(self):
        def molecule_from_formula(formula):
            symbol = formula.strip().capitalize()
            if symbol not in {"C", "H"}:
                raise ValueError(f"Unknown element symbol: {formula.strip()}")
            return types.SimpleNamespace(
                atomic_number=[6] if symbol == "C" else [1],
                charge=0,
                multiplicity=1,
                name=formula,
                number_of_atoms=1,
            )

        class FakeMolecule:
            @classmethod
            def from_xyz(cls, path):
                return types.SimpleNamespace(
                    atomic_number=[1, 1],
                    charge=0,
                    multiplicity=1,
                    name=Path(path).stem,
                    number_of_atoms=2,
                )

        molecule_mod = make_stub_module("pyar.Molecule", Molecule=FakeMolecule)
        aggregator_mod = make_stub_module(
            "pyar.aggregator",
            aggregate=lambda *args, **kwargs: None,
            solvate=lambda *args, **kwargs: None,
            generate_molecule_from_formula=molecule_from_formula,
        )
        reactor_mod = make_stub_module("pyar.reactor", react=lambda *args, **kwargs: None)
        scan_mod = make_stub_module("pyar.scan", scan_distance=lambda *args, **kwargs: None)

        for name, module in (
            ("pyar.Molecule", molecule_mod),
            ("pyar.aggregator", aggregator_mod),
            ("pyar.reactor", reactor_mod),
            ("pyar.scan", scan_mod),
        ):
            sys.modules[name] = module

        pyar_pkg = sys.modules.get("pyar")
        if pyar_pkg is not None:
            pyar_pkg.Molecule = molecule_mod
            pyar_pkg.aggregator = aggregator_mod
            pyar_pkg.reactor = reactor_mod
            pyar_pkg.scan = scan_mod

    def test_formula_rejects_xyz_inputs(self):
        sys.argv = [
            "pyar-cli",
            "-a",
            "-N",
            "8",
            "--formula",
            "C5H4",
            "seed.xyz",
        ]

        with self.assertRaises(SystemExit) as ctx:
            self.cli.main()

        self.assertEqual(str(ctx.exception), "Do not provide XYZ input files when using --formula")

    def test_formula_accepts_charge_and_multiplicity(self):
        sys.argv = [
            "pyar-cli",
            "-a",
            "-N",
            "8",
            "--formula",
            "C",
            "-c",
            "0",
            "-m",
            "1",
            "--software",
            "aiqm1_mlatom",
        ]

        self.cli.main()

    def test_react_requires_exactly_two_xyz_inputs(self):
        sys.argv = [
            "pyar-cli",
            "-r",
            "-N",
            "8",
            "seed.xyz",
        ]

        with self.assertRaises(SystemExit) as ctx:
            self.cli.main()

        self.assertEqual(str(ctx.exception), "Reactor requires exactly two XYZ input files")

    def test_aggregate_accepts_element_symbols(self):
        sys.argv = [
            "pyar-cli",
            "-a",
            "C",
            "H",
            "-as",
            "1",
            "4",
            "-N",
            "8",
        ]

        self.cli.main()

    def test_aggregate_without_software_logs_geometry_only_mode(self):
        sys.argv = [
            "pyar-cli",
            "-a",
            "C",
            "H",
            "-as",
            "1",
            "4",
            "-N",
            "8",
        ]

        self.cli.main()

        self.assertIn(
            "No --software specified: aggregate mode will generate trial "
            "geometries only; no quantum-chemistry optimization will be run.",
            Path("pyar.log").read_text(),
        )

    def test_aggregate_rejects_unknown_element_symbol(self):
        sys.argv = [
            "pyar-cli",
            "-a",
            "Xx",
            "-N",
            "8",
        ]

        with self.assertRaises(SystemExit) as ctx:
            self.cli.main()

        self.assertEqual(str(ctx.exception), "Unknown element symbol: Xx")


if __name__ == "__main__":
    unittest.main()
