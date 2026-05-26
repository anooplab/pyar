import importlib
import os
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest import mock

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from pyar.backend_capabilities import supported_geometry_backends
from pyar.reaction_state import ReactionStateError


class CliRealWorkflowImportsTests(unittest.TestCase):
    def test_formula_helper_uses_real_aggregate_workflow_module(self):
        cli = importlib.import_module("pyar.cli")

        self.assertEqual(cli._expand_formula_inputs("CH4"), (["C", "H"], [1, 4]))


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
            for name in (
                "pyar.Molecule",
                "pyar.aggregator",
                "pyar.reactor",
                "pyar.scan",
                "pyar.workflows",
                "pyar.workflows.aggregate",
                "pyar.workflows.solvation",
            )
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
            for attr in ("Molecule", "aggregator", "workflows", "reactor", "scan"):
                if hasattr(pyar_pkg, attr):
                    delattr(pyar_pkg, attr)

    def _install_stub_modules(self):
        def expand_formula_to_aggregate_inputs(formula):
            token = formula.strip()
            if token.lower() == "c":
                return ["C"], [1]
            if token.lower() == "h":
                return ["H"], [1]
            if token.lower() == "ch4":
                return ["C", "H"], [1, 4]
            if token.lower() == "c5h4":
                return ["C", "H"], [5, 4]
            raise ValueError(f"Unknown element symbol: {formula.strip()}")

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
        aggregate_workflow_mod = make_stub_module(
            "pyar.workflows.aggregate",
            aggregate=lambda *args, **kwargs: None,
            generate_molecule_from_formula=molecule_from_formula,
            expand_formula_to_aggregate_inputs=expand_formula_to_aggregate_inputs,
        )
        solvation_workflow_mod = make_stub_module(
            "pyar.workflows.solvation",
            solvate=lambda *args, **kwargs: None,
        )
        workflows_mod = make_stub_module(
            "pyar.workflows",
            aggregate=aggregate_workflow_mod,
            solvation=solvation_workflow_mod,
        )
        reactor_mod = make_stub_module("pyar.reactor", react=lambda *args, **kwargs: None)
        scan_mod = make_stub_module("pyar.scan", scan_distance=lambda *args, **kwargs: None)

        for name, module in (
            ("pyar.Molecule", molecule_mod),
            ("pyar.aggregator", make_stub_module("pyar.aggregator")),
            ("pyar.workflows", workflows_mod),
            ("pyar.workflows.aggregate", aggregate_workflow_mod),
            ("pyar.workflows.solvation", solvation_workflow_mod),
            ("pyar.reactor", reactor_mod),
            ("pyar.scan", scan_mod),
        ):
            sys.modules[name] = module

        pyar_pkg = sys.modules.get("pyar")
        if pyar_pkg is not None:
            pyar_pkg.Molecule = molecule_mod
            pyar_pkg.aggregator = sys.modules["pyar.aggregator"]
            pyar_pkg.workflows = workflows_mod
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

    def test_formula_expands_to_standard_aggregate_inputs(self):
        captured = {}

        def capture_aggregate(input_molecules, aggregate_sizes, hm_orientations, qc_params,
                              maximum_number_of_seeds, first_pathway, number_of_pathways, site):
            captured["names"] = [mol.name for mol in input_molecules]
            captured["sizes"] = aggregate_sizes
            captured["software"] = qc_params["software"]

        sys.modules["pyar.workflows.aggregate"].aggregate = capture_aggregate
        sys.argv = [
            "pyar-cli",
            "-a",
            "-N",
            "1",
            "--formula",
            "ch4",
            "--software",
            "aimnet_2",
        ]

        self.cli.main()

        self.assertEqual(captured["names"], ["C", "H"])
        self.assertEqual(captured["sizes"], [1, 4])
        self.assertEqual(captured["software"], "aimnet_2")

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

    def test_react_reports_restart_state_error_cleanly(self):
        def fail_resume(*args, **kwargs):
            raise ReactionStateError("Existing reaction state does not match this invocation")

        Path("a.xyz").touch()
        Path("b.xyz").touch()
        sys.modules["pyar.reactor"].react = fail_resume
        sys.argv = [
            "pyar-cli",
            "-r",
            "a.xyz",
            "b.xyz",
            "-N",
            "1",
            "-gmin",
            "100",
            "-gmax",
            "200",
            "--software",
            "xtb",
        ]

        with self.assertRaises(SystemExit) as ctx:
            self.cli.main()

        self.assertEqual(str(ctx.exception), "Existing reaction state does not match this invocation")

    def test_react_xtb_selects_geometric_afir_optimizer(self):
        captured = {}

        def capture_react(*args):
            captured["qc_params"] = args[5]

        Path("a.xyz").touch()
        Path("b.xyz").touch()
        sys.modules["pyar.reactor"].react = capture_react
        sys.argv = [
            "pyar-cli",
            "-r",
            "a.xyz",
            "b.xyz",
            "-N",
            "1",
            "-gmin",
            "100",
            "-gmax",
            "200",
            "--software",
            "xtb",
        ]

        self.cli.main()

        self.assertEqual(captured["qc_params"]["geometry_optimizer"], "geometric")
        self.assertEqual(captured["qc_params"]["opt_target"], "minimum")
        current_log = Path("pyar.log").read_text().rsplit("Run mode: react", 1)[-1]
        self.assertNotIn("ignores unsupported options: --gmin/--gmax", current_log)

    def test_react_xtb_rejects_native_optimizer_that_ignores_afir(self):
        Path("a.xyz").touch()
        Path("b.xyz").touch()
        sys.argv = [
            "pyar-cli",
            "-r",
            "a.xyz",
            "b.xyz",
            "-N",
            "1",
            "-gmin",
            "100",
            "-gmax",
            "200",
            "--software",
            "xtb",
            "--geometry-optimizer",
            "native",
        ]

        with self.assertRaises(SystemExit) as ctx:
            self.cli.main()

        self.assertEqual(
            str(ctx.exception),
            "AFIR reaction runs with "
            f"{', '.join(supported_geometry_backends())} require --geometry-optimizer geometric",
        )

    def test_react_xtb_rejects_unimplemented_transition_state_target(self):
        Path("a.xyz").touch()
        Path("b.xyz").touch()
        sys.argv = [
            "pyar-cli",
            "-r",
            "a.xyz",
            "b.xyz",
            "-N",
            "1",
            "-gmin",
            "100",
            "-gmax",
            "200",
            "--software",
            "xtb",
            "--opt-target",
            "ts",
        ]

        with self.assertRaises(SystemExit) as ctx:
            self.cli.main()

        self.assertEqual(
            str(ctx.exception),
            "Transition-state optimization is reserved for a future reaction-product workflow",
        )

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

    def test_react_subcommand_alias_dispatches_to_reactor(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cwd = os.getcwd()
            os.chdir(tmpdir)
            try:
                Path("a.xyz").write_text("1\na\nH 0.0 0.0 0.0\n")
                Path("b.xyz").write_text("1\nb\nH 0.0 0.0 0.0\n")
                sys.argv = [
                    "pyar-cli",
                    "react",
                    "a.xyz",
                    "b.xyz",
                    "-N",
                    "4",
                    "-gmin",
                    "0.1",
                    "-gmax",
                    "0.2",
                    "--software",
                    "xtb",
                ]

                with mock.patch.object(sys.modules["pyar.reactor"], "react") as react:
                    self.cli.main()
            finally:
                os.chdir(cwd)

        react.assert_called_once()

    def test_aggregate_infers_multiplicity_when_not_provided(self):
        captured = {}

        def capture_aggregate(input_molecules, aggregate_sizes, hm_orientations, qc_params,
                              maximum_number_of_seeds, first_pathway, number_of_pathways, site):
            captured["multiplicities"] = [mol.multiplicity for mol in input_molecules]
            captured["charges"] = [mol.charge for mol in input_molecules]

        sys.modules["pyar.workflows.aggregate"].aggregate = capture_aggregate
        sys.argv = [
            "pyar-cli",
            "-a",
            "C",
            "H",
            "-as",
            "1",
            "4",
            "-N",
            "1",
            "--software",
            "aimnet_2",
        ]

        self.cli.main()

        self.assertEqual(captured["charges"], [0, 0])
        self.assertEqual(captured["multiplicities"], [1, 2])

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
        self.assertIn("Backend family: none", Path("pyar.log").read_text())

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

    def test_aggregate_reports_missing_external_program_cleanly(self):
        def missing_orca(*args, **kwargs):
            raise FileNotFoundError("ORCA executable 'orca' was not found on PATH")

        sys.modules["pyar.workflows.aggregate"].aggregate = missing_orca
        sys.argv = [
            "pyar-cli",
            "-a",
            "--formula",
            "C",
            "-N",
            "8",
            "-m",
            "1",
            "--software",
            "orca",
        ]

        with self.assertRaises(SystemExit) as ctx:
            self.cli.main()

        self.assertEqual(str(ctx.exception), "ORCA executable 'orca' was not found on PATH")

    def test_aggregate_software_dependency_contract(self):
        python_only_backends = {"mlatom_aiqm1", "aiqm1_mlatom", "aimnet_2"}
        backend_messages = {
            "gaussian": "Gaussian executable 'g16' was not found on PATH",
            "mopac": "MOPAC executable 'mopac' was not found on PATH",
            "obabel": "OpenBabel executable 'obabel' was not found on PATH",
            "orca": "ORCA executable 'orca' was not found on PATH",
            "psi4": "Psi4 executable 'psi4' was not found on PATH",
            "turbomole": "Turbomole executable 'define' was not found on PATH",
            "xtb": "xTB executable 'xtb' was not found on PATH",
            "xtb_turbo": "Turbomole executable 'define' was not found on PATH",
            "xtb-aimnet2": "xTB executable 'xtb' was not found on PATH",
            "xtb-aiqm1": "xTB executable 'xtb' was not found on PATH",
        }

        def aggregate_backend_contract(input_molecules, aggregate_sizes, hm_orientations, qc_params,
                                       maximum_number_of_seeds, first_pathway, number_of_pathways, site):
            software = qc_params["software"]
            if software in python_only_backends:
                return None
            raise FileNotFoundError(backend_messages[software])

        sys.modules["pyar.workflows.aggregate"].aggregate = aggregate_backend_contract

        for software, expected in backend_messages.items():
            with self.subTest(software=software):
                sys.argv = [
                    "pyar-cli",
                    "-a",
                    "--formula",
                    "C",
                    "-N",
                    "8",
                    "-m",
                    "1",
                    "--software",
                    software,
                ]

                with self.assertRaises(SystemExit) as ctx:
                    self.cli.main()

                self.assertEqual(str(ctx.exception), expected)

        for software in sorted(python_only_backends):
            with self.subTest(software=software):
                sys.argv = [
                    "pyar-cli",
                    "-a",
                    "--formula",
                    "C",
                    "-N",
                    "8",
                    "-m",
                    "1",
                    "--software",
                    software,
                ]

                self.cli.main()

    def test_aimnet2_warns_for_unsupported_qc_flags(self):
        captured = {}

        def capture_aggregate(input_molecules, aggregate_sizes, hm_orientations, qc_params,
                              maximum_number_of_seeds, first_pathway, number_of_pathways, site):
            captured["qc_params"] = qc_params

        sys.modules["pyar.workflows.aggregate"].aggregate = capture_aggregate
        sys.argv = [
            "pyar-cli",
            "-a",
            "--formula",
            "C",
            "-N",
            "8",
            "-m",
            "1",
            "--software",
            "aimnet_2",
            "--basis",
            "def2-SVP",
            "--method",
            "B3LYP",
            "--scf-threshold",
            "tight",
            "--opt-cycles",
            "50",
            "--model",
            "requested-model",
            "--nprocs",
            "2",
        ]

        self.cli.main()
        log_text = Path("pyar.log").read_text()
        current_run_log = log_text.split("Backend 'aimnet_2' ignores unsupported options:")[-1]
        self.assertIn("Backend family: mlip", log_text)
        self.assertIn("ignores unsupported options", log_text)
        self.assertIn("--basis (basis set)", log_text)
        self.assertIn("--scf-threshold", log_text)
        self.assertIn("--opt-cycles", log_text)
        self.assertIn("--model", log_text)
        self.assertIn("--nprocs", log_text)
        self.assertIn("Ignored QC options:", log_text)
        self.assertIsNone(captured["qc_params"]["basis"])
        self.assertIsNone(captured["qc_params"]["method"])
        self.assertIsNone(captured["qc_params"]["scf_threshold"])
        self.assertIsNone(captured["qc_params"]["opt_threshold"])
        self.assertIsNone(captured["qc_params"]["opt_cycles"])
        self.assertIsNone(captured["qc_params"]["scf_cycles"])
        self.assertIsNone(captured["qc_params"]["model"])
        self.assertIsNone(captured["qc_params"]["nprocs"])
        self.assertNotIn("basis=def2-SVP", current_run_log)
        self.assertNotIn("method=B3LYP", current_run_log)
        self.assertNotIn("scf_cycles=", current_run_log)

    def test_orca_does_not_warn_for_supported_dft_flags(self):
        def aggregate_noop(*args, **kwargs):
            return None

        sys.modules["pyar.workflows.aggregate"].aggregate = aggregate_noop
        sys.argv = [
            "pyar-cli",
            "-a",
            "--formula",
            "C",
            "-N",
            "8",
            "-m",
            "1",
            "--software",
            "orca",
            "--basis",
            "def2-SVP",
            "--method",
            "B3LYP",
            "--scf-cycles",
            "200",
            "--nprocs",
            "4",
        ]

        self.cli.main()
        log_text = Path("pyar.log").read_text()
        self.assertIn("Backend family: dft_qc", log_text)
        self.assertIn("Ignored QC options: none", log_text)

    def test_xtb_uses_parallel_threads(self):
        captured = {}

        def capture_aggregate(input_molecules, aggregate_sizes, hm_orientations, qc_params,
                              maximum_number_of_seeds, first_pathway, number_of_pathways, site):
            captured["qc_params"] = qc_params

        sys.modules["pyar.workflows.aggregate"].aggregate = capture_aggregate
        sys.argv = [
            "pyar-cli",
            "-a",
            "--formula",
            "ch4",
            "-N",
            "8",
            "--software",
            "xtb",
            "--nprocs",
            "16",
        ]

        self.cli.main()
        log_text = Path("pyar.log").read_text()
        self.assertIn("xTB parallel threads: 16", log_text)
        self.assertEqual(captured["qc_params"]["nprocs"], 16)

    def test_capability_table_rejects_unwired_qc_options(self):
        _, psi4_ignored = self.cli._validate_backend_qc_options(
            "psi4",
            {"method", "basis", "custom_keywords", "nprocs"},
        )
        _, orca_ignored = self.cli._validate_backend_qc_options(
            "orca",
            {"custom_keywords", "opt_threshold", "nprocs"},
        )
        _, xtb_ignored = self.cli._validate_backend_qc_options(
            "xtb",
            {"method", "opt_threshold", "opt_cycles", "nprocs"},
        )

        self.assertEqual(psi4_ignored, ["basis", "custom_keywords", "method", "nprocs"])
        self.assertEqual(orca_ignored, ["custom_keywords", "opt_threshold"])
        self.assertEqual(xtb_ignored, ["method", "opt_cycles"])


if __name__ == "__main__":
    unittest.main()
