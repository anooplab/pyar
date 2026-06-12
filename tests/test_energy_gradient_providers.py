import unittest
import sys
from types import ModuleType
from types import SimpleNamespace
from pathlib import Path
from unittest import mock

import numpy as np
from ase import Atoms
from ase.units import Hartree

from pyar.backend_capabilities import (
    BackendCapabilities,
    backend_supports_geometry_optimization,
    get_backend_capabilities,
    register_backend_capabilities,
)
from pyar.energy_gradient_providers import (
    ENERGY_GRADIENT_PROVIDERS,
    EnergyGradientResult,
    get_energy_gradient_provider,
    register_energy_gradient_provider,
)


class EnergyGradientProviderTests(unittest.TestCase):
    def test_energy_gradient_result_rejects_nonfinite_energy(self):
        with self.assertRaisesRegex(ValueError, "finite scalar"):
            EnergyGradientResult(np.nan, np.zeros((1, 3)))

    def test_energy_gradient_result_rejects_invalid_gradient_shape(self):
        with self.assertRaisesRegex(ValueError, "shape \\(natoms, 3\\)"):
            EnergyGradientResult(1.0, np.zeros((3,)))

    def test_energy_gradient_result_rejects_nonfinite_gradient_values(self):
        gradient = np.zeros((2, 3))
        gradient[1, 2] = np.inf
        with self.assertRaisesRegex(ValueError, "finite values"):
            EnergyGradientResult(1.0, gradient)

    def test_registered_geometric_providers_are_accessible(self):
        self.assertIn("xtb", ENERGY_GRADIENT_PROVIDERS)
        self.assertIn("aimnet_2", ENERGY_GRADIENT_PROVIDERS)
        self.assertIn("gaussian", ENERGY_GRADIENT_PROVIDERS)
        self.assertIn("orca", ENERGY_GRADIENT_PROVIDERS)
        self.assertNotIn("psi4", ENERGY_GRADIENT_PROVIDERS)

        xtb_provider = get_energy_gradient_provider("xtb", {"charge": 0})
        aimnet_provider = get_energy_gradient_provider("aimnet_2", {"charge": 0})
        gaussian_provider = get_energy_gradient_provider(
            "gaussian",
            {"charge": 0, "multiplicity": 1, "method": "B97-D", "basis": "def2-SVP", "nprocs": 1},
        )
        orca_provider = get_energy_gradient_provider(
            "orca",
            {"charge": 0, "multiplicity": 1, "method": "B97-D", "basis": "def2-SVP", "nprocs": 1},
        )
        orca_alias_provider = get_energy_gradient_provider(
            "orca16",
            {"charge": 0, "multiplicity": 1, "method": "B97-D", "basis": "def2-SVP", "nprocs": 1},
        )
        self.assertTrue(hasattr(xtb_provider, "evaluate"))
        self.assertTrue(hasattr(aimnet_provider, "evaluate"))
        self.assertTrue(hasattr(gaussian_provider, "evaluate"))
        self.assertTrue(hasattr(orca_provider, "evaluate"))
        self.assertTrue(hasattr(orca_alias_provider, "evaluate"))

        with self.assertRaisesRegex(ValueError, "not 'psi4'"):
            get_energy_gradient_provider("psi4")

    def test_aimnet2_provider_loads_model_lazily(self):
        class FakeCalculator:
            def __init__(self, model, charge=0):
                self.model = model
                self.charge = charge

            def get_potential_energy(self, atoms=None):
                return 1.0

            def get_forces(self, atoms=None):
                return np.zeros((2, 3), dtype=float)

        provider = get_energy_gradient_provider("aimnet_2", {"charge": 1})
        molecule = Atoms(symbols=["H", "H"], positions=np.zeros((2, 3)))
        coordinates_bohr = np.zeros((2, 3), dtype=float)
        aimnet_backend = ModuleType("pyar.backends.aimnet_2")
        aimnet_backend.load_default_aimnet2_model = mock.Mock(return_value=object())
        aimnet_calculator = ModuleType("pyar.AIMNet2.calculators.aimnet2ase")
        aimnet_calculator.AIMNet2Calculator = FakeCalculator

        with (
            mock.patch.dict(
                sys.modules,
                {
                    "pyar.backends.aimnet_2": aimnet_backend,
                    "pyar.AIMNet2.calculators.aimnet2ase": aimnet_calculator,
                },
            )
        ):
            result = provider.evaluate(molecule, coordinates_bohr)

        aimnet_backend.load_default_aimnet2_model.assert_called_once_with()
        self.assertAlmostEqual(result.energy_hartree, 1.0 / Hartree)
        np.testing.assert_allclose(result.gradient_hartree_per_bohr, np.zeros((2, 3)))

    def test_orca_provider_runs_single_point_and_reads_engrad(self):
        provider = get_energy_gradient_provider(
            "orca",
            {
                "charge": 0,
                "multiplicity": 1,
                "method": "B97-D",
                "basis": "def2-SVP",
                "nprocs": 1,
                "scf_cycles": 50,
            },
        )
        molecule = Atoms(symbols=["H", "H"], positions=np.zeros((2, 3)))
        coordinates_bohr = np.zeros((2, 3), dtype=float)

        def fake_run(command, cwd=None, capture_output=None, text=None, check=None):
            inp_path = Path(command[1])
            engrad_path = inp_path.with_suffix(".engrad")
            engrad_path.write_text("2\n-1.234\n0.01\n0.02\n0.03\n0.04\n0.05\n0.06\n")
            return SimpleNamespace(returncode=0, stdout="FINAL SINGLE POINT ENERGY      -1.234000\n", stderr="")

        with mock.patch("pyar.energy_gradient_providers.require_executable", return_value="orca"), \
            mock.patch("pyar.energy_gradient_providers.subp.run", side_effect=fake_run):
            result = provider.evaluate(molecule, coordinates_bohr)

        self.assertAlmostEqual(result.energy_hartree, -1.234)
        np.testing.assert_allclose(result.gradient_hartree_per_bohr, np.asarray([[0.01, 0.02, 0.03], [0.04, 0.05, 0.06]]))

    def test_gaussian_provider_runs_force_and_reads_output(self):
        provider = get_energy_gradient_provider(
            "gaussian",
            {
                "charge": 0,
                "multiplicity": 1,
                "method": "B97-D",
                "basis": "def2-SVP",
                "nprocs": 1,
                "scf_cycles": 50,
            },
        )
        molecule = Atoms(symbols=["H", "H"], positions=np.zeros((2, 3)))
        coordinates_bohr = np.zeros((2, 3), dtype=float)

        gaussian_output = "\n".join(
            [
                " SCF Done:  E(RB97D) =  -1.234000     A.U. after   10 cycles",
                " Forces (Hartrees/Bohr)",
                " ---------------------------------------------------------------------",
                " 1 1  0.010000  0.020000  0.030000",
                " 2 1 -0.010000 -0.020000 -0.030000",
                " Normal termination of Gaussian 16",
            ]
        )

        def fake_run(command, cwd=None, capture_output=None, text=None, check=None):
            return SimpleNamespace(returncode=0, stdout=gaussian_output, stderr="")

        with mock.patch("pyar.energy_gradient_providers.require_executable", return_value="g16"), \
            mock.patch("pyar.energy_gradient_providers.subp.run", side_effect=fake_run):
            result = provider.evaluate(molecule, coordinates_bohr)

        self.assertAlmostEqual(result.energy_hartree, -1.234)
        np.testing.assert_allclose(
            result.gradient_hartree_per_bohr,
            np.asarray([[-0.01, -0.02, -0.03], [0.01, 0.02, 0.03]]),
        )

    def test_gaussian_provider_does_not_request_missing_pseudopotential_input(self):
        provider = get_energy_gradient_provider(
            "gaussian",
            {"method": "BP86", "basis": "def2-SVP", "nprocs": 1, "multiplicity": 1},
        )

        keyword = provider._build_keyword(Atoms("Fe", positions=[[0.0, 0.0, 0.0]]))

        self.assertNotIn("pseudo=read", keyword)

    def test_registration_adds_a_provider_without_touching_optimizer_code(self):
        class DummyProvider:
            def __init__(self, qc_params=None):
                self.qc_params = qc_params or {}

            def evaluate(self, molecule, coordinates_bohr):
                return EnergyGradientResult(1.0, np.full_like(coordinates_bohr, 2.0))

        original = dict(ENERGY_GRADIENT_PROVIDERS)
        original_capabilities = get_backend_capabilities("orca")
        try:
            register_energy_gradient_provider("orca", DummyProvider)
            provider = get_energy_gradient_provider("orca", {"charge": 0})
            self.assertIsInstance(provider, DummyProvider)
            self.assertTrue(backend_supports_geometry_optimization("orca"))
        finally:
            ENERGY_GRADIENT_PROVIDERS.clear()
            ENERGY_GRADIENT_PROVIDERS.update(original)
            register_backend_capabilities("orca", original_capabilities)

    def test_capability_without_provider_is_reported_cleanly(self):
        original_capabilities = get_backend_capabilities("missing_provider")
        try:
            register_backend_capabilities(
                "missing_provider",
                BackendCapabilities(energy_gradient=True),
            )
            with self.assertRaisesRegex(ValueError, "has no registered provider"):
                get_energy_gradient_provider("missing_provider")
        finally:
            register_backend_capabilities("missing_provider", original_capabilities)


if __name__ == "__main__":
    unittest.main()
