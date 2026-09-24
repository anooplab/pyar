import json
import os
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import numpy as np
from ase import Atoms
from ase.units import Bohr, Hartree

from pyar import optimiser
from pyar.biases import afir as restraints
from pyar.data.units import angstrom2bohr
from pyar.energy_gradient_providers import EnergyGradientResult


class GeometricOptimizerTests(unittest.TestCase):
    def setUp(self):
        self.molecule = SimpleNamespace(
            name="geom",
            title="geom",
            atoms_list=["C", "H", "H", "H", "H"],
            number_of_atoms=5,
            charge=0,
            multiplicity=1,
            scftype="rhf",
            coordinates=np.asarray(
                [
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [0.0, 0.0, 1.0],
                    [-1.0, 0.0, 0.0],
                ],
                dtype=float,
            ),
            fragments=[[0], [1, 2, 3, 4]],
        )

    def test_build_geometry_uses_geometric_wrapper_for_supported_backend(self):
        qc_params = {"software": "xtb", "geometry_optimizer": "geometric", "gamma": 0.0}

        with tempfile.TemporaryDirectory() as tmpdir:
            cwd = os.getcwd()
            os.chdir(tmpdir)
            try:
                with mock.patch("pyar.backends.geometric._find_geometric_executable", return_value="geometric-optimize"):
                    geometry = optimiser.build_geometry(self.molecule, qc_params)
            finally:
                os.chdir(cwd)

        from pyar.backends.geometric import Geometric

        self.assertIsInstance(geometry, Geometric)
        self.assertEqual(geometry.gamma, 0.0)
        self.assertEqual(geometry.software, "xtb")

    def test_build_geometry_rejects_unsupported_backend_for_geometric(self):
        qc_params = {"software": "mopac", "geometry_optimizer": "geometric", "gamma": 0.0}

        with self.assertRaisesRegex(ValueError, "does not expose Cartesian energy and gradients"):
            optimiser.build_geometry(self.molecule, qc_params)

    def test_geometric_calculator_adds_afir_only_when_gamma_nonzero(self):
        from pyar.backends.geometric import PyarGeometricCalculator

        atoms = Atoms(symbols=self.molecule.atoms_list, positions=self.molecule.coordinates)
        backend_energy = 1.5
        backend_forces = np.full((5, 3), 0.25)

        backend_result = EnergyGradientResult(
            backend_energy / Hartree,
            -backend_forces * Bohr / Hartree,
        )

        class DummyProvider:
            def evaluate(self, molecule, coordinates_bohr):
                return backend_result

        with mock.patch("pyar.backends.geometric._resolve_backend_evaluator", return_value=DummyProvider()):
            calculator = PyarGeometricCalculator(
                {"software": "xtb", "gamma": 0.0, "charge": 0},
                fragment_indices=self.molecule.fragments,
            )

        with tempfile.TemporaryDirectory() as tmpdir:
            cwd = os.getcwd()
            os.chdir(tmpdir)
            try:
                calculator.calculate(atoms=atoms, properties=["energy", "forces"])
            finally:
                os.chdir(cwd)

        self.assertAlmostEqual(calculator.results["energy"], backend_energy)
        np.testing.assert_allclose(calculator.results["forces"], backend_forces)

    def test_geometric_calculator_combines_backend_and_afir(self):
        from pyar.backends.geometric import PyarGeometricCalculator

        atoms = Atoms(symbols=self.molecule.atoms_list, positions=self.molecule.coordinates)
        backend_energy = 2.0
        backend_forces = np.full((5, 3), 0.25)
        afir_energy_hartree = 0.4
        afir_forces_hartree_per_bohr = np.full((5, 3), 0.1)

        backend_result = EnergyGradientResult(
            backend_energy / Hartree,
            -backend_forces * Bohr / Hartree,
        )

        class DummyProvider:
            def evaluate(self, molecule, coordinates_bohr):
                return backend_result

        with mock.patch("pyar.backends.geometric._resolve_backend_evaluator", return_value=DummyProvider()), \
            mock.patch("pyar.backends.geometric.restraints.isotropic", return_value=(afir_energy_hartree, afir_forces_hartree_per_bohr)) as isotropic:
            calculator = PyarGeometricCalculator(
                {"software": "xtb", "gamma": 37.5, "charge": 0},
                fragment_indices=self.molecule.fragments,
            )
            with tempfile.TemporaryDirectory() as tmpdir:
                cwd = os.getcwd()
                os.chdir(tmpdir)
                try:
                    calculator.calculate(atoms=atoms, properties=["energy", "forces"])
                finally:
                    os.chdir(cwd)

        isotropic.assert_called_once()
        self.assertAlmostEqual(calculator.results["energy"], backend_energy + afir_energy_hartree * Hartree)
        np.testing.assert_allclose(
            calculator.results["forces"],
            backend_forces + (afir_forces_hartree_per_bohr * Hartree / Bohr),
        )

    def test_geometric_calculator_selects_softmin_bias(self):
        from pyar.backends.geometric import PyarGeometricCalculator

        atoms = Atoms(symbols=self.molecule.atoms_list, positions=self.molecule.coordinates)
        backend_result = EnergyGradientResult(0.0, np.zeros((5, 3)))
        softmin_forces = np.full((5, 3), 0.1)

        class DummyProvider:
            def evaluate(self, molecule, coordinates_bohr):
                return backend_result

        with mock.patch("pyar.backends.geometric._resolve_backend_evaluator", return_value=DummyProvider()), \
            mock.patch("pyar.backends.geometric.softmin.softmin", return_value=(0.4, softmin_forces)) as selected_bias:
            calculator = PyarGeometricCalculator(
                {
                    "software": "xtb",
                    "gamma": 37.5,
                    "bias_potential": "softmin",
                    "softmin_beta": 2.5,
                    "charge": 0,
                },
                fragment_indices=self.molecule.fragments,
            )
            with tempfile.TemporaryDirectory() as tmpdir:
                cwd = os.getcwd()
                os.chdir(tmpdir)
                try:
                    calculator.calculate(atoms=atoms, properties=["energy", "forces"])
                finally:
                    os.chdir(cwd)

        selected_bias.assert_called_once()
        self.assertEqual(selected_bias.call_args.kwargs["beta"], 2.5)
        self.assertEqual(calculator.bias_potential, "softmin")
        self.assertEqual(calculator.softmin_beta, 2.5)
        self.assertAlmostEqual(calculator.results["energy"], 0.4 * Hartree)
        np.testing.assert_allclose(calculator.results["forces"], softmin_forces * Hartree / Bohr)

    def test_fixed_controller_preserves_direct_afir_and_softmin_objective(self):
        from pyar.backends.geometric import PyarGeometricCalculator
        from pyar.biases.softmin import softmin

        atoms = Atoms(symbols=self.molecule.atoms_list, positions=self.molecule.coordinates)
        backend_gradient = np.arange(15, dtype=float).reshape(5, 3) / 1000.0
        backend_energy_hartree = -0.25
        backend_result = EnergyGradientResult(backend_energy_hartree, backend_gradient)

        class DummyProvider:
            def evaluate(self, molecule, coordinates_bohr):
                return backend_result

        coordinates_bohr = angstrom2bohr(self.molecule.coordinates)
        gamma = 37.5
        for potential in ("afir", "softmin"):
            with self.subTest(potential=potential), tempfile.TemporaryDirectory() as tmpdir:
                cwd = os.getcwd()
                os.chdir(tmpdir)
                try:
                    with mock.patch("pyar.backends.geometric._resolve_backend_evaluator",
                                    return_value=DummyProvider()):
                        calculator = PyarGeometricCalculator(
                            {"software": "xtb", "gamma": gamma, "bias_potential": potential,
                             "softmin_beta": 1.7}, self.molecule.fragments
                        )
                        calculator.calculate(atoms=atoms, properties=["energy", "forces"])
                finally:
                    os.chdir(cwd)

            if potential == "afir":
                bias_energy_hartree, bias_forces = restraints.isotropic(
                    self.molecule.fragments, self.molecule.atoms_list, coordinates_bohr, gamma
                )
            else:
                bias_energy_hartree, bias_forces = softmin(
                    self.molecule.fragments, self.molecule.atoms_list, coordinates_bohr,
                    gamma, beta=1.7
                )
            expected_energy = (backend_energy_hartree + bias_energy_hartree) * Hartree
            expected_forces = (-backend_gradient + bias_forces) * Hartree / Bohr
            self.assertAlmostEqual(calculator.results["energy"], expected_energy, places=12)
            np.testing.assert_allclose(calculator.results["forces"], expected_forces, rtol=1e-13, atol=1e-13)
            self.assertEqual(calculator.bias_controller.decision.policy, "fixed")

    def test_geometric_calculator_rejects_unknown_bias_potential(self):
        from pyar.backends.geometric import PyarGeometricCalculator

        with self.assertRaisesRegex(ValueError, "Unsupported bias potential"):
            PyarGeometricCalculator({"software": "xtb", "bias_potential": "unknown"})

    def test_afir_term_points_fragments_toward_each_other(self):
        coordinates_bohr = angstrom2bohr(np.asarray([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]))
        _, afir_force = restraints.isotropic([[0], [1]], ["C", "H"], coordinates_bohr, 100.0)

        self.assertGreater(afir_force[0, 0], 0.0)
        self.assertLess(afir_force[1, 0], 0.0)

    def test_geometric_preserves_molecular_charge_and_spin_for_backend(self):
        from pyar.backends.geometric import Geometric

        molecule = SimpleNamespace(**self.molecule.__dict__)
        molecule.charge = -1
        molecule.multiplicity = 2
        molecule.scftype = "uhf"
        with tempfile.TemporaryDirectory() as tmpdir:
            cwd = os.getcwd()
            os.chdir(tmpdir)
            try:
                with mock.patch("pyar.backends.geometric._find_geometric_executable", return_value="geometric-optimize"):
                    geometry = Geometric(molecule, {"software": "xtb", "gamma": 0.0})
            finally:
                os.chdir(cwd)

        self.assertEqual(geometry.qc_params["charge"], -1)
        self.assertEqual(geometry.qc_params["multiplicity"], 2)
        self.assertEqual(geometry.qc_params["scftype"], "uhf")

    def test_geometric_uses_valid_convergence_preset_syntax(self):
        from pyar.backends.geometric import Geometric

        with tempfile.TemporaryDirectory() as tmpdir:
            cwd = os.getcwd()
            os.chdir(tmpdir)
            try:
                with mock.patch("pyar.backends.geometric._find_geometric_executable", return_value="geometric-optimize"):
                    geometry = Geometric(
                        self.molecule,
                        {"software": "xtb", "gamma": 0.0, "opt_threshold": "tight"},
                    )
                    command = geometry._build_command()
            finally:
                os.chdir(cwd)

        convergence_index = command.index("--converge")
        self.assertEqual(command[convergence_index:convergence_index + 3], ["--converge", "set", "GAU_TIGHT"])
        self.assertLess(command.index(geometry.start_xyz_file), convergence_index)

    def test_geometric_returns_cycle_exceeded_when_endpoint_is_preserved(self):
        from pyar.backends import write_xyz
        from pyar.backends.geometric import Geometric

        with tempfile.TemporaryDirectory() as tmpdir:
            cwd = os.getcwd()
            os.chdir(tmpdir)
            try:
                with mock.patch("pyar.backends.geometric._find_geometric_executable",
                                return_value="geometric-optimize"):
                    geometry = Geometric(
                        self.molecule,
                        {"software": "xtb", "gamma": 0.0,
                         "geometry_optimizer": "geometric"},
                    )

                def fake_run(command, **kwargs):
                    kwargs["stdout"].write("Maximum iterations reached (1000); increase --maxiter for more\n")
                    write_xyz(self.molecule.atoms_list, self.molecule.coordinates,
                              "trial_geom_optim.xyz")
                    Path("pyar_geometric_state.json").write_text(json.dumps({
                        "backend_energy_hartree": -1.25,
                        "optimization_status": "cycle_exceeded",
                        "positions_angstrom": self.molecule.coordinates.tolist(),
                    }))
                    return mock.Mock(returncode=1)

                with mock.patch("pyar.backends.geometric.subp.run", side_effect=fake_run):
                    self.assertEqual(geometry.optimize(), "CycleExceeded")
                self.assertEqual(geometry.energy, -1.25)
            finally:
                os.chdir(cwd)

    def test_geometric_executable_lookup_keeps_virtual_environment_path(self):
        from pyar.backends.geometric import _find_geometric_executable

        with tempfile.TemporaryDirectory() as tmpdir:
            bin_dir = Path(tmpdir) / "bin"
            bin_dir.mkdir()
            system_python = Path(tmpdir) / "system-python"
            system_python.touch()
            environment_python = bin_dir / "python"
            environment_python.symlink_to(system_python)
            executable = bin_dir / "geometric-optimize"
            executable.touch()

            with mock.patch.object(sys, "executable", str(environment_python)), \
                mock.patch("pyar.backends.geometric.require_executable") as fallback:
                self.assertEqual(_find_geometric_executable(), str(executable))

            fallback.assert_not_called()

    def test_transition_target_is_reserved_for_future_implementation(self):
        from pyar.backends.geometric import Geometric

        with tempfile.TemporaryDirectory() as tmpdir:
            cwd = os.getcwd()
            os.chdir(tmpdir)
            try:
                with mock.patch("pyar.backends.geometric._find_geometric_executable", return_value="geometric-optimize"):
                    geometry = Geometric(
                        self.molecule,
                        {"software": "xtb", "gamma": 0.0, "opt_target": "ts"},
                    )
                    with self.assertRaisesRegex(NotImplementedError, "reserved for a future"):
                        geometry._build_command()
            finally:
                os.chdir(cwd)


if __name__ == "__main__":
    unittest.main()
