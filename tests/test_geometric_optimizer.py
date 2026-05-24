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
from pyar.afir import restraints
from pyar.data.units import angstrom2bohr


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
                with mock.patch("pyar.interface.geometric._find_geometric_executable", return_value="geometric-optimize"):
                    geometry = optimiser.build_geometry(self.molecule, qc_params)
            finally:
                os.chdir(cwd)

        from pyar.interface.geometric import Geometric

        self.assertIsInstance(geometry, Geometric)
        self.assertEqual(geometry.gamma, 0.0)
        self.assertEqual(geometry.software, "xtb")

    def test_build_geometry_rejects_unsupported_backend_for_geometric(self):
        qc_params = {"software": "orca", "geometry_optimizer": "geometric", "gamma": 0.0}

        with self.assertRaisesRegex(ValueError, "supports only 'xtb' and 'aimnet_2'"):
            optimiser.build_geometry(self.molecule, qc_params)

    def test_geometric_calculator_adds_afir_only_when_gamma_nonzero(self):
        from pyar.interface.geometric import PyarGeometricCalculator

        atoms = Atoms(symbols=self.molecule.atoms_list, positions=self.molecule.coordinates)
        backend_energy = 1.5
        backend_forces = np.full((5, 3), 0.25)

        with mock.patch("pyar.interface.geometric._resolve_backend_evaluator", return_value=lambda a, p: (backend_energy, backend_forces)):
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
        from pyar.interface.geometric import PyarGeometricCalculator

        atoms = Atoms(symbols=self.molecule.atoms_list, positions=self.molecule.coordinates)
        backend_energy = 2.0
        backend_forces = np.full((5, 3), 0.25)
        afir_energy_hartree = 0.4
        afir_forces_hartree_per_bohr = np.full((5, 3), 0.1)

        with mock.patch("pyar.interface.geometric._resolve_backend_evaluator", return_value=lambda a, p: (backend_energy, backend_forces)), \
            mock.patch("pyar.interface.geometric.restraints.isotropic", return_value=(afir_energy_hartree, afir_forces_hartree_per_bohr)) as isotropic:
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

    def test_afir_term_points_fragments_toward_each_other(self):
        coordinates_bohr = angstrom2bohr(np.asarray([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]))
        _, afir_force = restraints.isotropic([[0], [1]], ["C", "H"], coordinates_bohr, 100.0)

        self.assertGreater(afir_force[0, 0], 0.0)
        self.assertLess(afir_force[1, 0], 0.0)

    def test_geometric_preserves_molecular_charge_and_spin_for_backend(self):
        from pyar.interface.geometric import Geometric

        molecule = SimpleNamespace(**self.molecule.__dict__)
        molecule.charge = -1
        molecule.multiplicity = 2
        molecule.scftype = "uhf"
        with tempfile.TemporaryDirectory() as tmpdir:
            cwd = os.getcwd()
            os.chdir(tmpdir)
            try:
                with mock.patch("pyar.interface.geometric._find_geometric_executable", return_value="geometric-optimize"):
                    geometry = Geometric(molecule, {"software": "xtb", "gamma": 0.0})
            finally:
                os.chdir(cwd)

        self.assertEqual(geometry.qc_params["charge"], -1)
        self.assertEqual(geometry.qc_params["multiplicity"], 2)
        self.assertEqual(geometry.qc_params["scftype"], "uhf")

    def test_geometric_uses_valid_convergence_preset_syntax(self):
        from pyar.interface.geometric import Geometric

        with tempfile.TemporaryDirectory() as tmpdir:
            cwd = os.getcwd()
            os.chdir(tmpdir)
            try:
                with mock.patch("pyar.interface.geometric._find_geometric_executable", return_value="geometric-optimize"):
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

    def test_geometric_executable_lookup_keeps_virtual_environment_path(self):
        from pyar.interface.geometric import _find_geometric_executable

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
                mock.patch("pyar.interface.geometric.require_executable") as fallback:
                self.assertEqual(_find_geometric_executable(), str(executable))

            fallback.assert_not_called()

    def test_transition_target_is_reserved_for_future_implementation(self):
        from pyar.interface.geometric import Geometric

        with tempfile.TemporaryDirectory() as tmpdir:
            cwd = os.getcwd()
            os.chdir(tmpdir)
            try:
                with mock.patch("pyar.interface.geometric._find_geometric_executable", return_value="geometric-optimize"):
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
