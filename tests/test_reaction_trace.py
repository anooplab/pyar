import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np
from ase import Atoms
from ase.units import Bohr, Hartree

from pyar.energy_gradient_providers import EnergyGradientResult
from pyar.reaction_analysis import analyse_reaction_trace
from pyar.reaction_trace import bond_changes, infer_bonds, load_trace_records


class ReactionTraceTests(unittest.TestCase):
    def test_bond_change_heuristic_detects_forming_and_breaking_bonds(self):
        symbols = ["C", "H"]
        far_geometry = np.asarray([[0.0, 0.0, 0.0], [2.2, 0.0, 0.0]])
        near_geometry = np.asarray([[0.0, 0.0, 0.0], [1.05, 0.0, 0.0]])
        rebound_geometry = np.asarray([[0.0, 0.0, 0.0], [2.4, 0.0, 0.0]])

        far_bonds = infer_bonds(symbols, far_geometry)
        near_bonds = infer_bonds(symbols, near_geometry, far_bonds)
        rebound_bonds = infer_bonds(symbols, rebound_geometry, near_bonds)

        formed, broken = bond_changes(far_bonds, near_bonds)
        rebound_formed, rebound_broken = bond_changes(near_bonds, rebound_bonds)

        self.assertEqual(far_bonds, set())
        self.assertEqual(near_bonds, {(0, 1)})
        self.assertEqual(formed, [(0, 1)])
        self.assertEqual(broken, [])
        self.assertEqual(rebound_bonds, set())
        self.assertEqual(rebound_formed, [])
        self.assertEqual(rebound_broken, [(0, 1)])

    def test_geometric_calculator_records_trace_files(self):
        from pyar.backends.geometric import PyarGeometricCalculator

        atoms = Atoms(symbols=["C", "H"], positions=np.asarray([[0.0, 0.0, 0.0], [1.1, 0.0, 0.0]]))
        backend_forces = np.asarray([[0.10, 0.00, 0.00], [-0.10, 0.00, 0.00]])
        afir_forces = np.asarray([[0.03, 0.00, 0.00], [-0.03, 0.00, 0.00]])
        backend_result = EnergyGradientResult(
            1.2,
            -backend_forces * Bohr / Hartree,
        )

        class DummyProvider:
            def evaluate(self, molecule, coordinates_bohr):
                return backend_result

        with tempfile.TemporaryDirectory() as tmpdir:
            cwd = os.getcwd()
            os.chdir(tmpdir)
            try:
                with mock.patch("pyar.backends.geometric._resolve_backend_evaluator", return_value=DummyProvider()), \
                    mock.patch("pyar.backends.geometric.restraints.isotropic", return_value=(0.4, afir_forces * Bohr / Hartree)):
                    calculator = PyarGeometricCalculator(
                        {
                            "software": "xtb",
                            "gamma": 50.0,
                            "trace_enabled": True,
                            "charge": 0,
                        },
                        fragment_indices=[[0], [1]],
                    )
                    calculator.calculate(atoms=atoms, properties=["energy", "forces"])
                    atoms.positions[1, 0] = 1.2
                    calculator.calculate(atoms=atoms, properties=["energy", "forces"])
            finally:
                os.chdir(cwd)

            trace_records = load_trace_records(Path(tmpdir) / "reaction_trace")
            self.assertEqual(len(trace_records), 2)
            self.assertAlmostEqual(trace_records[0]["backend_energy_hartree"], 1.2)
            self.assertAlmostEqual(trace_records[0]["total_energy_hartree"], 1.6)
            self.assertIn("current_bonds", trace_records[0])
            self.assertTrue((Path(tmpdir) / "reaction_trace" / "steps" / "step_000000.xyz").exists())
            self.assertTrue((Path(tmpdir) / "reaction_trace" / "steps" / "step_000001.xyz").exists())

    def test_trace_pipeline_smoke_flow_uses_mocked_calculator_calls(self):
        from pyar.backends.geometric import PyarGeometricCalculator

        atoms = Atoms(symbols=["C", "H"], positions=np.asarray([[0.0, 0.0, 0.0], [1.8, 0.0, 0.0]]))

        backend_results = [
            EnergyGradientResult(0.6, np.zeros((2, 3))),
            EnergyGradientResult(1.7, np.zeros((2, 3))),
            EnergyGradientResult(0.9, np.zeros((2, 3))),
        ]
        afir_records = [
            (0.1, np.zeros((2, 3))),
            (0.2, np.zeros((2, 3))),
            (0.3, np.zeros((2, 3))),
        ]

        class DummyProvider:
            def __init__(self):
                self.calls = 0

            def evaluate(self, molecule, coordinates_bohr):
                result = backend_results[self.calls]
                self.calls += 1
                return result

        with tempfile.TemporaryDirectory() as tmpdir:
            cwd = os.getcwd()
            os.chdir(tmpdir)
            try:
                with mock.patch("pyar.backends.geometric._resolve_backend_evaluator", return_value=DummyProvider()), \
                    mock.patch("pyar.backends.geometric.restraints.isotropic", side_effect=afir_records):
                    calculator = PyarGeometricCalculator(
                        {
                            "software": "xtb",
                            "gamma": 75.0,
                            "trace_enabled": True,
                            "charge": 0,
                        },
                        fragment_indices=[[0], [1]],
                    )

                    calculator.calculate(atoms=atoms, properties=["energy", "forces"])
                    atoms.positions[1, 0] = 1.35
                    calculator.calculate(atoms=atoms, properties=["energy", "forces"])
                    atoms.positions[1, 0] = 1.05
                    calculator.calculate(atoms=atoms, properties=["energy", "forces"])
            finally:
                os.chdir(cwd)

            summary = analyse_reaction_trace(Path(tmpdir))
            trace_records = load_trace_records(Path(tmpdir) / "reaction_trace")
            self.assertEqual(len(trace_records), 3)
            self.assertIsNotNone(summary)
            self.assertTrue((Path(tmpdir) / "path_summary.csv").exists())
            self.assertTrue((Path(tmpdir) / "candidate_ts" / "highest_backend_energy.xyz").exists())
            self.assertTrue((Path(tmpdir) / "candidate_ts" / "pre_product_geometry.xyz").exists())
            self.assertTrue((Path(tmpdir) / "candidate_ts" / "max_bond_change.xyz").exists())
            self.assertTrue((Path(tmpdir) / "candidate_ts" / "highest_total_energy.xyz").exists())
            self.assertEqual(summary["highest_backend_energy_index"], 1)
            self.assertIn("2", Path(tmpdir, "candidate_ts", "highest_backend_energy.xyz").read_text())

    def test_reaction_trace_analysis_selects_backend_energy_not_total_energy(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            job_dir = Path(tmpdir)
            trace_dir = job_dir / "reaction_trace"
            trace_dir.mkdir()

            records = [
                {
                    "step_index": 0,
                    "symbols": ["C", "H"],
                    "coordinates_angstrom": [[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]],
                    "backend_energy_hartree": 0.5,
                    "afir_energy_hartree": 0.1,
                    "total_energy_hartree": 0.6,
                    "backend_force_norm": 1.0,
                    "afir_force_norm": 0.2,
                    "total_force_norm": 1.2,
                    "max_force": 0.8,
                    "current_bonds": [],
                    "formed_bonds": [],
                    "broken_bonds": [],
                    "bond_change_count": 0,
                    "min_interfragment_distance_angstrom": 1.5,
                },
                {
                    "step_index": 1,
                    "symbols": ["C", "H"],
                    "coordinates_angstrom": [[0.0, 0.0, 0.0], [1.1, 0.0, 0.0]],
                    "backend_energy_hartree": 2.5,
                    "afir_energy_hartree": 0.5,
                    "total_energy_hartree": 3.0,
                    "backend_force_norm": 2.0,
                    "afir_force_norm": 0.3,
                    "total_force_norm": 2.3,
                    "max_force": 1.7,
                    "current_bonds": [[0, 1]],
                    "formed_bonds": [[0, 1]],
                    "broken_bonds": [],
                    "bond_change_count": 1,
                    "min_interfragment_distance_angstrom": 1.1,
                },
            ]

            with (trace_dir / "trace.jsonl").open("w", encoding="utf-8") as fp:
                for record in records:
                    json.dump(record, fp)
                    fp.write("\n")

            summary = analyse_reaction_trace(job_dir)

            self.assertIsNotNone(summary)
            self.assertTrue((job_dir / "path_summary.csv").exists())
            self.assertTrue((job_dir / "candidate_ts" / "highest_backend_energy.xyz").exists())
            self.assertTrue((job_dir / "candidate_ts" / "highest_total_energy.xyz").exists())
            self.assertTrue((job_dir / "candidate_ts" / "pre_product_geometry.xyz").exists())
            self.assertTrue((job_dir / "candidate_ts" / "max_bond_change.xyz").exists())

            highest_backend = Path(job_dir / "candidate_ts" / "highest_backend_energy.xyz").read_text()
            highest_total = Path(job_dir / "candidate_ts" / "highest_total_energy.xyz").read_text()
            self.assertIn("2.5", highest_backend)
            self.assertIn("3.0", highest_total)
