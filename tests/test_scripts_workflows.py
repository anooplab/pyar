import importlib
import os
import sys
import tempfile
import unittest
from types import SimpleNamespace
from unittest import mock

import numpy as np

from pyar.Molecule import Molecule
from pyar.reaction_state import ReactionStateError


class StandaloneWorkflowScriptTests(unittest.TestCase):
    def _import_in_tempdir(self, module_name):
        sys.modules.pop(module_name, None)
        return importlib.import_module(module_name)

    def test_explore_uses_current_trial_generation_signature(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cwd = os.getcwd()
            os.chdir(tmpdir)
            try:
                explore = self._import_in_tempdir("pyar.scripts.explore")
                seed = SimpleNamespace()
                monomer = SimpleNamespace(atoms_list=["H"])
                merged = object()
                with mock.patch.object(explore, "generate_points", return_value=[[1, 0, 0, 0, 0, 0]]) as generator:
                    with mock.patch.object(explore, "merge_two_molecules", return_value=merged):
                        result = explore.create_composite_molecule_wrapper(
                            seed,
                            monomer,
                            ["H"],
                            sequence_offset=3,
                        )
            finally:
                os.chdir(cwd)

        generator.assert_called_once_with(1, sequence_offset=3)
        self.assertIs(result, merged)

    def test_react_uses_current_reactor_signature(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cwd = os.getcwd()
            os.chdir(tmpdir)
            try:
                script = self._import_in_tempdir("pyar.scripts.react")
                arguments = SimpleNamespace(
                    input_files=["a.xyz", "b.xyz"],
                    how_many_orientations="4",
                    gmin=0.1,
                    gmax=0.5,
                    software="xtb",
                    index=0,
                )
                with mock.patch.object(script, "argument_parse", return_value=arguments):
                    with mock.patch.object(script.Molecule, "from_xyz", side_effect=[object(), object()]):
                        with mock.patch.object(script.reactor, "react") as react:
                            script.main()
            finally:
                os.chdir(cwd)

        react.assert_called_once()
        self.assertEqual(len(react.call_args.args), 8)
        self.assertEqual(react.call_args.args[-2:], (None, 2.3))

    def test_react_reports_restart_state_error_cleanly(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cwd = os.getcwd()
            os.chdir(tmpdir)
            try:
                script = self._import_in_tempdir("pyar.scripts.react")
                arguments = SimpleNamespace(
                    input_files=["a.xyz", "b.xyz"],
                    how_many_orientations="4",
                    gmin=0.1,
                    gmax=0.5,
                    software="xtb",
                    index=0,
                )
                with mock.patch.object(script, "argument_parse", return_value=arguments), \
                    mock.patch.object(script.Molecule, "from_xyz", side_effect=[object(), object()]), \
                    mock.patch.object(
                        script.reactor,
                        "react",
                        side_effect=ReactionStateError("Restart geometry snapshot is unavailable"),
                    ):
                    with self.assertRaises(SystemExit) as context:
                        script.main()
            finally:
                os.chdir(cwd)

        self.assertEqual(str(context.exception), "Restart geometry snapshot is unavailable")

    def test_reactor_handles_failed_optimization_without_copying_invalid_geometry(self):
        import pyar.reactor as reactor

        molecule = Molecule(
            ["H", "H"],
            np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.74]]),
            name="000_geom",
        )

        def fail_optimization(current, _params):
            current.coordinates = None
            current.energy = None
            return None

        with tempfile.TemporaryDirectory() as tmpdir:
            cwd = os.getcwd()
            os.chdir(tmpdir)
            try:
                with mock.patch.object(
                    reactor.pyar.interface.babel,
                    "make_inchi_string_from_xyz",
                    return_value="start-inchi",
                ):
                    with mock.patch.object(
                        reactor.pyar.interface.babel,
                        "make_smile_string_from_xyz",
                        return_value="start-smile",
                    ):
                        with mock.patch.object(reactor, "optimise", side_effect=fail_optimization):
                            result = reactor.optimize_all(
                                "gamma",
                                [molecule],
                                None,
                                tmpdir,
                                {"gamma": 1.0},
                            )
            finally:
                os.chdir(cwd)

        self.assertEqual(result, [])

    def test_reactor_with_gamma_copies_parameters(self):
        import pyar.reactor as reactor

        qc_params = {"software": "xtb", "gamma": 100.0, "basis": "def2-SVP"}
        updated = reactor.with_gamma(qc_params, 125.0)

        self.assertEqual(qc_params["gamma"], 100.0)
        self.assertEqual(updated["gamma"], 125.0)
        self.assertEqual(updated["software"], "xtb")
        self.assertIsNot(updated, qc_params)

    def test_reactor_relaxes_bonded_candidates_without_afir_bias(self):
        import pyar.reactor as reactor

        molecule = Molecule(
            ["H", "H"],
            np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.5]]),
            name="000_geom",
            fragments=[[0], [1]],
        )
        optimization_parameters = []

        def succeed(current, params):
            optimization_parameters.append(dict(params))
            current.energy = -1.0
            return True

        with tempfile.TemporaryDirectory() as tmpdir:
            cwd = os.getcwd()
            os.chdir(tmpdir)
            try:
                with mock.patch.object(reactor.file_manager, "make_directories", side_effect=lambda path: os.makedirs(path, exist_ok=True)), \
                    mock.patch.object(reactor.pyar.interface.babel, "make_inchi_string_from_xyz", return_value="same-inchi"), \
                    mock.patch.object(reactor.pyar.interface.babel, "make_smile_string_from_xyz", return_value="same-smile"), \
                    mock.patch.object(reactor, "optimise", side_effect=succeed):
                    reactor.optimize_all(
                        "gamma",
                        [molecule],
                        None,
                        tmpdir,
                        {"software": "xtb", "geometry_optimizer": "geometric", "gamma": 100.0},
                    )
            finally:
                os.chdir(cwd)

        self.assertEqual(optimization_parameters[0]["gamma"], 100.0)
        self.assertEqual(optimization_parameters[1]["gamma"], 0.0)

    def test_reactor_writes_disconnected_product_reference(self):
        import pyar.reactor as reactor

        molecule = Molecule(
            ["C", "H"],
            np.array([[0.0, 0.0, 0.0], [1.1, 0.0, 0.0]]),
            name="candidate",
            fragments=[[0], [1]],
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            output = os.path.join(tmpdir, "reference.xyz")
            reactor.write_disconnected_reference(molecule, output)
            reference = Molecule.from_xyz(output)

        self.assertAlmostEqual(reference.coordinates[1, 0] - reference.coordinates[0, 0], 101.1)
        self.assertAlmostEqual(molecule.coordinates[1, 0] - molecule.coordinates[0, 0], 1.1)

    def test_reactor_build_gamma_schedule_returns_numeric_values(self):
        import pyar.reactor as reactor

        schedule = reactor.build_gamma_schedule(0.1, 0.2, steps=3)

        self.assertEqual(len(schedule), 3)
        self.assertAlmostEqual(float(schedule[0]), 0.1)
        self.assertAlmostEqual(float(schedule[1]), 0.15)
        self.assertAlmostEqual(float(schedule[2]), 0.2)
        self.assertTrue(all(isinstance(value, np.floating) or isinstance(value, float) for value in schedule))

    def test_reactor_rejects_invalid_gamma_schedule(self):
        import pyar.reactor as reactor

        with self.assertRaisesRegex(ValueError, "non-negative"):
            reactor.build_gamma_schedule(-1.0, 100.0)
        with self.assertRaisesRegex(ValueError, "greater than or equal"):
            reactor.build_gamma_schedule(100.0, 10.0)

    def test_reactor_equal_gamma_limits_run_one_cycle(self):
        import pyar.reactor as reactor

        schedule = reactor.build_gamma_schedule(100.0, 100.0)

        self.assertEqual(schedule.tolist(), [100.0])

    def test_reactor_format_gamma_id_uses_zero_padding(self):
        import pyar.reactor as reactor

        self.assertEqual(reactor.format_gamma_id(0.1), "0000p1")
        self.assertEqual(reactor.format_gamma_id(12.9), "0012p9")
        self.assertEqual(reactor.format_gamma_id(100.0), "0100")

    def test_reactor_initialize_reaction_run_prepares_fresh_state(self):
        import pyar.reactor as reactor

        reactant_a = Molecule(["H"], np.array([[0.0, 0.0, 0.0]]), name="a")
        reactant_b = Molecule(["H"], np.array([[0.0, 0.0, 1.0]]), name="b")
        orientation = Molecule(
            ["H", "H"],
            np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]]),
            name="g1",
            fragments=[[0], [1]],
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            cwd = os.getcwd()
            os.chdir(tmpdir)
            try:
                with mock.patch.object(reactor.file_manager, "make_directories", side_effect=lambda path: os.makedirs(path, exist_ok=True)), \
                    mock.patch.object(reactor.trial_generation, "create_trial_geometries", return_value=[orientation]):
                    workdir, cwd_after, run_state, gamma_list, orientations, product_dir = reactor.initialize_reaction_run(
                        reactant_a,
                        reactant_b,
                        10.0,
                        20.0,
                        2,
                        {"software": "xtb"},
                        None,
                        2.3,
                    )
                self.assertTrue(run_state.state_file.exists())
                self.assertEqual(len(run_state.pending_molecules()), 1)
            finally:
                os.chdir(cwd)

        self.assertEqual(workdir, tmpdir)
        self.assertTrue(cwd_after.endswith("reaction"))
        self.assertEqual(len(gamma_list), 10)
        self.assertEqual(len(orientations), 1)
        self.assertTrue(product_dir.endswith("/products"))

    def test_reactor_does_not_delete_existing_directory_without_restart_state(self):
        import pyar.reactor as reactor

        reactant_a = Molecule(["H"], np.array([[0.0, 0.0, 0.0]]), name="a")
        reactant_b = Molecule(["H"], np.array([[0.0, 0.0, 1.0]]), name="b")
        with tempfile.TemporaryDirectory() as tmpdir:
            reaction_dir = os.path.join(tmpdir, "reaction")
            os.makedirs(reaction_dir)
            marker = os.path.join(reaction_dir, "previous_result.xyz")
            with open(marker, "w") as fp:
                fp.write("existing result\n")
            cwd = os.getcwd()
            os.chdir(tmpdir)
            try:
                with self.assertRaisesRegex(ReactionStateError, "has no resumable state"):
                    reactor.initialize_reaction_run(
                        reactant_a,
                        reactant_b,
                        10.0,
                        20.0,
                        2,
                        {"software": "xtb"},
                        None,
                        2.3,
                    )
            finally:
                os.chdir(cwd)

            self.assertTrue(os.path.isfile(marker))

    def test_reactor_uses_numeric_gamma_schedule(self):
        import pyar.reactor as reactor

        with tempfile.TemporaryDirectory() as tmpdir:
            reaction_dir = os.path.join(tmpdir, "reaction")
            os.makedirs(reaction_dir)
            cwd = os.getcwd()
            os.chdir(tmpdir)
            try:
                run_state = mock.Mock()
                with mock.patch.object(
                    reactor,
                    "initialize_reaction_run",
                    return_value=(tmpdir, reaction_dir, run_state, [0.1], [], os.path.join(reaction_dir, "products")),
                ), \
                    mock.patch.object(reactor, "optimize_all", return_value=[]) as optimize_all, \
                    mock.patch.object(reactor.os, "chdir"):
                    reactor.react(
                        object(),
                        object(),
                        0.1,
                        0.2,
                        1,
                        {"software": "xtb"},
                        None,
                        2.3,
                    )
            finally:
                os.chdir(cwd)

        self.assertGreaterEqual(optimize_all.call_count, 1)
        first_qc_params = optimize_all.call_args_list[0].args[-1]
        self.assertIsInstance(first_qc_params["gamma"], float)
        self.assertEqual(first_qc_params["gamma"], 0.1)
        run_state.complete_cycle.assert_called_once_with(0.1, [])
        run_state.finish.assert_called_once_with("completed_no_candidates")

    def test_trial_generation_make_composite_uses_population_offsets(self):
        script = self._import_in_tempdir("pyar.scripts.trial_generation")
        result = SimpleNamespace(title=None, mol_to_xyz=lambda _path: None)

        with tempfile.TemporaryDirectory() as tmpdir:
            cwd = os.getcwd()
            os.chdir(tmpdir)
            try:
                with mock.patch(
                    "sys.argv",
                    [
                        "pyar-trial-generation",
                        "-N",
                        "2",
                        "-i",
                        "seed.xyz",
                        "monomer.xyz",
                        "--make-composite",
                    ],
                ):
                    with mock.patch.object(script.Molecule, "from_xyz", return_value=object()):
                        with mock.patch.object(script.generation, "generate_points", return_value=[]) as generate:
                            with mock.patch.object(
                                script.generation,
                                "create_composite_molecule",
                                return_value=result,
                            ):
                                script.main()
            finally:
                os.chdir(cwd)

        self.assertEqual(generate.call_args_list[1].kwargs["sequence_offset"], 0)
        self.assertEqual(generate.call_args_list[2].kwargs["sequence_offset"], 1)


if __name__ == "__main__":
    unittest.main()
