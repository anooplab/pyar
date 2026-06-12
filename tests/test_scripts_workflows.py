import importlib
import os
import sys
import tempfile
import unittest
from types import SimpleNamespace
from pathlib import Path
from unittest import mock

import numpy as np

from pyar.core.molecule import Molecule
from pyar.reaction_state import ReactionStateError
from pyar.workflow_results import ReactionResult


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
                        with mock.patch.object(script.reaction_workflow, "react") as react:
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
                        script.reaction_workflow,
                        "react",
                        side_effect=ReactionStateError("Restart geometry snapshot is unavailable"),
                    ):
                    with self.assertRaises(SystemExit) as context:
                        script.main()
            finally:
                os.chdir(cwd)

        self.assertEqual(str(context.exception), "Restart geometry snapshot is unavailable")

    def test_react_selects_geometric_from_backend_capability(self):
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
                    software="future_provider",
                    geometry_optimizer=None,
                    opt_target="minimum",
                    index=0,
                )
                with mock.patch.object(script, "argument_parse", return_value=arguments), \
                    mock.patch.object(script.Molecule, "from_xyz", side_effect=[object(), object()]), \
                    mock.patch.object(script, "backend_supports_geometry_optimization", return_value=True) as supports, \
                    mock.patch.object(script.reaction_workflow, "react") as react:
                    script.main()
            finally:
                os.chdir(cwd)

        supports.assert_called_once_with("future_provider")
        self.assertEqual(react.call_args.args[5]["geometry_optimizer"], "geometric")
        self.assertEqual(react.call_args.args[5]["method"], "BP86")
        self.assertEqual(react.call_args.args[5]["basis"], "def2-SVP")
        self.assertEqual(react.call_args.args[5]["scf_cycles"], 1000)
        self.assertEqual(react.call_args.args[5]["nprocs"], 8)

    def test_conformer_script_calls_workflow(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cwd = os.getcwd()
            os.chdir(tmpdir)
            try:
                script = self._import_in_tempdir("pyar.scripts.conformer")
                result = SimpleNamespace(
                    status="completed",
                    run_directory=str(Path(tmpdir) / "conformers"),
                    selected_paths=("conformers/selected/conf_0000.xyz",),
                )
                with mock.patch.object(script, "conformer_search", return_value=result) as search:
                    script.main(["CCO", "--num-conformers", "5", "--top-n", "2"])
            finally:
                os.chdir(cwd)

        search.assert_called_once()
        self.assertEqual(search.call_args.args[0], "CCO")
        self.assertEqual(search.call_args.kwargs["num_conformers"], 5)
        self.assertEqual(search.call_args.kwargs["top_n"], 2)
        self.assertIsNone(search.call_args.kwargs["qc_params"])

    def test_conformer_script_preflights_backend_refinement(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cwd = os.getcwd()
            os.chdir(tmpdir)
            try:
                script = self._import_in_tempdir("pyar.scripts.conformer")
                result = SimpleNamespace(
                    status="completed",
                    run_directory=str(Path(tmpdir) / "conformers"),
                    selected_paths=(),
                )
                with mock.patch("pyar.cli._preflight_cli_requirements") as preflight:
                    with mock.patch.object(script, "conformer_search", return_value=result) as search:
                        script.main(["CCO", "--software", "xtb"])
            finally:
                os.chdir(cwd)

        preflight.assert_called_once_with("conformer", "xtb", "native")
        self.assertEqual(search.call_args.kwargs["qc_params"]["software"], "xtb")

    def test_conformer_script_reports_workflow_error_cleanly(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cwd = os.getcwd()
            os.chdir(tmpdir)
            try:
                script = self._import_in_tempdir("pyar.scripts.conformer")
                with mock.patch.object(
                    script,
                    "conformer_search",
                    side_effect=script.ConformerWorkflowError("RDKit did not generate any conformers"),
                ):
                    with self.assertRaises(SystemExit) as context:
                        script.main(["CCO"])
            finally:
                os.chdir(cwd)

        self.assertEqual(str(context.exception), "RDKit did not generate any conformers")

    def test_reaction_trace_script_analyzes_and_plots(self):
        import json

        with tempfile.TemporaryDirectory() as tmpdir:
            cwd = os.getcwd()
            os.chdir(tmpdir)
            try:
                script = self._import_in_tempdir("pyar.scripts.reaction_trace")
                trace_dir = Path(tmpdir) / "reaction_trace"
                trace_dir.mkdir()
                record = {
                    "step_index": 0,
                    "symbols": ["H", "H"],
                    "coordinates_angstrom": [[0.0, 0.0, 0.0], [0.7, 0.0, 0.0]],
                    "backend_energy_hartree": -1.0,
                    "afir_energy_hartree": 0.0,
                    "total_energy_hartree": -1.0,
                    "backend_force_norm": 0.1,
                    "afir_force_norm": 0.0,
                    "total_force_norm": 0.1,
                    "max_force": 0.1,
                    "current_bonds": [[0, 1]],
                    "formed_bonds": [[0, 1]],
                    "broken_bonds": [],
                    "bond_change_count": 1,
                    "min_interfragment_distance_angstrom": 0.7,
                }
                with (trace_dir / "trace.jsonl").open("w", encoding="utf-8") as fp:
                    json.dump(record, fp)
                    fp.write("\n")

                with mock.patch.object(script, "argument_parse", return_value=SimpleNamespace(
                    path=tmpdir,
                    plot=True,
                    plot_directory=None,
                )):
                    script.main()

                self.assertTrue((Path(tmpdir) / "path_summary.csv").exists())
                self.assertTrue((Path(tmpdir) / "trace_plots" / "reaction_trace_energies.png").exists())
                self.assertTrue((Path(tmpdir) / "trace_plots" / "reaction_trace_metrics.png").exists())
            finally:
                os.chdir(cwd)

    def test_reaction_trace_script_plot_only_skips_summary_generation(self):
        import json

        with tempfile.TemporaryDirectory() as tmpdir:
            cwd = os.getcwd()
            os.chdir(tmpdir)
            try:
                script = self._import_in_tempdir("pyar.scripts.reaction_trace")
                trace_dir = Path(tmpdir) / "reaction_trace"
                trace_dir.mkdir()
                record = {
                    "step_index": 0,
                    "symbols": ["H", "H"],
                    "coordinates_angstrom": [[0.0, 0.0, 0.0], [0.7, 0.0, 0.0]],
                    "backend_energy_hartree": -1.0,
                    "afir_energy_hartree": 0.0,
                    "total_energy_hartree": -1.0,
                    "backend_force_norm": 0.1,
                    "afir_force_norm": 0.0,
                    "total_force_norm": 0.1,
                    "max_force": 0.1,
                    "current_bonds": [[0, 1]],
                    "formed_bonds": [[0, 1]],
                    "broken_bonds": [],
                    "bond_change_count": 1,
                    "min_interfragment_distance_angstrom": 0.7,
                }
                with (trace_dir / "trace.jsonl").open("w", encoding="utf-8") as fp:
                    json.dump(record, fp)
                    fp.write("\n")

                with mock.patch.object(script, "argument_parse", return_value=SimpleNamespace(
                    path=tmpdir,
                    plot=False,
                    plot_only=True,
                    plot_directory=None,
                )):
                    script.main()

                self.assertFalse((Path(tmpdir) / "path_summary.csv").exists())
                self.assertTrue((Path(tmpdir) / "trace_plots" / "reaction_trace_energies.png").exists())
                self.assertTrue((Path(tmpdir) / "trace_plots" / "reaction_trace_metrics.png").exists())
            finally:
                os.chdir(cwd)

    def test_react_rejects_geometric_backend_without_provider_before_run(self):
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
                    software="mopac",
                    geometry_optimizer="geometric",
                    opt_target="minimum",
                    index=0,
                )
                with mock.patch.object(script, "argument_parse", return_value=arguments), \
                    mock.patch.object(script.Molecule, "from_xyz", side_effect=[object(), object()]), \
                    mock.patch.object(script.reaction_workflow, "react") as react:
                    with self.assertRaisesRegex(
                        SystemExit,
                        "does not expose Cartesian energy and gradients",
                    ):
                        script.main()
            finally:
                os.chdir(cwd)

        react.assert_not_called()

    def test_reactor_handles_failed_optimization_without_copying_invalid_geometry(self):
        from pyar.workflows import reaction as reactor

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
                    reactor,
                    "molecule_identity_from_xyz",
                    return_value={"inchi": "start-inchi", "smiles": "start-smile"},
                ), \
                    mock.patch.object(reactor, "optimise", side_effect=fail_optimization):
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
        from pyar.workflows import reaction as reactor

        qc_params = {"software": "xtb", "gamma": 100.0, "basis": "def2-SVP"}
        updated = reactor.with_gamma(qc_params, 125.0)

        self.assertEqual(qc_params["gamma"], 100.0)
        self.assertEqual(updated["gamma"], 125.0)
        self.assertEqual(updated["software"], "xtb")
        self.assertIsNot(updated, qc_params)

    def test_reactor_relaxes_bonded_candidates_without_afir_bias(self):
        from pyar.workflows import reaction as reactor

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
                    mock.patch.object(
                        reactor,
                        "molecule_identity_from_xyz",
                        return_value={"inchi": "same-inchi", "smiles": "same-smile"},
                    ), \
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

    def test_reactor_deduplicates_by_canonical_inchi_not_smiles_format(self):
        from pyar.workflows import reaction as reactor

        original_registry = dict(reactor.saved_product_identities)
        try:
            reactor.saved_product_identities.clear()
            reactor.saved_product_identities["existing"] = {
                "inchi": "same-inchi",
                "smiles": "smiles-a",
            }
            self.assertTrue(
                reactor._is_known_product({"inchi": "same-inchi", "smiles": "smiles-b"})
            )
            self.assertFalse(
                reactor._is_known_product({"inchi": "other-inchi", "smiles": "smiles-a"})
            )
        finally:
            reactor.saved_product_identities.clear()
            reactor.saved_product_identities.update(original_registry)

    def test_unbiased_relaxation_restores_job_name_after_failure(self):
        from pyar.workflows import reaction as reactor

        molecule = Molecule(
            ["C", "H"],
            np.array([[0.0, 0.0, 0.0], [1.1, 0.0, 0.0]]),
            name="candidate",
            fragments=[[0], [1]],
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            cwd = os.getcwd()
            os.chdir(tmpdir)
            try:
                with mock.patch.object(reactor, "optimise", side_effect=RuntimeError("failed")):
                    with self.assertRaisesRegex(RuntimeError, "failed"):
                        reactor.relax_without_afir_bias(molecule, {"gamma": 100.0})
            finally:
                os.chdir(cwd)

        self.assertEqual(molecule.name, "candidate")

    def test_reactor_writes_disconnected_product_reference(self):
        import pyar.reaction_identity as reaction_identity

        molecule = Molecule(
            ["C", "H"],
            np.array([[0.0, 0.0, 0.0], [1.1, 0.0, 0.0]]),
            name="candidate",
            fragments=[[0], [1]],
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            output = os.path.join(tmpdir, "reference.xyz")
            reaction_identity.write_disconnected_reference(molecule, output)
            reference = Molecule.from_xyz(output)

        self.assertAlmostEqual(reference.coordinates[1, 0] - reference.coordinates[0, 0], 101.1)
        self.assertAlmostEqual(molecule.coordinates[1, 0] - molecule.coordinates[0, 0], 1.1)

    def test_reactor_build_gamma_schedule_returns_numeric_values(self):
        from pyar.workflows import reaction as reactor

        schedule = reactor.build_gamma_schedule(0.1, 0.2, steps=3)

        self.assertEqual(len(schedule), 3)
        self.assertAlmostEqual(float(schedule[0]), 0.1)
        self.assertAlmostEqual(float(schedule[1]), 0.15)
        self.assertAlmostEqual(float(schedule[2]), 0.2)
        self.assertTrue(all(isinstance(value, np.floating) or isinstance(value, float) for value in schedule))

    def test_reactor_rejects_invalid_gamma_schedule(self):
        from pyar.workflows import reaction as reactor

        with self.assertRaisesRegex(ValueError, "non-negative"):
            reactor.build_gamma_schedule(-1.0, 100.0)
        with self.assertRaisesRegex(ValueError, "greater than or equal"):
            reactor.build_gamma_schedule(100.0, 10.0)

    def test_reactor_equal_gamma_limits_run_one_cycle(self):
        from pyar.workflows import reaction as reactor

        schedule = reactor.build_gamma_schedule(100.0, 100.0)

        self.assertEqual(schedule.tolist(), [100.0])

    def test_reactor_format_gamma_id_uses_zero_padding(self):
        from pyar.workflows import reaction as reactor

        self.assertEqual(reactor.format_gamma_id(0.1), "0000p1")
        self.assertEqual(reactor.format_gamma_id(12.9), "0012p9")
        self.assertEqual(reactor.format_gamma_id(100.0), "0100")

    def test_reactor_initialize_reaction_run_prepares_fresh_state(self):
        from pyar.workflows import reaction as reactor

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
                    mock.patch.object(reactor.trial_generation, "create_trial_geometries", return_value=[orientation]), \
                    mock.patch.object(
                        reactor,
                        "separated_reactant_identity",
                        return_value={"inchi": "reactants-inchi", "smiles": "reactants-smiles"},
                    ):
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
        from pyar.workflows import reaction as reactor

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
        from pyar.workflows import reaction as reactor

        original_registry = dict(reactor.saved_product_identities)
        with tempfile.TemporaryDirectory() as tmpdir:
            reaction_dir = os.path.join(tmpdir, "reaction")
            os.makedirs(reaction_dir)
            cwd = os.getcwd()
            os.chdir(tmpdir)
            try:
                reactor.saved_product_identities.clear()
                run_state = mock.Mock()
                run_state.data = {"products": []}
                with mock.patch.object(
                    reactor,
                    "initialize_reaction_run",
                    return_value=(tmpdir, reaction_dir, run_state, [0.1], [], os.path.join(reaction_dir, "products")),
                ), \
                    mock.patch.object(reactor, "optimize_all", return_value=[]) as optimize_all, \
                    mock.patch.object(reactor.os, "chdir"):
                    result = reactor.react(
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
                reactor.saved_product_identities.clear()
                reactor.saved_product_identities.update(original_registry)

        self.assertGreaterEqual(optimize_all.call_count, 1)
        first_qc_params = optimize_all.call_args_list[0].args[-1]
        self.assertIsInstance(first_qc_params["gamma"], float)
        self.assertEqual(first_qc_params["gamma"], 0.1)
        run_state.complete_cycle.assert_called_once_with(0.1, [])
        run_state.finish.assert_called_once_with("completed_no_candidates")
        self.assertIsInstance(result, ReactionResult)
        self.assertEqual(result.workflow, "reaction")
        self.assertEqual(result.status, "completed_no_candidates")

    def test_reactor_reports_product_terminal_state(self):
        from pyar.workflows import reaction as reactor

        original_registry = dict(reactor.saved_product_identities)
        try:
            reactor.saved_product_identities.clear()
            reactor.saved_product_identities["product"] = {
                "inchi": "product-inchi",
                "smiles": "product-smiles",
            }
            run_state = mock.Mock()
            run_state.data = {"products": [{"path": "products/product.xyz"}]}
            with tempfile.TemporaryDirectory() as tmpdir:
                with mock.patch.object(
                    reactor,
                    "initialize_reaction_run",
                    return_value=(tmpdir, tmpdir, run_state, [100.0], [], tmpdir),
                ), \
                    mock.patch.object(reactor, "optimize_all", return_value=[]), \
                    mock.patch.object(reactor.os, "chdir"):
                    result = reactor.react(object(), object(), 100.0, 100.0, 1, {"software": "xtb"}, None, 2.3)

            run_state.complete_cycle.assert_called_once_with(100.0, [])
            run_state.finish.assert_called_once_with("completed_products_found")
            self.assertIsInstance(result, ReactionResult)
            self.assertEqual(result.status, "completed_products_found")
            self.assertEqual(result.selected_paths, ("products/product.xyz",))
        finally:
            reactor.saved_product_identities.clear()
            reactor.saved_product_identities.update(original_registry)

    def test_reactor_continues_gamma_schedule_when_products_exist_and_cycle_is_empty(self):
        reactor = importlib.import_module("pyar.workflows.reaction")

        original_registry = dict(reactor.saved_product_identities)
        try:
            reactor.saved_product_identities.clear()
            self.assertFalse(reactor._should_continue_after_product([]))
            reactor.saved_product_identities["product"] = {
                "inchi": "product-inchi",
                "smiles": "product-smiles",
            }
            self.assertTrue(reactor._should_continue_after_product([]))
            self.assertFalse(reactor._should_continue_after_product([object()]))
        finally:
            reactor.saved_product_identities.clear()
            reactor.saved_product_identities.update(original_registry)

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
