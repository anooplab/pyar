import json
import os
import importlib
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from pyar import aggregator
from pyar.aggregate_state import AggregateStateError
from pyar.workflows import _growth as growth
from pyar.workflow_results import AggregateResult

aggregate_workflow = importlib.import_module("pyar.workflows.aggregate")


class DummyMolecule:
    def __init__(self, name, n_atoms=2, atoms_list=None):
        self.name = name
        self.atoms_list = atoms_list if atoms_list is not None else ["H"] * n_atoms
        self.coordinates = [[0.0, 0.0, 0.0] for _ in range(n_atoms)]
        self.number_of_atoms = n_atoms
        self.charge = 0
        self.multiplicity = 1
        self.scftype = "rhf"
        self.atomic_number = [1] * n_atoms
        self.fragments = [0] * n_atoms

    def __len__(self):
        return len(self.atoms_list)

    def mol_to_xyz(self, filename):
        Path(filename).write_text(
            f"{len(self.atoms_list)}\n{self.name}\n"
            + "\n".join("H 0.0 0.0 0.0" for _ in self.atoms_list)
            + "\n"
        )


class AggregatorTests(unittest.TestCase):
    def test_single_stage_backend_skips_refinement(self):
        seed = DummyMolecule("seed", n_atoms=2)
        monomer = DummyMolecule("monomer", n_atoms=1)
        trial_a = DummyMolecule("trial_a", n_atoms=1)
        trial_a.coordinates = [[0.0, 0.0, 0.0]]
        trial_a.energy = 0.0
        trial_b = DummyMolecule("trial_b", n_atoms=1)
        trial_b.coordinates = [[1.0, 0.0, 0.0]]
        trial_b.energy = 1.0
        thresholds = []

        def fake_optimize(molecule, qc_params):
            thresholds.append(qc_params.get("opt_threshold"))
            return True

        with tempfile.TemporaryDirectory() as tmpdir:
            cwd = os.getcwd()
            os.chdir(tmpdir)
            try:
                qc_params = {"software": "aimnet_2", "opt_threshold": None, "_two_layer_optimization": False}
                with mock.patch.object(growth.trial_generation, "create_trial_geometries", return_value=[trial_a, trial_b]):
                    with mock.patch.object(growth, "optimise", side_effect=fake_optimize):
                        with mock.patch.object(growth.clustering, "choose_geometries", return_value=[trial_a, trial_b]) as chooser:
                            result = aggregator.add_one(
                                aggregate_id="ag_test",
                                seeds=[seed],
                                monomer=monomer,
                                hm_orientations=2,
                                qc_params=qc_params,
                                maximum_number_of_seeds=2,
                                site=None,
                            )
                snapshot_exists = Path(tmpdir, "selected", "result_trial_a.xyz").is_file()
                second_snapshot_exists = Path(tmpdir, "selected", "result_trial_b.xyz").is_file()
                nested_selected_exists = Path(tmpdir, "selected", "stoichiometry_H").exists()
            finally:
                os.chdir(cwd)

        self.assertEqual(thresholds, [None, None])
        chooser.assert_called_once()
        self.assertTrue(chooser.call_args.kwargs["persist_basin_memory"])
        self.assertFalse(chooser.call_args.kwargs["group_basin_by_stoichiometry"])
        self.assertEqual([m.name for m in result], ["trial_a", "trial_b"])
        self.assertTrue(snapshot_exists)
        self.assertTrue(second_snapshot_exists)
        self.assertFalse(nested_selected_exists)

    def test_two_stage_backend_uses_loose_then_normal(self):
        seed = DummyMolecule("seed", n_atoms=2)
        monomer = DummyMolecule("monomer", n_atoms=1)
        trial_a = DummyMolecule("trial_a", n_atoms=1)
        trial_a.coordinates = [[0.0, 0.0, 0.0]]
        trial_a.energy = 0.0
        trial_b = DummyMolecule("trial_b", n_atoms=1)
        trial_b.coordinates = [[1.0, 0.0, 0.0]]
        trial_b.energy = 1.0
        thresholds = []

        def fake_optimize(molecule, qc_params):
            thresholds.append(qc_params.get("opt_threshold"))
            return True

        with tempfile.TemporaryDirectory() as tmpdir:
            cwd = os.getcwd()
            os.chdir(tmpdir)
            try:
                qc_params = {"software": "xtb", "opt_threshold": "tight", "_two_layer_optimization": True}
                with mock.patch.object(growth.trial_generation, "create_trial_geometries", return_value=[trial_a, trial_b]):
                    with mock.patch.object(growth, "optimise", side_effect=fake_optimize):
                        with mock.patch.object(growth.shutil, "copy", return_value=None):
                            with mock.patch.object(growth.clustering, "choose_geometries", return_value=[trial_a]) as chooser:
                                result = aggregator.add_one(
                                    aggregate_id="ag_test",
                                    seeds=[seed],
                                    monomer=monomer,
                                    hm_orientations=2,
                                    qc_params=qc_params,
                                    maximum_number_of_seeds=2,
                                    site=None,
                                )
                snapshot_exists = Path(tmpdir, "selected", "result_trial_a.xyz").is_file()
                registry_exists = Path(tmpdir, "selected", "basin_registry.json").is_file()
                nested_selected_exists = Path(tmpdir, "selected", "stoichiometry_H").exists()
            finally:
                os.chdir(cwd)

        self.assertEqual(thresholds, ["loose", "loose", "normal"])
        self.assertEqual(qc_params["opt_threshold"], "tight")
        chooser.assert_called_once()
        self.assertFalse(chooser.call_args.kwargs["persist_basin_memory"])
        self.assertFalse(chooser.call_args.kwargs["group_basin_by_stoichiometry"])
        self.assertEqual([m.name for m in result], ["trial_a"])
        self.assertTrue(snapshot_exists)
        self.assertFalse(registry_exists)
        self.assertFalse(nested_selected_exists)

    def test_two_stage_refinement_does_not_skip_after_failure(self):
        seed = DummyMolecule("seed", n_atoms=2)
        monomer = DummyMolecule("monomer", n_atoms=1)
        trial_a = DummyMolecule("trial_a", n_atoms=1)
        trial_a.coordinates = [[0.0, 0.0, 0.0]]
        trial_a.energy = 0.0
        trial_b = DummyMolecule("trial_b", n_atoms=1)
        trial_b.coordinates = [[1.0, 0.0, 0.0]]
        trial_b.energy = 1.0
        refined_names = []

        def fake_optimize(molecule, qc_params):
            if qc_params.get("opt_threshold") == "normal":
                refined_names.append(molecule.name)
                return molecule.name == "trial_b"
            return True

        with tempfile.TemporaryDirectory() as tmpdir:
            cwd = os.getcwd()
            os.chdir(tmpdir)
            try:
                qc_params = {"software": "xtb", "opt_threshold": "normal", "_two_layer_optimization": True}
                with mock.patch.object(growth.trial_generation, "create_trial_geometries", return_value=[trial_a, trial_b]):
                    with mock.patch.object(growth, "optimise", side_effect=fake_optimize):
                        with mock.patch.object(growth.shutil, "copy", return_value=None):
                            with mock.patch.object(growth.clustering, "choose_geometries", return_value=[trial_a, trial_b]):
                                result = aggregator.add_one(
                                    aggregate_id="ag_test",
                                    seeds=[seed],
                                    monomer=monomer,
                                    hm_orientations=2,
                                    qc_params=qc_params,
                                    maximum_number_of_seeds=2,
                                    site=None,
                                )
            finally:
                os.chdir(cwd)

        self.assertEqual(refined_names, ["trial_a", "trial_b"])
        self.assertEqual([m.name for m in result], ["trial_b"])

    def test_add_one_restores_cwd_when_seed_optimization_raises(self):
        seed = DummyMolecule("seed", n_atoms=2)
        monomer = DummyMolecule("monomer", n_atoms=1)
        trial = DummyMolecule("trial", n_atoms=1)

        with tempfile.TemporaryDirectory() as tmpdir:
            cwd = os.getcwd()
            os.chdir(tmpdir)
            try:
                qc_params = {"software": "xtb", "_two_layer_optimization": True}
                with mock.patch.object(growth.trial_generation, "create_trial_geometries", return_value=[trial]):
                    with mock.patch.object(growth, "optimise", side_effect=RuntimeError("seed optimization failed")):
                        with self.assertRaisesRegex(RuntimeError, "seed optimization failed"):
                            aggregator.add_one(
                                aggregate_id="ag_test",
                                seeds=[seed],
                                monomer=monomer,
                                hm_orientations=1,
                                qc_params=qc_params,
                                maximum_number_of_seeds=1,
                                site=None,
                            )
                self.assertEqual(os.getcwd(), tmpdir)
            finally:
                os.chdir(cwd)

    def test_add_one_restores_cwd_when_refinement_raises(self):
        seed = DummyMolecule("seed", n_atoms=2)
        monomer = DummyMolecule("monomer", n_atoms=1)
        trial = DummyMolecule("trial", n_atoms=1)
        trial.energy = 0.0

        def fail_during_refinement(molecule, qc_params):
            if qc_params.get("opt_threshold") == "normal":
                raise RuntimeError("refinement failed")
            return True

        with tempfile.TemporaryDirectory() as tmpdir:
            cwd = os.getcwd()
            os.chdir(tmpdir)
            try:
                qc_params = {"software": "xtb", "_two_layer_optimization": True}
                with mock.patch.object(growth.trial_generation, "create_trial_geometries", return_value=[trial]):
                    with mock.patch.object(growth, "optimise", side_effect=fail_during_refinement):
                        with mock.patch.object(growth.clustering, "choose_geometries", return_value=[trial]):
                            with self.assertRaisesRegex(RuntimeError, "refinement failed"):
                                aggregator.add_one(
                                    aggregate_id="ag_test",
                                    seeds=[seed],
                                    monomer=monomer,
                                    hm_orientations=1,
                                    qc_params=qc_params,
                                    maximum_number_of_seeds=1,
                                    site=None,
                                )
                self.assertEqual(os.getcwd(), tmpdir)
            finally:
                os.chdir(cwd)

    def test_geometry_only_add_one_restores_cwd_when_generation_raises(self):
        seed = DummyMolecule("seed", n_atoms=2)
        monomer = DummyMolecule("monomer", n_atoms=1)

        with tempfile.TemporaryDirectory() as tmpdir:
            cwd = os.getcwd()
            os.chdir(tmpdir)
            try:
                with mock.patch.object(
                    growth.trial_generation,
                    "create_trial_geometries",
                    side_effect=RuntimeError("orientation generation failed"),
                ):
                    with self.assertRaisesRegex(RuntimeError, "orientation generation failed"):
                        aggregator.add_one(
                            aggregate_id="ag_test",
                            seeds=[seed],
                            monomer=monomer,
                            hm_orientations=1,
                            qc_params={},
                            maximum_number_of_seeds=1,
                            site=None,
                        )
                self.assertEqual(os.getcwd(), tmpdir)
            finally:
                os.chdir(cwd)

    def test_restart_finished_jobs_mutates_remaining_molecules_safely(self):
        done = DummyMolecule("done", n_atoms=1)
        done.energy = 0.0
        pending = DummyMolecule("pending", n_atoms=1)
        pending.energy = 1.0
        molecules = [done, pending]

        with tempfile.TemporaryDirectory() as tmpdir:
            cwd = os.getcwd()
            os.chdir(tmpdir)
            try:
                Path("selected", "job_done").mkdir(parents=True)
                Path("selected", "result_done.xyz").write_text("1\ndone\nH 0.0 0.0 0.0\n")
                result = aggregator.check_for_the_finished_jobs_on_restart(molecules, tmpdir)
                self.assertEqual(os.getcwd(), tmpdir)
            finally:
                os.chdir(cwd)

        self.assertEqual([m.name for m in result], ["done"])
        self.assertEqual([m.name for m in molecules], ["pending"])

    def test_legacy_path_conversion_preserves_repeated_fragment_counts(self):
        fragment_a = DummyMolecule("a", n_atoms=1)
        fragment_b = DummyMolecule("b", n_atoms=1)

        pathways = aggregator.old_path_to_new_path(
            [fragment_a, fragment_b, fragment_b],
            ["  abb  "],
        )

        self.assertEqual([[molecule.name for molecule in path] for path in pathways], [["a", "b", "b"]])

    def test_selected_snapshot_is_grouped_by_stoichiometry(self):
        selected = [
            DummyMolecule("sel_1", n_atoms=5, atoms_list=["C", "H", "H", "H", "H"]),
            DummyMolecule("sel_2", n_atoms=5, atoms_list=["C", "H", "H", "H", "H"]),
        ]

        with tempfile.TemporaryDirectory() as tmpdir:
            cwd = os.getcwd()
            os.chdir(tmpdir)
            try:
                snapshot_dir = aggregator._snapshot_selected_geometries(
                    selected,
                    group_by_stoichiometry=True,
                )
                self.assertEqual(snapshot_dir, "selected/stoichiometry_CH4")
                self.assertTrue(Path(tmpdir, "selected", "stoichiometry_CH4", "result_sel_1.xyz").is_file())
                self.assertTrue(Path(tmpdir, "selected", "stoichiometry_CH4", "result_sel_2.xyz").is_file())
            finally:
                os.chdir(cwd)

    def test_pathway_snapshot_replaces_obsolete_flat_results(self):
        previous = DummyMolecule("previous", n_atoms=2)
        current = DummyMolecule("current", n_atoms=2)

        with tempfile.TemporaryDirectory() as tmpdir:
            cwd = os.getcwd()
            os.chdir(tmpdir)
            try:
                aggregator._snapshot_selected_geometries(
                    [previous],
                    output_root="selected",
                    group_by_stoichiometry=False,
                )
                aggregator._snapshot_selected_geometries(
                    [current],
                    output_root="selected",
                    group_by_stoichiometry=False,
                )
                result_files = sorted(path.name for path in Path("selected").glob("result_*.xyz"))
            finally:
                os.chdir(cwd)

        self.assertEqual(result_files, ["result_current.xyz"])

    def test_final_selected_geometries_cluster_across_pathways(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cwd = os.getcwd()
            os.chdir(tmpdir)
            try:
                Path("ag_a_001_000", "selected").mkdir(parents=True)
                Path("ag_a_002_000", "selected").mkdir(parents=True)
                Path("ag_a_001_000", "selected", "result_one.xyz").write_text(
                    "5\none: -1.0\n"
                    "C 0 0 0\n"
                    "H 0 0 1\n"
                    "H 0 1 0\n"
                    "H 1 0 0\n"
                    "H 0 0 -1\n"
                )
                Path("ag_a_002_000", "selected", "result_two.xyz").write_text(
                    "5\ntwo: -2.0\n"
                    "C 0 0 0\n"
                    "H 0 0 1\n"
                    "H 0 1 0\n"
                    "H 1 0 0\n"
                    "H 0 0 -1\n"
                )
                Path("ag_a_002_000", "selected", "stoichiometry_CH4").mkdir()
                Path(
                    "ag_a_002_000", "selected", "stoichiometry_CH4", "result_old_layout.xyz"
                ).write_text(
                    "5\nold layout: -10.0\n"
                    "C 0 0 0\n"
                    "H 0 0 1\n"
                    "H 0 1 0\n"
                    "H 1 0 0\n"
                    "H 0 0 -1\n"
                )

                with mock.patch.object(
                    growth.clustering,
                    "choose_geometries",
                    side_effect=lambda molecules, **kwargs: [molecules[0]],
                ) as chooser:
                    stale_final_dir = Path(tmpdir, "selected", "stoichiometry_CH4")
                    stale_final_dir.mkdir(parents=True)
                    stale_result = stale_final_dir / "result_stale.xyz"
                    stale_result.write_text("1\nstale: -3.0\nH 0 0 0\n")
                    selected = aggregator._finalize_selected_geometries(
                        aggregate_root='.',
                        maximum_number_of_seeds=1,
                        algorithm='maxmin',
                    )
                final_snapshot_exists = Path(tmpdir, "selected", "stoichiometry_CH4").exists()
                final_readme_exists = Path(tmpdir, "selected", "stoichiometry_CH4", "README.txt").is_file()
                stale_final_result_exists = stale_result.exists()
            finally:
                os.chdir(cwd)

        self.assertEqual(len(selected), 1)
        self.assertTrue(final_snapshot_exists)
        self.assertTrue(final_readme_exists)
        self.assertFalse(stale_final_result_exists)
        self.assertFalse(chooser.call_args.kwargs["apply_basin_memory"])
        self.assertEqual(len(chooser.call_args.args[0]), 2)

    def test_aggregate_ignores_stale_restart_paths_without_aggregates_directory(self):
        molecules = [DummyMolecule("a", n_atoms=1), DummyMolecule("b", n_atoms=1)]

        with tempfile.TemporaryDirectory() as tmpdir:
            cwd = os.getcwd()
            os.chdir(tmpdir)
            try:
                with mock.patch.object(aggregate_workflow, "read_old_path", return_value=["a"]):
                    with mock.patch.object(aggregate_workflow, "old_path_to_new_path") as restart_paths:
                        with mock.patch.object(aggregate_workflow, "select_pathways", return_value=[[]]) as new_paths:
                            with mock.patch.object(aggregate_workflow, "add_one", return_value=[]):
                                with self.assertLogs("pyar.aggregator", level="INFO") as captured:
                                    with self.assertWarnsRegex(DeprecationWarning, "pyar.aggregator.aggregate"):
                                        aggregator.aggregate(
                                            molecules=molecules,
                                            aggregate_sizes=[1, 1],
                                            hm_orientations=2,
                                            qc_params={"software": None},
                                            maximum_number_of_seeds=2,
                                            first_pathway=0,
                                            number_of_pathways=1,
                                            site=None,
                                        )

                self.assertTrue(Path(tmpdir, "aggregates").is_dir())
                restart_paths.assert_not_called()
                new_paths.assert_called_once()
                self.assertTrue(
                    any("Ignoring restart markers" in message for message in captured.output)
                )
            finally:
                os.chdir(cwd)

    def test_aggregate_imports_legacy_path_markers_into_state(self):
        molecules = [DummyMolecule("input_a", n_atoms=1), DummyMolecule("input_b", n_atoms=1)]

        with tempfile.TemporaryDirectory() as tmpdir:
            cwd = os.getcwd()
            os.chdir(tmpdir)
            try:
                Path("aggregates").mkdir()
                with mock.patch.object(aggregate_workflow, "read_old_path", return_value=["  ab  "]):
                    with mock.patch.object(aggregate_workflow, "select_pathways") as new_paths:
                        with mock.patch.object(aggregate_workflow, "add_one", return_value=[]):
                            with self.assertLogs("pyar.aggregator", level="WARNING"):
                                aggregate_workflow.aggregate(
                                    molecules=molecules,
                                    aggregate_sizes=[1, 1],
                                    hm_orientations=2,
                                    qc_params={"software": None},
                                    maximum_number_of_seeds=2,
                                    first_pathway=0,
                                    number_of_pathways=1,
                                    site=None,
                                )
                with Path("aggregates", "state.json").open() as fp:
                    state = json.load(fp)
            finally:
                os.chdir(cwd)

        new_paths.assert_not_called()
        self.assertEqual(state["legacy_import"], "pyar.log")
        self.assertEqual(state["pathways"], ["ab"])
        self.assertEqual(state["status"], "completed")

    def test_aggregate_api_starts_without_existing_log(self):
        molecule = DummyMolecule("seed", n_atoms=1)

        with tempfile.TemporaryDirectory() as tmpdir:
            cwd = os.getcwd()
            os.chdir(tmpdir)
            try:
                with self.assertWarnsRegex(DeprecationWarning, "pyar.aggregator.aggregate"):
                    result = aggregator.aggregate(
                        molecules=[molecule],
                        aggregate_sizes=[1],
                        hm_orientations=2,
                        qc_params={"software": None},
                        maximum_number_of_seeds=2,
                        first_pathway=0,
                        number_of_pathways=1,
                        site=None,
                    )

                self.assertTrue(Path(tmpdir, "aggregates").is_dir())
                with Path(tmpdir, "aggregates", "state.json").open() as fp:
                    state = json.load(fp)
            finally:
                os.chdir(cwd)

        self.assertIsInstance(result, AggregateResult)
        self.assertEqual(result.workflow, "aggregate")
        self.assertEqual(result.status, "completed")
        self.assertTrue(result.state_path.endswith("aggregates/state.json"))
        self.assertEqual(state["workflow"], "aggregate")
        self.assertEqual(state["status"], "completed")
        self.assertEqual(state["request"]["fragments"][0]["scftype"], "rhf")
        self.assertEqual(state["request"]["fragments"][0]["fragment_definition"], [0])

    def test_aggregate_refuses_existing_output_without_resumable_state(self):
        molecule = DummyMolecule("seed", n_atoms=1)

        with tempfile.TemporaryDirectory() as tmpdir:
            cwd = os.getcwd()
            os.chdir(tmpdir)
            try:
                Path("aggregates", "archived_output").mkdir(parents=True)
                with mock.patch.object(aggregate_workflow, "read_old_path", return_value=None):
                    with self.assertRaisesRegex(AggregateStateError, "has no resumable state"):
                        aggregate_workflow.aggregate(
                            molecules=[molecule],
                            aggregate_sizes=[1],
                            hm_orientations=2,
                            qc_params={"software": None},
                            maximum_number_of_seeds=2,
                            first_pathway=0,
                            number_of_pathways=1,
                            site=None,
                        )
            finally:
                os.chdir(cwd)

    def test_aggregate_resumes_only_uncompleted_persisted_pathways(self):
        molecule_a = DummyMolecule("input_a", n_atoms=1)
        molecule_b = DummyMolecule("input_b", n_atoms=1)
        calls = []

        def stop_on_second_path(*args, **kwargs):
            calls.append(args[0])
            if len(calls) == 2:
                raise RuntimeError("interrupted second pathway")
            return []

        with tempfile.TemporaryDirectory() as tmpdir:
            cwd = os.getcwd()
            os.chdir(tmpdir)
            try:
                with mock.patch.object(
                    aggregate_workflow,
                    "select_pathways",
                    return_value=[(molecule_a, molecule_b), (molecule_b, molecule_a)],
                ):
                    with mock.patch.object(aggregate_workflow, "add_one", side_effect=stop_on_second_path):
                        with self.assertRaisesRegex(RuntimeError, "interrupted second pathway"):
                            aggregate_workflow.aggregate(
                                molecules=[molecule_a, molecule_b],
                                aggregate_sizes=[1, 1],
                                hm_orientations=2,
                                qc_params={"software": None},
                                maximum_number_of_seeds=2,
                                first_pathway=0,
                                number_of_pathways=2,
                                site=None,
                            )
                with Path("aggregates", "state.json").open() as fp:
                    interrupted_state = json.load(fp)

                with mock.patch.object(aggregate_workflow, "select_pathways") as selected_paths:
                    with mock.patch.object(aggregate_workflow, "add_one", return_value=[]) as resumed_add:
                        aggregate_workflow.aggregate(
                            molecules=[molecule_a, molecule_b],
                            aggregate_sizes=[1, 1],
                            hm_orientations=2,
                            qc_params={"software": None},
                            maximum_number_of_seeds=2,
                            first_pathway=0,
                            number_of_pathways=2,
                            site=None,
                        )
                with Path("aggregates", "state.json").open() as fp:
                    completed_state = json.load(fp)
            finally:
                os.chdir(cwd)

        self.assertEqual(interrupted_state["pathways"], ["ab", "ba"])
        self.assertEqual(len(interrupted_state["completed_pathways"]), 1)
        selected_paths.assert_not_called()
        resumed_add.assert_called_once()
        self.assertEqual(completed_state["status"], "completed")
        self.assertEqual(len(completed_state["completed_pathways"]), 2)


if __name__ == "__main__":
    unittest.main()
