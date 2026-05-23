import os
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from pyar import aggregator


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
                with mock.patch.object(aggregator.tabu, "create_trial_geometries", return_value=[trial_a, trial_b]):
                    with mock.patch.object(aggregator, "optimise", side_effect=fake_optimize):
                        with mock.patch.object(aggregator.clustering, "choose_geometries", return_value=[trial_a, trial_b]) as chooser:
                            result = aggregator.add_one(
                                aggregate_id="ag_test",
                                seeds=[seed],
                                monomer=monomer,
                                hm_orientations=2,
                                qc_params=qc_params,
                                maximum_number_of_seeds=2,
                                tabu_on=True,
                                grid_on=True,
                                site=None,
                            )
            finally:
                os.chdir(cwd)

        self.assertEqual(thresholds, [None, None])
        chooser.assert_called_once()
        self.assertTrue(chooser.call_args.kwargs["persist_basin_memory"])
        self.assertEqual([m.name for m in result], ["trial_a", "trial_b"])

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
                qc_params = {"software": "xtb", "opt_threshold": "normal", "_two_layer_optimization": True}
                with mock.patch.object(aggregator.tabu, "create_trial_geometries", return_value=[trial_a, trial_b]):
                    with mock.patch.object(aggregator, "optimise", side_effect=fake_optimize):
                        with mock.patch.object(aggregator.shutil, "copy", return_value=None):
                            with mock.patch.object(aggregator.clustering, "choose_geometries", return_value=[trial_a]) as chooser:
                                result = aggregator.add_one(
                                    aggregate_id="ag_test",
                                seeds=[seed],
                                monomer=monomer,
                                hm_orientations=2,
                                qc_params=qc_params,
                                maximum_number_of_seeds=2,
                                tabu_on=True,
                                grid_on=True,
                                site=None,
                            )
                snapshot_exists = Path(tmpdir, "selected", "stoichiometry_H", "result_trial_a.xyz").is_file()
                registry_exists = Path(tmpdir, "selected", "stoichiometry_H", "basin_registry.json").is_file()
                nested_selected_exists = Path(tmpdir, "selected", "selected").exists()
            finally:
                os.chdir(cwd)

        self.assertEqual(thresholds, ["loose", "loose", "normal"])
        chooser.assert_called_once()
        self.assertFalse(chooser.call_args.kwargs["persist_basin_memory"])
        self.assertEqual([m.name for m in result], ["trial_a"])
        self.assertTrue(snapshot_exists)
        self.assertTrue(registry_exists)
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
                with mock.patch.object(aggregator.tabu, "create_trial_geometries", return_value=[trial_a, trial_b]):
                    with mock.patch.object(aggregator, "optimise", side_effect=fake_optimize):
                        with mock.patch.object(aggregator.shutil, "copy", return_value=None):
                            with mock.patch.object(aggregator.clustering, "choose_geometries", return_value=[trial_a, trial_b]):
                                result = aggregator.add_one(
                                    aggregate_id="ag_test",
                                    seeds=[seed],
                                    monomer=monomer,
                                    hm_orientations=2,
                                    qc_params=qc_params,
                                    maximum_number_of_seeds=2,
                                    tabu_on=True,
                                    grid_on=True,
                                    site=None,
                                )
            finally:
                os.chdir(cwd)

        self.assertEqual(refined_names, ["trial_a", "trial_b"])
        self.assertEqual([m.name for m in result], ["trial_b"])

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

    def test_selected_snapshot_is_grouped_by_stoichiometry(self):
        selected = [
            DummyMolecule("sel_1", n_atoms=5, atoms_list=["C", "H", "H", "H", "H"]),
            DummyMolecule("sel_2", n_atoms=5, atoms_list=["C", "H", "H", "H", "H"]),
        ]

        with tempfile.TemporaryDirectory() as tmpdir:
            cwd = os.getcwd()
            os.chdir(tmpdir)
            try:
                snapshot_dir = aggregator._snapshot_selected_geometries(selected)
                self.assertEqual(snapshot_dir, "selected/stoichiometry_CH4")
                self.assertTrue(Path(tmpdir, "selected", "stoichiometry_CH4", "result_sel_1.xyz").is_file())
                self.assertTrue(Path(tmpdir, "selected", "stoichiometry_CH4", "result_sel_2.xyz").is_file())
            finally:
                os.chdir(cwd)

    def test_aggregate_ignores_stale_restart_paths_without_aggregates_directory(self):
        molecules = [DummyMolecule("a", n_atoms=1), DummyMolecule("b", n_atoms=1)]

        with tempfile.TemporaryDirectory() as tmpdir:
            cwd = os.getcwd()
            os.chdir(tmpdir)
            try:
                with mock.patch.object(aggregator, "read_old_path", return_value=["a"]):
                    with mock.patch.object(aggregator, "old_path_to_new_path") as restart_paths:
                        with mock.patch.object(aggregator, "select_pathways", return_value=[[]]) as new_paths:
                            with self.assertLogs("pyar.aggregator", level="INFO") as captured:
                                aggregator.aggregate(
                                    molecules=molecules,
                                    aggregate_sizes=[1, 1],
                                    hm_orientations=2,
                                    qc_params={"software": None},
                                    maximum_number_of_seeds=2,
                                    first_pathway=0,
                                    number_of_pathways=1,
                                    tabu_on=True,
                                    grid_on=True,
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

    def test_aggregate_api_starts_without_existing_log(self):
        molecule = DummyMolecule("seed", n_atoms=1)

        with tempfile.TemporaryDirectory() as tmpdir:
            cwd = os.getcwd()
            os.chdir(tmpdir)
            try:
                aggregator.aggregate(
                    molecules=[molecule],
                    aggregate_sizes=[1],
                    hm_orientations=2,
                    qc_params={"software": None},
                    maximum_number_of_seeds=2,
                    first_pathway=0,
                    number_of_pathways=1,
                    tabu_on=True,
                    grid_on=True,
                    site=None,
                )

                self.assertTrue(Path(tmpdir, "aggregates").is_dir())
            finally:
                os.chdir(cwd)


if __name__ == "__main__":
    unittest.main()
