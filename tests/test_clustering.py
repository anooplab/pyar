import json
import os
import tempfile
import unittest
from types import SimpleNamespace
from unittest import mock

import numpy as np

from pyar.data_analysis import clustering


class ClusteringTests(unittest.TestCase):
    def setUp(self):
        clustering._MBTR_RUNTIME_DISABLED = False
        clustering._MBTR_DISABLE_REASON = None

    def test_choose_geometries_falls_back_when_mbtr_fails(self):
        molecules = [
            SimpleNamespace(name="a", atoms_list=["C"], coordinates=[[0.0, 0.0, 0.0]], energy=0.0),
            SimpleNamespace(name="b", atoms_list=["H"], coordinates=[[1.0, 0.0, 0.0]], energy=1.0),
        ]

        with mock.patch("pyar.representations.mbtr_descriptor", side_effect=TypeError("boom")):
            with mock.patch("pyar.data_analysis.clustering.remove_similar", return_value=molecules) as remover:
                result = clustering.choose_geometries(molecules, maximum_number_of_seeds=2)

        remover.assert_called_once_with(molecules)
        self.assertEqual(result, molecules)

    def test_choose_geometries_returns_available_set_when_underfull(self):
        molecules = [
            SimpleNamespace(name="a", atoms_list=["H"], coordinates=[[0.0, 0.0, 0.0]], energy=0.0),
        ]

        result = clustering.choose_geometries(molecules, maximum_number_of_seeds=4)

        self.assertEqual(result, molecules)

    def test_mbtr_is_disabled_after_first_runtime_failure(self):
        molecules = [
            SimpleNamespace(name="a", atoms_list=["C"], coordinates=[[0.0, 0.0, 0.0]], energy=0.0),
            SimpleNamespace(name="b", atoms_list=["H"], coordinates=[[1.0, 0.0, 0.0]], energy=1.0),
        ]

        with mock.patch("pyar.representations.mbtr_descriptor", side_effect=TypeError("mbtr boom")) as mbtr_mock:
            with mock.patch("pyar.data_analysis.clustering.remove_similar", return_value=molecules):
                clustering.choose_geometries(molecules, maximum_number_of_seeds=1)

        self.assertEqual(mbtr_mock.call_count, 1)

        with mock.patch("pyar.representations.mbtr_descriptor") as mbtr_mock:
            with mock.patch("pyar.data_analysis.clustering.remove_similar", return_value=molecules):
                result = clustering.choose_geometries(molecules, maximum_number_of_seeds=1)

        mbtr_mock.assert_not_called()
        self.assertEqual(len(result), 1)

    def test_fallback_enforces_maximum_seed_limit(self):
        molecules = [
            SimpleNamespace(name=f"m{i}", atoms_list=["H"], coordinates=[[float(i), 0.0, 0.0]], energy=float(i))
            for i in range(5)
        ]

        with mock.patch("pyar.representations.mbtr_descriptor", side_effect=TypeError("mbtr boom")):
            with mock.patch("pyar.data_analysis.clustering.remove_similar", return_value=molecules):
                result = clustering.choose_geometries(molecules, maximum_number_of_seeds=2)

        self.assertEqual(len(result), 2)
        self.assertEqual([m.name for m in result], ["m0", "m1"])

    def test_algorithm_failure_enforces_maximum_seed_limit(self):
        molecules = [
            SimpleNamespace(name=f"m{i}", atoms_list=["H"], coordinates=[[float(i), 0.0, 0.0]], energy=float(i))
            for i in range(4)
        ]

        with self.assertLogs("pyar.cluster", level="WARNING") as captured:
            with mock.patch("pyar.representations.mbtr_descriptor", return_value=[0.0]):
                with mock.patch("pyar.data_analysis.clustering.generate_labels", side_effect=RuntimeError("cluster fail")):
                    result = clustering.choose_geometries(molecules, maximum_number_of_seeds=2)

        self.assertEqual(len(result), 2)
        self.assertEqual([m.name for m in result], ["m0", "m1"])
        fallback_messages = [
            message for message in captured.output
            if "Clustering algorithm hybrid failed" in message
        ]
        self.assertEqual(len(fallback_messages), 1)
        self.assertNotIn("Traceback", fallback_messages[0])

    def test_algorithm_failure_persists_basin_registry(self):
        molecules = [
            SimpleNamespace(name=f"m{i}", atoms_list=["H"], coordinates=[[float(i), 0.0, 0.0]], energy=float(i))
            for i in range(4)
        ]

        with tempfile.TemporaryDirectory() as tmpdir:
            cwd = os.getcwd()
            os.chdir(tmpdir)
            try:
                os.mkdir("selected")
                with mock.patch("pyar.representations.mbtr_descriptor", return_value=[0.0]):
                    with mock.patch("pyar.data_analysis.clustering.generate_labels", side_effect=RuntimeError("cluster fail")):
                        clustering.choose_geometries(molecules, maximum_number_of_seeds=2)

                registry_path = os.path.join("selected", "stoichiometry_H", "basin_registry.json")
                self.assertTrue(os.path.exists(registry_path))
            finally:
                os.chdir(cwd)

    def test_persist_basin_memory_false_does_not_write_registry(self):
        molecules = [
            SimpleNamespace(name="m0", atoms_list=["H"], coordinates=[[0.0, 0.0, 0.0]], energy=0.0),
        ]

        with tempfile.TemporaryDirectory() as tmpdir:
            cwd = os.getcwd()
            os.chdir(tmpdir)
            try:
                os.mkdir("selected")
                clustering.choose_geometries(
                    molecules,
                    maximum_number_of_seeds=2,
                    persist_basin_memory=False,
                )
                registry_path = os.path.join("selected", "stoichiometry_H", "basin_registry.json")
                self.assertFalse(os.path.exists(registry_path))
            finally:
                os.chdir(cwd)

    def test_remove_similar_handles_overlapping_atoms(self):
        molecules = [
            SimpleNamespace(name="a", atoms_list=["H", "H"], coordinates=[[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]], energy=0.0),
            SimpleNamespace(name="b", atoms_list=["H", "H"], coordinates=[[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]], energy=1.0),
        ]

        result = clustering.remove_similar(molecules)

        self.assertGreaterEqual(len(result), 1)

    def test_remove_similar_keeps_lower_energy_duplicate(self):
        molecules = [
            SimpleNamespace(
                name="high",
                atoms_list=["H", "H"],
                coordinates=[[1.0, 0.0, 0.0], [1.5, 0.0, 0.0]],
                energy=1e-6,
            ),
            SimpleNamespace(
                name="low",
                atoms_list=["H", "H"],
                coordinates=[[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]],
                energy=0.0,
            ),
        ]

        with mock.patch("pyar.data_analysis.clustering.calc_fingerprint_distance", return_value=0.5):
            result = clustering.remove_similar(molecules)

        self.assertEqual([m.name for m in result], ["low"])

    def test_remove_similar_keeps_same_energy_distinct_structures(self):
        molecules = [
            SimpleNamespace(name="a", atoms_list=["H"], coordinates=[[0.0, 0.0, 0.0]], energy=0.0),
            SimpleNamespace(name="b", atoms_list=["H"], coordinates=[[10.0, 0.0, 0.0]], energy=0.0),
        ]

        with mock.patch("pyar.data_analysis.clustering.calc_fingerprint_distance", return_value=5.0):
            result = clustering.remove_similar(molecules)

        self.assertEqual([m.name for m in result], ["a", "b"])

    def test_choose_geometries_discards_disconnected_candidates_when_connected_options_exist(self):
        molecules = [
            SimpleNamespace(
                name="connected",
                atoms_list=["H", "H"],
                coordinates=[[0.0, 0.0, 0.0], [0.7, 0.0, 0.0]],
                energy=0.0,
            ),
            SimpleNamespace(
                name="broken",
                atoms_list=["H", "H"],
                coordinates=[[0.0, 0.0, 0.0], [5.0, 0.0, 0.0]],
                energy=0.1,
            ),
        ]

        with mock.patch("pyar.data_analysis.clustering.remove_similar", return_value=molecules):
            with mock.patch("pyar.tabu.broken", side_effect=lambda mol: mol.name == "broken"):
                result = clustering.choose_geometries(molecules, maximum_number_of_seeds=1)

        self.assertEqual([m.name for m in result], ["connected"])

    def test_rmsd_duplicate_test_detects_same_geometry_under_translation(self):
        a = SimpleNamespace(
            name="a",
            atoms_list=["H", "H"],
            coordinates=[[0.0, 0.0, 0.0], [0.7, 0.0, 0.0]],
            energy=0.0,
        )
        b = SimpleNamespace(
            name="b",
            atoms_list=["H", "H"],
            coordinates=[[5.0, -2.0, 1.0], [5.7, -2.0, 1.0]],
            energy=0.0,
        )

        with mock.patch("pyar.data_analysis.clustering.calc_fingerprint_distance", return_value=0.5):
            self.assertLess(clustering._rmsd_after_alignment(a, b), 1e-6)

    def test_rmsd_duplicate_test_detects_same_geometry_under_rotation(self):
        a = SimpleNamespace(
            name="a",
            atoms_list=["H", "H"],
            coordinates=[[0.0, 0.0, 0.0], [0.7, 0.0, 0.0]],
            energy=0.0,
        )
        b = SimpleNamespace(
            name="b",
            atoms_list=["H", "H"],
            coordinates=[[0.0, 0.0, 0.0], [0.0, 0.7, 0.0]],
            energy=0.0,
        )

        self.assertLess(clustering._rmsd_after_alignment(a, b), 1e-6)

    def test_rmsd_duplicate_test_matches_equivalent_atoms_after_reordering(self):
        a = SimpleNamespace(
            name="a",
            atoms_list=["C", "H", "H", "H", "H"],
            coordinates=[
                [0.0, 0.0, 0.0],
                [1.0, 1.0, 1.0],
                [1.0, -1.0, -1.0],
                [-1.0, 1.0, -1.0],
                [-1.0, -1.0, 1.0],
            ],
            energy=0.0,
        )
        b = SimpleNamespace(
            name="b",
            atoms_list=["H", "C", "H", "H", "H"],
            coordinates=[
                [-1.0, -1.0, -1.0],
                [0.0, 0.0, 0.0],
                [-1.0, 1.0, 1.0],
                [1.0, -1.0, 1.0],
                [1.0, 1.0, -1.0],
            ],
            energy=0.0,
        )

        self.assertLess(clustering._rmsd_after_alignment(a, b), 1e-6)

    def test_rmsd_duplicate_test_keeps_changed_bond_length(self):
        a = SimpleNamespace(
            name="a",
            atoms_list=["H", "H"],
            coordinates=[[0.0, 0.0, 0.0], [0.7, 0.0, 0.0]],
            energy=0.0,
        )
        b = SimpleNamespace(
            name="b",
            atoms_list=["H", "H"],
            coordinates=[[0.0, 0.0, 0.0], [1.2, 0.0, 0.0]],
            energy=0.0,
        )

        self.assertGreater(clustering._rmsd_after_alignment(a, b), 0.1)

    def test_rmsd_duplicate_test_rejects_different_compositions(self):
        a = SimpleNamespace(
            name="a",
            atoms_list=["C", "H"],
            coordinates=[[0.0, 0.0, 0.0], [0.7, 0.0, 0.0]],
            energy=0.0,
        )
        b = SimpleNamespace(
            name="b",
            atoms_list=["N", "H"],
            coordinates=[[0.0, 0.0, 0.0], [0.7, 0.0, 0.0]],
            energy=0.0,
        )

        self.assertEqual(clustering._rmsd_after_alignment(a, b), float("inf"))

    def test_adaptive_duplicate_threshold_stays_conservative(self):
        molecules = [
            SimpleNamespace(name="a", atoms_list=["H", "H"], coordinates=[[0.0, 0.0, 0.0], [0.7, 0.0, 0.0]], energy=0.0),
            SimpleNamespace(name="b", atoms_list=["H", "H"], coordinates=[[0.0, 0.0, 0.0], [0.71, 0.0, 0.0]], energy=0.1),
            SimpleNamespace(name="c", atoms_list=["H", "H"], coordinates=[[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]], energy=0.2),
        ]

        with mock.patch("pyar.data_analysis.clustering.calc_fingerprint_distance", return_value=0.5):
            threshold = clustering._adaptive_duplicate_rmsd_threshold(molecules)

        self.assertGreaterEqual(threshold, 0.05)
        self.assertLessEqual(threshold, 0.15)

    def test_maxmin_selects_diverse_seed_set(self):
        molecules = [
            SimpleNamespace(name="m0", atoms_list=["H"], coordinates=[[0.0, 0.0, 0.0]], energy=0.0),
            SimpleNamespace(name="m1", atoms_list=["H"], coordinates=[[1.0, 0.0, 0.0]], energy=1.0),
            SimpleNamespace(name="m2", atoms_list=["H"], coordinates=[[10.0, 0.0, 0.0]], energy=2.0),
        ]

        def descriptor_from_coordinate(_atoms, coordinates):
            return [coordinates[0][0]]

        with mock.patch.dict("os.environ", {"PYAR_CLUSTERING_ALGORITHM": "maxmin"}):
            with mock.patch("pyar.representations.mbtr_descriptor", side_effect=descriptor_from_coordinate):
                result = clustering.choose_geometries(molecules, maximum_number_of_seeds=2)

        self.assertEqual(len(result), 2)
        self.assertEqual([m.name for m in result], ["m0", "m2"])

    def test_maxmin_prunes_similar_molecules_before_descriptor_selection(self):
        molecules = [
            SimpleNamespace(name="duplicate", atoms_list=["H"], coordinates=[[0.1, 0.0, 0.0]], energy=0.1),
            SimpleNamespace(name="m0", atoms_list=["H"], coordinates=[[0.0, 0.0, 0.0]], energy=0.0),
            SimpleNamespace(name="m1", atoms_list=["H"], coordinates=[[1.0, 0.0, 0.0]], energy=1.0),
            SimpleNamespace(name="m2", atoms_list=["H"], coordinates=[[10.0, 0.0, 0.0]], energy=2.0),
        ]
        pruned = molecules[1:]
        descriptor_calls = []

        def descriptor_from_coordinate(_atoms, coordinates):
            descriptor_calls.append(coordinates[0][0])
            return [coordinates[0][0]]

        with mock.patch.dict("os.environ", {"PYAR_CLUSTERING_ALGORITHM": "maxmin"}):
            with mock.patch("pyar.data_analysis.clustering.remove_similar", return_value=pruned):
                with mock.patch("pyar.representations.mbtr_descriptor", side_effect=descriptor_from_coordinate):
                    result = clustering.choose_geometries(molecules, maximum_number_of_seeds=2)

        self.assertEqual(descriptor_calls, [0.0, 1.0, 10.0])
        self.assertEqual([m.name for m in result], ["m0", "m2"])

    def test_hybrid_clusters_then_fills_with_maxmin(self):
        molecules = [
            SimpleNamespace(name="m0", atoms_list=["H"], coordinates=[[0.0, 0.0, 0.0]], energy=0.0),
            SimpleNamespace(name="m1", atoms_list=["H"], coordinates=[[1.0, 0.0, 0.0]], energy=1.0),
            SimpleNamespace(name="m2", atoms_list=["H"], coordinates=[[10.0, 0.0, 0.0]], energy=2.0),
            SimpleNamespace(name="m3", atoms_list=["H"], coordinates=[[20.0, 0.0, 0.0]], energy=3.0),
        ]

        def descriptor_from_coordinate(_atoms, coordinates):
            return [coordinates[0][0]]

        with mock.patch.dict("os.environ", {"PYAR_CLUSTERING_ALGORITHM": "hybrid"}):
            with mock.patch("pyar.representations.mbtr_descriptor", side_effect=descriptor_from_coordinate):
                with mock.patch("pyar.data_analysis.clustering.generate_labels", return_value=[0, 0, 1, 1]):
                    result = clustering.choose_geometries(molecules, maximum_number_of_seeds=3)

        self.assertEqual([m.name for m in result], ["m0", "m2", "m3"])

    def test_hybrid_trims_cluster_minima_with_maxmin_when_too_many_clusters(self):
        molecules = [
            SimpleNamespace(name="m0", atoms_list=["H"], coordinates=[[0.0, 0.0, 0.0]], energy=0.0),
            SimpleNamespace(name="m1", atoms_list=["H"], coordinates=[[1.0, 0.0, 0.0]], energy=1.0),
            SimpleNamespace(name="m2", atoms_list=["H"], coordinates=[[10.0, 0.0, 0.0]], energy=2.0),
            SimpleNamespace(name="m3", atoms_list=["H"], coordinates=[[20.0, 0.0, 0.0]], energy=3.0),
        ]

        def descriptor_from_coordinate(_atoms, coordinates):
            return [coordinates[0][0]]

        with mock.patch.dict("os.environ", {"PYAR_CLUSTERING_ALGORITHM": "hybrid"}):
            with mock.patch("pyar.representations.mbtr_descriptor", side_effect=descriptor_from_coordinate):
                    with mock.patch("pyar.data_analysis.clustering.generate_labels", return_value=[0, 1, 2, 3]):
                        result = clustering.choose_geometries(molecules, maximum_number_of_seeds=2)

        self.assertEqual([m.name for m in result], ["m0", "m3"])

    def test_noise_points_keep_best_representative(self):
        molecules = [
            SimpleNamespace(name="noise_low", atoms_list=["H"], coordinates=[[0.0, 0.0, 0.0]], energy=0.0),
            SimpleNamespace(name="cluster", atoms_list=["H"], coordinates=[[1.0, 0.0, 0.0]], energy=1.0),
            SimpleNamespace(name="noise_high", atoms_list=["H"], coordinates=[[2.0, 0.0, 0.0]], energy=2.0),
        ]

        best = clustering.select_best_from_each_cluster([-1, 0, -1], molecules)

        self.assertEqual([m.name for m in best], ["cluster", "noise_low"])

    def test_generate_labels_supports_optics_and_spectral(self):
        dt = np.array([[0.0], [1.0], [2.0], [10.0]])

        optics_labels = clustering.generate_labels(dt, "optics", maximum_number_of_seeds=2)
        spectral_labels = clustering.generate_labels(dt, "spectral", maximum_number_of_seeds=2)

        self.assertEqual(len(optics_labels), 4)
        self.assertEqual(len(spectral_labels), 4)

    def test_basin_memory_reduces_candidate_pool_before_mbtr(self):
        molecules = [
            SimpleNamespace(name=f"m{i}", atoms_list=["H"], coordinates=[[float(i), 0.0, 0.0]], energy=float(i))
            for i in range(6)
        ]
        basin_entries = [
            {"fingerprint": [0.0, 6.0]},
            {"fingerprint": [1.0, 5.0]},
        ]
        descriptor_calls = []

        def fingerprint_from_coordinate(_atoms, coordinates):
            value = coordinates[0][0]
            return [value, 6.0 - value]

        def descriptor_from_coordinate(_atoms, coordinates):
            descriptor_calls.append(coordinates[0][0])
            return [coordinates[0][0]]

        with mock.patch.dict("os.environ", {"PYAR_CLUSTERING_ALGORITHM": "maxmin"}):
            with mock.patch("pyar.data_analysis.clustering.remove_similar", return_value=molecules):
                with mock.patch("pyar.data_analysis.clustering._load_basin_registry", return_value=basin_entries):
                    with mock.patch("pyar.representations.fingerprint", side_effect=fingerprint_from_coordinate):
                        with mock.patch("pyar.representations.mbtr_descriptor", side_effect=descriptor_from_coordinate):
                            result = clustering.choose_geometries(molecules, maximum_number_of_seeds=2)

        self.assertEqual(set(descriptor_calls), {0.0, 1.0, 4.0, 5.0})
        self.assertEqual(len(descriptor_calls), 4)
        self.assertEqual(len(result), 2)

    def test_basin_registry_persists_unique_entries(self):
        molecules = [
            SimpleNamespace(name="m0", atoms_list=["H"], coordinates=[[0.0, 0.0, 0.0]], energy=0.0),
            SimpleNamespace(name="m1", atoms_list=["H"], coordinates=[[1.0, 0.0, 0.0]], energy=1.0),
        ]

        def fingerprint_from_coordinate(_atoms, coordinates):
            return [coordinates[0][0]]

        with tempfile.TemporaryDirectory() as tmpdir:
            registry_path = os.path.join(tmpdir, "selected", "stoichiometry_H", "basin_registry.json")
            with mock.patch("pyar.representations.fingerprint", side_effect=fingerprint_from_coordinate):
                clustering._persist_basin_registry(registry_path, molecules)
                clustering._persist_basin_registry(registry_path, molecules)

            self.assertTrue(os.path.exists(registry_path))
            with open(registry_path, "r") as fp:
                payload = json.load(fp)

        self.assertEqual(len(payload["entries"]), 2)

    def test_read_energy_from_xyz_file_uses_trailing_energy_token(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            xyz_path = os.path.join(tmpdir, "example.xyz")
            with open(xyz_path, "w") as fp:
                fp.write("1\n")
                fp.write("trial orientation 002_002_ag_a_001_b_004: -40.5243989267341\n")
                fp.write("H 0.0 0.0 0.0\n")

            energy = clustering.read_energy_from_xyz_file(xyz_path)

        self.assertAlmostEqual(energy, -40.5243989267341)


if __name__ == "__main__":
    unittest.main()
