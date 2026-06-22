import builtins
import csv
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from pyar.core.molecule import Molecule
from pyar.workflows import conformer


class ConformerWorkflowTests(unittest.TestCase):
    def _molecule_factory(self):
        def build_molecule(rdkit_molecule, conformer_id, *, name, seed, charge, multiplicity, scftype, energy):
            return Molecule(
                ["C", "H"],
                [[0.0, 0.0, float(conformer_id)], [1.0, 0.0, float(conformer_id)]],
                name=name,
                title=f"seed {seed} conformer {conformer_id}",
                charge=0 if charge is None else charge,
                multiplicity=multiplicity,
                scftype=scftype,
                energy=energy,
            )

        return build_molecule

    def test_conformer_search_writes_ranked_rdkit_outputs(self):
        records = [
            conformer.ConformerRecord(1, 2, 0.5, "converged", "mmff"),
            conformer.ConformerRecord(1, 1, 1.5, "converged", "mmff"),
            conformer.ConformerRecord(1, 3, 2.5, "not_converged", "mmff"),
        ]
        with tempfile.TemporaryDirectory() as tmpdir:
            with mock.patch.object(conformer, "_rdkit_modules", return_value=(object(), object())), \
                mock.patch.object(conformer, "_load_rdkit_molecule", return_value=(object(), "smiles")), \
                mock.patch.object(conformer, "_embed_and_rank_conformers", return_value=records), \
                mock.patch.object(
                    conformer,
                    "_molecule_from_rdkit_conformer",
                    side_effect=self._molecule_factory(),
                ):
                result = conformer.conformer_search(
                    "CCO",
                    num_conformers=3,
                    top_n=2,
                    num_seeds=1,
                    use_random_coords=True,
                    root_directory=tmpdir,
                )

            run_directory = Path(result.run_directory)
            state = json.loads((run_directory / "state.json").read_text())
            summary_rows = list(csv.DictReader((run_directory / "summary.csv").open()))
            rdkit_file_exists = (run_directory / "rdkit" / "seed_000001_conf_0002.xyz").is_file()
            second_rdkit_file_exists = (run_directory / "rdkit" / "seed_000001_conf_0001.xyz").is_file()
            selected_file_exists = (run_directory / "selected" / "conf_0000.xyz").is_file()

        self.assertEqual(result.status, "completed")
        self.assertEqual(len(result.selected_paths), 2)
        self.assertTrue(rdkit_file_exists)
        self.assertTrue(second_rdkit_file_exists)
        self.assertTrue(selected_file_exists)
        self.assertEqual(state["generated_conformers"], 3)
        self.assertEqual(state["source_format"], "smiles")
        self.assertEqual(state["request"]["torsion_rounds"], 2)
        self.assertEqual(state["request"]["torsion_mode"], "random")
        self.assertEqual(state["request"]["torsion_kicks_per_conformer"], 6)
        self.assertEqual(state["request"]["generation_dedup_rms"], 0.5)
        self.assertEqual(state["request"]["torsion_dedup_rms"], 0.5)
        self.assertEqual(summary_rows[0]["source_conf_id"], "2")
        self.assertEqual(summary_rows[0]["generation_stage"], "etkdg")
        self.assertEqual(summary_rows[0]["selected_path"], result.selected_paths[0])

    def test_conformer_search_dispatches_stratified_random_torsion_kicks(self):
        records = [conformer.ConformerRecord(1, 0, 0.5, "converged", "mmff")]
        with tempfile.TemporaryDirectory() as tmpdir:
            with mock.patch.object(conformer, "_rdkit_modules", return_value=(object(), object())), \
                mock.patch.object(conformer, "_load_rdkit_molecule", return_value=(object(), "smiles")), \
                mock.patch.object(conformer, "_embed_and_rank_conformers", return_value=records), \
                mock.patch.object(conformer, "_generate_torsion_random_records", return_value=[]) as random_kicks, \
                mock.patch.object(
                    conformer,
                    "_molecule_from_rdkit_conformer",
                    side_effect=self._molecule_factory(),
                ):
                conformer.conformer_search(
                    "CCO",
                    num_conformers=1,
                    top_n=1,
                    num_seeds=1,
                    torsion_mode="torsion-kick",
                    torsion_rounds=1,
                    root_directory=tmpdir,
                )

        random_kicks.assert_called_once()
        self.assertEqual(random_kicks.call_args.kwargs["kicks_per_conformer"], 6)
        self.assertEqual(random_kicks.call_args.kwargs["max_bonds"], 3)
        self.assertNotIn("p_xgm_initial", random_kicks.call_args.kwargs)

    def test_conformer_search_rejects_removed_torsion_modes(self):
        for mode in ("adaptive", "basinhop", "bonobo", "evolve", "mc", "grid"):
            with self.subTest(mode=mode):
                with self.assertRaisesRegex(conformer.ConformerWorkflowError, "--torsion-mode must be 'random'"):
                    conformer.conformer_search("CCO", torsion_mode=mode)

    def test_conformer_search_refines_only_top_ranked_records_with_backend(self):
        records = [
            conformer.ConformerRecord(1, 0, 1.0, "converged", "mmff"),
            conformer.ConformerRecord(1, 1, 2.0, "converged", "mmff"),
            conformer.ConformerRecord(1, 2, 3.0, "converged", "mmff"),
        ]
        optimized = []

        def fake_optimise(molecule, qc_params):
            optimized.append(molecule.name)
            job_dir = Path(f"job_{molecule.name}")
            job_dir.mkdir()
            (job_dir / f"result_{molecule.name}.xyz").write_text(
                "2\noptimized\nC 0 0 0\nH 1 0 0\n"
            )
            molecule.energy = 5.0 if molecule.name.endswith("conf_0001") else 10.0
            return True

        with tempfile.TemporaryDirectory() as tmpdir:
            with mock.patch.object(conformer, "_rdkit_modules", return_value=(object(), object())), \
                mock.patch.object(conformer, "_load_rdkit_molecule", return_value=(object(), "smiles")), \
                mock.patch.object(conformer, "_embed_and_rank_conformers", return_value=records), \
                mock.patch.object(
                    conformer,
                    "_molecule_from_rdkit_conformer",
                    side_effect=self._molecule_factory(),
                ), \
                mock.patch.object(conformer, "optimise", side_effect=fake_optimise):
                result = conformer.conformer_search(
                    "CCO",
                    num_conformers=3,
                    top_n=2,
                    backend_top_n=3,
                    num_seeds=1,
                    use_random_coords=True,
                    qc_params={"software": "xtb"},
                    root_directory=tmpdir,
                )

            state = json.loads((Path(result.run_directory) / "state.json").read_text())

        self.assertEqual(
            optimized,
            [
                "seed_000001_conf_0000",
                "seed_000001_conf_0001",
                "seed_000001_conf_0002",
            ],
        )
        self.assertEqual(result.status, "completed")
        self.assertEqual(len(result.selected_paths), 2)
        self.assertEqual(state["records"][1]["backend_energy"], 5.0)
        self.assertTrue(result.selected_paths[0].endswith("selected/conf_0000.xyz"))

    def test_refinement_selection_adds_geometric_diversity(self):
        def make_record(name, source_conf_id, energy, coordinates):
            record = conformer.ConformerRecord(
                1,
                source_conf_id,
                energy,
                "converged",
                "mmff",
            )
            record.name = name
            record.molecule = Molecule(
                ["C", "C", "C"],
                coordinates,
                name=name,
                energy=energy,
            )
            return record

        records = [
            make_record("low", 0, 0.0, [[0, 0, 0], [1, 0, 0], [0, 1, 0]]),
            make_record("near", 1, 1.0, [[0, 0, 0], [1.02, 0, 0], [0, 1, 0]]),
            make_record("moderate", 2, 2.0, [[0, 0, 0], [1, 0, 0], [0, 1.5, 0]]),
            make_record("far", 3, 9.0, [[0, 0, 0], [3, 0, 0], [0, 3, 0]]),
        ]

        selected = conformer._select_refinement_records(
            records,
            limit=3,
            top_n=1,
            diversity_fraction=0.5,
            compactness_fraction=0.5,
        )

        selected_names = {record.name for record in selected}
        self.assertEqual(selected_names, {"low", "near", "moderate"})
        self.assertEqual(records[0].refinement_reason, "energy")
        self.assertEqual(records[1].refinement_reason, "contact_compact")
        self.assertEqual(records[2].refinement_reason, "contact_compact")
        self.assertIsNotNone(records[1].radius_gyration)

    def test_refinement_selection_protects_compact_and_open_basins(self):
        def make_record(name, source_conf_id, energy, coordinates):
            record = conformer.ConformerRecord(1, source_conf_id, energy, "converged", "mmff")
            record.name = name
            record.molecule = Molecule(["C", "C", "C"], coordinates, name=name, energy=energy)
            return record

        records = [
            make_record("energy", 0, 0.0, [[0, 0, 0], [1, 0, 0], [0, 1, 0]]),
            make_record("near_energy", 1, 0.1, [[0, 0, 0], [1.05, 0, 0], [0, 1.05, 0]]),
            make_record("compact", 2, 5.0, [[0, 0, 0], [0.4, 0, 0], [0, 0.4, 0]]),
            make_record("open", 3, 6.0, [[0, 0, 0], [4, 0, 0], [0, 4, 0]]),
            make_record("diverse", 4, 7.0, [[0, 0, 0], [2, 0, 0], [0, 5, 0]]),
        ]

        selected = conformer._select_refinement_records(
            records,
            limit=5,
            top_n=1,
            diversity_fraction=0.2,
            compactness_fraction=0.2,
        )

        reasons = {record.name: record.refinement_reason for record in selected}
        self.assertEqual(reasons["compact"], "contact_compact")
        self.assertEqual(reasons["open"], "open")
        self.assertIn("diversity", reasons.values())

    def test_deduplicate_records_keeps_low_energy_unique_conformers(self):
        def make_record(name, source_conf_id, energy, coordinates):
            record = conformer.ConformerRecord(1, source_conf_id, energy, "converged", "mmff")
            record.name = name
            record.molecule = Molecule(["C", "C", "C"], coordinates, name=name, energy=energy)
            return record

        records = [
            make_record("low", 0, 0.0, [[0, 0, 0], [1, 0, 0], [0, 1, 0]]),
            make_record("duplicate", 1, 1.0, [[0, 0, 0], [1.05, 0, 0], [0, 1.05, 0]]),
            make_record("unique", 2, 2.0, [[0, 0, 0], [3, 0, 0], [0, 3, 0]]),
        ]

        selected = conformer._deduplicate_records(records, 0.5)

        self.assertEqual([record.name for record in selected], ["low", "unique"])

    def test_generation_pool_collapses_without_torsion_kicks(self):
        records = [
            conformer.ConformerRecord(1, 0, 0.0, "converged", "mmff"),
            conformer.ConformerRecord(1, 1, 1.0, "converged", "mmff"),
            conformer.ConformerRecord(1, 2, 2.0, "converged", "mmff"),
        ]
        coordinates = {
            0: [[0, 0, 0], [1, 0, 0], [0, 1, 0]],
            1: [[0, 0, 0], [1.03, 0, 0], [0, 1.03, 0]],
            2: [[0, 0, 0], [3, 0, 0], [0, 3, 0]],
        }
        for record in records:
            record.molecule = Molecule(
                ["C", "C", "C"],
                coordinates[record.source_conf_id],
                name=f"record_{record.source_conf_id}",
                energy=record.rdkit_energy,
            )

        def molecule_factory(rdkit_molecule, conformer_id, *, name, seed, charge, multiplicity, scftype, energy):
            return Molecule(
                ["C", "C", "C"],
                coordinates[int(conformer_id)],
                name=name,
                energy=energy,
            )

        with tempfile.TemporaryDirectory() as tmpdir:
            with mock.patch.object(conformer, "_rdkit_modules", return_value=(object(), object())), \
                mock.patch.object(conformer, "_load_rdkit_molecule", return_value=(object(), "smiles")), \
                mock.patch.object(conformer, "_embed_and_rank_conformers", return_value=records), \
                mock.patch.object(
                    conformer,
                    "_molecule_from_rdkit_conformer",
                    side_effect=molecule_factory,
                ):
                result = conformer.conformer_search(
                    "CCO",
                    num_conformers=3,
                    top_n=3,
                    num_seeds=1,
                    torsion_kicks=False,
                    root_directory=tmpdir,
                )

            summary_rows = list(csv.DictReader((Path(result.run_directory) / "summary.csv").open()))

        self.assertEqual([row["source_conf_id"] for row in summary_rows], ["0", "2"])

    def test_select_unique_records_removes_backend_collapsed_duplicates(self):
        def make_record(name, source_conf_id, energy, coordinates):
            record = conformer.ConformerRecord(1, source_conf_id, energy, "converged", "mmff")
            record.name = name
            record.backend_energy = energy
            record.molecule = Molecule(["C", "C", "C"], coordinates, name=name, energy=energy)
            return record

        records = [
            make_record("best", 0, -3.0, [[0, 0, 0], [1, 0, 0], [0, 1, 0]]),
            make_record("collapsed", 1, -2.9, [[0, 0, 0], [1.03, 0, 0], [0, 1.03, 0]]),
            make_record("distinct", 2, -2.0, [[0, 0, 0], [3, 0, 0], [0, 1, 0]]),
        ]

        selected = conformer._select_unique_records(records, 2, 0.75)

        self.assertEqual([record.name for record in selected], ["best", "distinct"])

    def test_torsion_kick_move_count_mixes_local_and_broad_moves(self):
        counts = [
            conformer._torsion_kick_move_count(num_torsions=5, max_bonds=4, kick_index=index)
            for index in range(8)
        ]

        self.assertEqual(counts, [1, 1, 2, 2, 4, 1, 1, 3])

    def test_stratified_torsion_kick_plan_cycles_bonds_and_angle_bins(self):
        rng = conformer.np.random.default_rng(7)

        torsion_indices, angle_indices = conformer._stratified_torsion_kick_plan(
            4,
            5,
            parent_index=1,
            round_index=1,
            kick_index=2,
            kicks_per_conformer=4,
            max_bonds=3,
            rng=rng,
        )

        self.assertEqual(torsion_indices[0], 2)
        self.assertEqual(len(torsion_indices), 2)
        self.assertEqual(angle_indices[0], 1)
        self.assertEqual(len(angle_indices), len(torsion_indices))

    def test_torsion_random_generates_minimized_rdkit_records(self):
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem
        except ImportError:
            self.skipTest("RDKit is not available")

        molecule = Chem.AddHs(Chem.MolFromSmiles("CCCCCC"))
        embedded = conformer._embed_and_rank_conformers(
            molecule,
            AllChem,
            num_conformers=1,
            prune_rms_threshold=0.5,
            force_field="mmff",
            seed=1,
            num_threads=0,
            use_random_coords=True,
            max_iterations=200,
        )

        kicked = conformer._generate_torsion_random_records(
            embedded,
            Chem,
            AllChem,
            enabled=True,
            kicks_per_conformer=4,
            max_bonds=3,
            seed=1,
            round_index=1,
            max_iterations=200,
        )

        self.assertGreaterEqual(len(kicked), 1)
        self.assertEqual(kicked[0].generation_stage, "torsion_kick")
        self.assertEqual(kicked[0].generation_round, 1)
        self.assertEqual(kicked[0].parent_name, "seed_000001_conf_0000")
        self.assertTrue(kicked[0].torsion_moves)
        self.assertTrue(kicked[0].rdkit_status in {"converged", "not_converged"})

    def test_torsion_random_uses_unique_ids_with_many_kicks(self):
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem
        except ImportError:
            self.skipTest("RDKit is not available")

        molecule = Chem.AddHs(Chem.MolFromSmiles("CCCCCC"))
        embedded = conformer._embed_and_rank_conformers(
            molecule,
            AllChem,
            num_conformers=1,
            prune_rms_threshold=0.5,
            force_field="mmff",
            seed=1,
            num_threads=0,
            use_random_coords=True,
            max_iterations=100,
        )

        kicked = conformer._generate_torsion_random_records(
            embedded,
            Chem,
            AllChem,
            enabled=True,
            kicks_per_conformer=12,
            max_bonds=3,
            seed=1,
            round_index=1,
            max_iterations=100,
        )

        source_ids = [record.source_conf_id for record in kicked]
        self.assertGreater(len(source_ids), 1)
        self.assertEqual(len(source_ids), len(set(source_ids)))

    def test_conformer_search_runs_multiple_torsion_rounds(self):
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem
        except ImportError:
            self.skipTest("RDKit is not available")

        with tempfile.TemporaryDirectory() as tmpdir:
            result = conformer.conformer_search(
                "CCCCCC",
                num_conformers=2,
                top_n=2,
                backend_top_n=4,
                num_seeds=1,
                torsion_rounds=2,
                torsion_kicks_per_conformer=2,
                torsion_max_bonds=1,
                torsion_dedup_rms=0.1,
                force_field="mmff",
                root_directory=tmpdir,
            )

            state = json.loads((Path(result.run_directory) / "state.json").read_text())

        rounds = {
            record["generation_round"]
            for record in state["records"]
            if record["generation_stage"] == "torsion_kick"
        }
        self.assertTrue({1, 2}.issubset(rounds))

    def test_conformer_search_spreads_across_multiple_seeds(self):
        records = {
            10: [conformer.ConformerRecord(10, 0, 1.0, "converged", "mmff")],
            11: [conformer.ConformerRecord(11, 1, 0.5, "converged", "mmff")],
        }
        calls = []

        def fake_embed_and_rank(
            rdkit_molecule,
            AllChem,
            *,
            seed,
            num_conformers,
            prune_rms_threshold,
            force_field,
            num_threads,
            use_random_coords,
            max_iterations,
        ):
            calls.append(seed)
            return records[seed]

        with tempfile.TemporaryDirectory() as tmpdir:
            with mock.patch.object(conformer, "_rdkit_modules", return_value=(object(), object())), \
                mock.patch.object(conformer, "_load_rdkit_molecule", return_value=(object(), "smiles")), \
                mock.patch.object(conformer, "_embed_and_rank_conformers", side_effect=fake_embed_and_rank), \
                mock.patch.object(
                    conformer,
                    "_molecule_from_rdkit_conformer",
                    side_effect=self._molecule_factory(),
                ):
                result = conformer.conformer_search(
                    "CCO",
                    num_conformers=1,
                    top_n=2,
                    num_seeds=2,
                    seed=10,
                    compactness_fraction=0.75,
                    use_random_coords=False,
                    root_directory=tmpdir,
                )

        self.assertEqual(calls, [10, 11])
        self.assertEqual(len(result.selected_paths), 2)

    def test_backend_selection_uses_larger_default_refinement_pool(self):
        records = [
            conformer.ConformerRecord(1, 0, 0.1, "converged", "mmff"),
            conformer.ConformerRecord(1, 1, 0.2, "converged", "mmff"),
            conformer.ConformerRecord(1, 2, 0.3, "converged", "mmff"),
            conformer.ConformerRecord(1, 3, 0.4, "converged", "mmff"),
            conformer.ConformerRecord(1, 4, 0.5, "converged", "mmff"),
            conformer.ConformerRecord(1, 5, 0.6, "converged", "mmff"),
        ]
        selected_names = []

        def fake_optimise(molecule, qc_params):
            selected_names.append(molecule.name)
            job_dir = Path(f"job_{molecule.name}")
            job_dir.mkdir()
            (job_dir / f"result_{molecule.name}.xyz").write_text(
                "2\noptimized\nC 0 0 0\nH 1 0 0\n"
            )
            molecule.energy = 0.1
            return True

        with tempfile.TemporaryDirectory() as tmpdir:
            with mock.patch.object(conformer, "_rdkit_modules", return_value=(object(), object())), \
                mock.patch.object(conformer, "_load_rdkit_molecule", return_value=(object(), "smiles")), \
                mock.patch.object(conformer, "_embed_and_rank_conformers", return_value=records), \
                mock.patch.object(
                    conformer,
                    "_molecule_from_rdkit_conformer",
                    side_effect=self._molecule_factory(),
                ), \
                mock.patch.object(conformer, "optimise", side_effect=fake_optimise):
                result = conformer.conformer_search(
                    "CCO",
                    num_conformers=6,
                    top_n=2,
                    num_seeds=1,
                    qc_params={"software": "xtb"},
                    root_directory=tmpdir,
                )

        self.assertGreater(len(selected_names), 2)
        self.assertEqual(len(result.selected_paths), 2)

    def test_missing_rdkit_uses_conformer_extra_hint(self):
        original_import = builtins.__import__

        def block_rdkit(name, *args, **kwargs):
            if name == "rdkit" or name.startswith("rdkit."):
                raise ImportError("blocked rdkit")
            return original_import(name, *args, **kwargs)

        with mock.patch("builtins.__import__", side_effect=block_rdkit):
            with self.assertRaises(ImportError) as context:
                conformer._rdkit_modules()

        self.assertIn("conformer", str(context.exception))
        self.assertIn("pyar-chem[conformer]", str(context.exception))

    def test_xyz_input_uses_openbabel_conversion_before_rdkit_load(self):
        loaded_molecule = object()

        class FakeChem:
            @staticmethod
            def SDMolSupplier(path, removeHs=False):
                return [loaded_molecule]

            @staticmethod
            def SanitizeMol(molecule):
                return None

            @staticmethod
            def AddHs(molecule, addCoords=True):
                return molecule

        def fake_run_command(command, stdout_path=None, stderr_path=None):
            output_path = Path(command[-1])
            output_path.write_text("sdf")
            return 0

        with tempfile.TemporaryDirectory() as tmpdir:
            xyz_path = Path(tmpdir) / "input.xyz"
            xyz_path.write_text("1\nh\nH 0 0 0\n")
            with mock.patch.object(conformer, "require_executable", return_value="/usr/bin/obabel"), \
                mock.patch.object(conformer, "run_command", side_effect=fake_run_command) as run_command:
                molecule, source_format = conformer._load_rdkit_molecule(
                    str(xyz_path),
                    "xyz",
                    FakeChem,
                    Path(tmpdir) / "conformers",
                )

        self.assertIs(molecule, loaded_molecule)
        self.assertEqual(source_format, "xyz")
        command = run_command.call_args.args[0]
        self.assertEqual(command[:3], ["/usr/bin/obabel", "-ixyz", xyz_path])
        self.assertIn("-osdf", command)
        self.assertIn("-O", command)

    def test_auto_force_field_prefers_mmff_and_falls_back_to_uff(self):
        class Params:
            pass

        class FakeAllChem:
            mmff = True

            @staticmethod
            def ETKDGv3():
                return Params()

            @classmethod
            def MMFFHasAllMoleculeParams(cls, molecule):
                return cls.mmff

            @staticmethod
            def UFFHasAllMoleculeParams(molecule):
                return True

            @staticmethod
            def EmbedMultipleConfs(molecule, numConfs, params):
                return [7]

            @staticmethod
            def MMFFOptimizeMoleculeConfs(molecule, numThreads, maxIters):
                return [(0, -1.0)]

            @staticmethod
            def UFFOptimizeMoleculeConfs(molecule, numThreads, maxIters):
                return [(0, -0.5)]

        records = conformer._embed_and_rank_conformers(
            object(),
            FakeAllChem,
            num_conformers=1,
            prune_rms_threshold=0.5,
            force_field="auto",
            seed=1,
            num_threads=0,
            use_random_coords=False,
            max_iterations=200,
        )
        FakeAllChem.mmff = False
        fallback_records = conformer._embed_and_rank_conformers(
            object(),
            FakeAllChem,
            num_conformers=1,
            prune_rms_threshold=0.5,
            force_field="auto",
            seed=1,
            num_threads=0,
            use_random_coords=False,
            max_iterations=200,
        )

        self.assertEqual(records[0].force_field, "mmff")
        self.assertEqual(fallback_records[0].force_field, "uff")

    def test_embed_parameters_sets_random_coords_and_pruning_threshold(self):
        class Params:
            pass

        class FakeAllChem:
            @staticmethod
            def ETKDGv3():
                return Params()

        params = conformer._embed_parameters(
            FakeAllChem,
            seed=7,
            prune_rms_threshold=0.25,
            num_threads=4,
            use_random_coords=True,
        )

        self.assertEqual(params.randomSeed, 7)
        self.assertEqual(params.pruneRmsThresh, 0.25)
        self.assertEqual(params.numThreads, 4)
        self.assertTrue(params.useRandomCoords)


if __name__ == "__main__":
    unittest.main()
