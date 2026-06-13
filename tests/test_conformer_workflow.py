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
        self.assertEqual(summary_rows[0]["source_conf_id"], "2")
        self.assertEqual(summary_rows[0]["selected_path"], result.selected_paths[0])

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
        )

        selected_names = {record.name for record in selected}
        self.assertEqual(selected_names, {"low", "near", "far"})
        self.assertEqual(records[0].refinement_reason, "energy")
        self.assertEqual(records[1].refinement_reason, "energy")
        self.assertEqual(records[3].refinement_reason, "diversity")
        self.assertGreater(records[3].refinement_diversity, records[2].refinement_diversity or 0.0)

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
                    use_random_coords=False,
                    root_directory=tmpdir,
                )

        self.assertEqual(calls, [10, 11])
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
