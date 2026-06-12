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
        def build_molecule(rdkit_molecule, conformer_id, *, name, charge, multiplicity, scftype, energy):
            return Molecule(
                ["C", "H"],
                [[0.0, 0.0, float(conformer_id)], [1.0, 0.0, float(conformer_id)]],
                name=name,
                title=f"conformer {conformer_id}",
                charge=0 if charge is None else charge,
                multiplicity=multiplicity,
                scftype=scftype,
                energy=energy,
            )

        return build_molecule

    def test_conformer_search_writes_ranked_rdkit_outputs(self):
        records = [
            conformer.ConformerRecord(2, 0.5, "converged", "mmff"),
            conformer.ConformerRecord(1, 1.5, "converged", "mmff"),
            conformer.ConformerRecord(3, 2.5, "not_converged", "mmff"),
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
                    root_directory=tmpdir,
                )

            run_directory = Path(result.run_directory)
            state = json.loads((run_directory / "state.json").read_text())
            summary_rows = list(csv.DictReader((run_directory / "summary.csv").open()))
            rdkit_file_exists = (run_directory / "rdkit" / "conf_0000.xyz").is_file()
            second_rdkit_file_exists = (run_directory / "rdkit" / "conf_0001.xyz").is_file()
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
            conformer.ConformerRecord(0, 1.0, "converged", "mmff"),
            conformer.ConformerRecord(1, 2.0, "converged", "mmff"),
            conformer.ConformerRecord(2, 3.0, "converged", "mmff"),
        ]
        optimized = []

        def fake_optimise(molecule, qc_params):
            optimized.append(molecule.name)
            job_dir = Path(f"job_{molecule.name}")
            job_dir.mkdir()
            (job_dir / f"result_{molecule.name}.xyz").write_text(
                "2\noptimized\nC 0 0 0\nH 1 0 0\n"
            )
            molecule.energy = 5.0 if molecule.name == "conf_0001" else 10.0
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
                    qc_params={"software": "xtb"},
                    root_directory=tmpdir,
                )

            state = json.loads((Path(result.run_directory) / "state.json").read_text())

        self.assertEqual(optimized, ["conf_0000", "conf_0001"])
        self.assertEqual(result.status, "completed")
        self.assertEqual(len(result.selected_paths), 2)
        self.assertEqual(state["records"][1]["backend_energy"], 5.0)
        self.assertTrue(result.selected_paths[0].endswith("selected/conf_0000.xyz"))

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
            rms_threshold=0.5,
            force_field="auto",
            seed=1,
            num_threads=0,
            max_iterations=200,
        )
        FakeAllChem.mmff = False
        fallback_records = conformer._embed_and_rank_conformers(
            object(),
            FakeAllChem,
            num_conformers=1,
            rms_threshold=0.5,
            force_field="auto",
            seed=1,
            num_threads=0,
            max_iterations=200,
        )

        self.assertEqual(records[0].force_field, "mmff")
        self.assertEqual(fallback_records[0].force_field, "uff")


if __name__ == "__main__":
    unittest.main()
