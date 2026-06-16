import csv
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from pyar.benchmarks import conformer as benchmark
from pyar.workflow_results import ConformerResult


def write_xyz(path, atoms, coordinates, energy=0.0):
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [str(len(atoms)), f"energy: {energy}"]
    lines.extend(
        f"{atom} {coord[0]:.6f} {coord[1]:.6f} {coord[2]:.6f}"
        for atom, coord in zip(atoms, coordinates)
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_benchmark(path, reference_path):
    path.write_text(
        json.dumps(
            {
                "name": "one-case",
                "method": "reference",
                "cases": [
                    {
                        "id": "case1",
                        "input": "CCO",
                        "input_format": "smiles",
                        "charge": 0,
                        "multiplicity": 1,
                        "reference_xyz": str(reference_path),
                        "reference_energy": None,
                        "notes": "synthetic",
                    }
                ],
            }
        ),
        encoding="utf-8",
    )


class ConformerBenchmarkTests(unittest.TestCase):
    def make_case(self, reference_path=None):
        return benchmark.ConformerBenchmarkCase(
            id="case1",
            input="CCO",
            input_format="smiles",
            charge=0,
            multiplicity=1,
            reference_xyz=None if reference_path is None else str(reference_path),
            notes="synthetic",
        )

    def make_result(self, run_directory, *, backend_requested=False, generated=None, selected=None):
        run_directory = Path(run_directory)
        (run_directory / "rdkit").mkdir(parents=True)
        (run_directory / "selected").mkdir(parents=True)
        state = {
            "status": "completed",
            "backend_requested": backend_requested,
        }
        (run_directory / "state.json").write_text(json.dumps(state), encoding="utf-8")
        rows = []
        generated = generated or []
        selected = selected or []
        selected_by_name = {}
        selected_paths = []
        for item in selected:
            selected_name, selected_atoms, selected_coordinates, selected_energy = item
            selected_file = run_directory / "selected" / f"{selected_name}.xyz"
            write_xyz(selected_file, selected_atoms, selected_coordinates, energy=selected_energy)
            selected_by_name[selected_name] = str(selected_file)
            selected_paths.append(str(selected_file))
        for index, item in enumerate(generated):
            name, atoms, coordinates, energy = item[:4]
            backend_energy = item[4] if len(item) > 4 else ""
            refinement_reason = item[5] if len(item) > 5 else ""
            path = run_directory / "rdkit" / f"{name}.xyz"
            write_xyz(path, atoms, coordinates, energy=energy)
            selected_path = selected_by_name.get(name, "")
            rows.append(
                {
                    "rank": index,
                    "name": name,
                    "rdkit_energy": energy,
                    "rdkit_path": str(path),
                    "backend_energy": backend_energy,
                    "selected_path": selected_path,
                    "refinement_reason": refinement_reason,
                    "generation_stage": "test",
                    "backend_status": "completed" if backend_energy != "" else "",
                }
            )
        with (run_directory / "summary.csv").open("w", newline="", encoding="utf-8") as fp:
            writer = csv.DictWriter(
                fp,
                fieldnames=[
                    "rank",
                    "name",
                    "rdkit_energy",
                    "rdkit_path",
                    "backend_energy",
                    "selected_path",
                    "refinement_reason",
                    "generation_stage",
                    "backend_status",
                ],
            )
            writer.writeheader()
            writer.writerows(rows)
        return ConformerResult(
            workflow="conformer",
            status="completed",
            run_directory=str(run_directory),
            state_path=str(run_directory / "state.json"),
            selected_paths=tuple(selected_paths),
        )

    def test_load_json_benchmark(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            spec_path = Path(tmpdir) / "benchmark.json"
            reference_path = Path(tmpdir) / "ref.xyz"
            write_benchmark(spec_path, reference_path)

            spec = benchmark.load_conformer_benchmark(spec_path)

        self.assertEqual(spec.name, "one-case")
        self.assertEqual(len(spec.cases), 1)
        self.assertEqual(spec.cases[0].id, "case1")
        self.assertEqual(spec.cases[0].input, "CCO")

    def test_load_json_requires_case_fields(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            spec_path = Path(tmpdir) / "benchmark.json"
            spec_path.write_text(
                json.dumps({"name": "bad", "method": "reference", "cases": [{"id": "missing-input"}]}),
                encoding="utf-8",
            )

            with self.assertRaisesRegex(benchmark.ConformerBenchmarkError, "input"):
                benchmark.load_conformer_benchmark(spec_path)

    def test_diagnosis_found(self):
        atoms = ["C", "C", "O"]
        reference_coordinates = [[0, 0, 0], [1, 0, 0], [0, 1, 0]]
        with tempfile.TemporaryDirectory() as tmpdir:
            reference = Path(tmpdir) / "ref.xyz"
            write_xyz(reference, atoms, reference_coordinates)
            result = self.make_result(
                Path(tmpdir) / "conformers",
                generated=[("conf_0000", atoms, reference_coordinates, 0.0)],
                selected=[("conf_0000", atoms, reference_coordinates, 0.0)],
            )

            diagnosis = benchmark.diagnose_conformer_case(self.make_case(reference), result)

        self.assertEqual(diagnosis["failure_class"], "found")
        self.assertTrue(diagnosis["reference_found_generated"])
        self.assertTrue(diagnosis["reference_found_selected"])
        self.assertEqual(diagnosis["reference_selection_reason"], "selected")
        self.assertEqual(diagnosis["reference_rdkit_rank"], 0)
        self.assertEqual(diagnosis["lowest_rdkit_energy_native"], 0.0)

    def test_diagnosis_generated_lost_selection(self):
        atoms = ["C", "C", "O"]
        reference_coordinates = [[0, 0, 0], [1, 0, 0], [0, 1, 0]]
        far_coordinates = [[0, 0, 0], [3, 0, 0], [0, 3, 0]]
        with tempfile.TemporaryDirectory() as tmpdir:
            reference = Path(tmpdir) / "ref.xyz"
            write_xyz(reference, atoms, reference_coordinates)
            result = self.make_result(
                Path(tmpdir) / "conformers",
                generated=[("gen0", atoms, reference_coordinates, 0.0)],
                selected=[("conf_0000", atoms, far_coordinates, 1.0)],
            )

            diagnosis = benchmark.diagnose_conformer_case(self.make_case(reference), result)

        self.assertEqual(diagnosis["failure_class"], "generated_lost_selection")
        self.assertEqual(diagnosis["reference_selection_reason"], "not_in_selected_outputs")

    def test_diagnosis_keeps_rdkit_energy_window_in_native_units(self):
        atoms = ["C", "C", "O"]
        reference_coordinates = [[0, 0, 0], [1, 0, 0], [0, 1, 0]]
        far_coordinates = [[0, 0, 0], [3, 0, 0], [0, 3, 0]]
        with tempfile.TemporaryDirectory() as tmpdir:
            reference = Path(tmpdir) / "ref.xyz"
            write_xyz(reference, atoms, reference_coordinates)
            result = self.make_result(
                Path(tmpdir) / "conformers",
                generated=[
                    ("low_energy_far", atoms, far_coordinates, 0.0),
                    ("reference_like", atoms, reference_coordinates, 0.01),
                ],
                selected=[("conf_0000", atoms, far_coordinates, 0.0)],
            )

            diagnosis = benchmark.diagnose_conformer_case(
                self.make_case(reference),
                result,
                energy_window=5.0,
            )

        self.assertEqual(diagnosis["failure_class"], "generated_lost_selection")
        self.assertEqual(diagnosis["reference_rdkit_energy_native"], 0.01)
        self.assertEqual(diagnosis["reference_rdkit_rank"], 1)

    def test_diagnosis_wrong_ranking_uses_rdkit_native_energy_window(self):
        atoms = ["C", "C", "O"]
        reference_coordinates = [[0, 0, 0], [1, 0, 0], [0, 1, 0]]
        far_coordinates = [[0, 0, 0], [3, 0, 0], [0, 3, 0]]
        with tempfile.TemporaryDirectory() as tmpdir:
            reference = Path(tmpdir) / "ref.xyz"
            write_xyz(reference, atoms, reference_coordinates)
            result = self.make_result(
                Path(tmpdir) / "conformers",
                generated=[
                    ("low_energy_far", atoms, far_coordinates, 0.0),
                    ("reference_like", atoms, reference_coordinates, 10.0),
                ],
                selected=[("conf_0000", atoms, far_coordinates, 0.0)],
            )

            diagnosis = benchmark.diagnose_conformer_case(
                self.make_case(reference),
                result,
                energy_window=5.0,
            )

        self.assertEqual(diagnosis["failure_class"], "wrong_ranking")

    def test_diagnosis_backend_loss_records_reference_metadata(self):
        atoms = ["C", "C", "O"]
        reference_coordinates = [[0, 0, 0], [1, 0, 0], [0, 1, 0]]
        far_coordinates = [[0, 0, 0], [3, 0, 0], [0, 3, 0]]
        with tempfile.TemporaryDirectory() as tmpdir:
            reference = Path(tmpdir) / "ref.xyz"
            write_xyz(reference, atoms, reference_coordinates)
            result = self.make_result(
                Path(tmpdir) / "conformers",
                backend_requested=True,
                generated=[
                    ("reference_like", atoms, reference_coordinates, 0.0, -1.0, "energy"),
                    ("low_backend_far", atoms, far_coordinates, 1.0, -1.1, "diversity"),
                ],
                selected=[("conf_0000", atoms, far_coordinates, 1.0)],
            )

            diagnosis = benchmark.diagnose_conformer_case(self.make_case(reference), result)

        self.assertEqual(diagnosis["failure_class"], "generated_lost_backend_or_dedup")
        self.assertTrue(diagnosis["reference_lost_after_dedup"])
        self.assertEqual(diagnosis["reference_refinement_reason"], "energy")
        self.assertEqual(diagnosis["reference_backend_rank"], 1)

    def test_diagnosis_never_generated(self):
        atoms = ["C", "C", "O"]
        reference_coordinates = [[0, 0, 0], [1, 0, 0], [0, 1, 0]]
        far_coordinates = [[0, 0, 0], [3, 0, 0], [0, 3, 0]]
        with tempfile.TemporaryDirectory() as tmpdir:
            reference = Path(tmpdir) / "ref.xyz"
            write_xyz(reference, atoms, reference_coordinates)
            result = self.make_result(
                Path(tmpdir) / "conformers",
                generated=[("gen0", atoms, far_coordinates, 0.0)],
                selected=[("conf_0000", atoms, far_coordinates, 0.0)],
            )

            diagnosis = benchmark.diagnose_conformer_case(self.make_case(reference), result)

        self.assertEqual(diagnosis["failure_class"], "never_generated")

    def test_diagnosis_input_chemistry_issue_for_missing_reference(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            result = self.make_result(Path(tmpdir) / "conformers")

            diagnosis = benchmark.diagnose_conformer_case(self.make_case(Path(tmpdir) / "missing.xyz"), result)

        self.assertEqual(diagnosis["failure_class"], "input_chemistry_issue")
        self.assertIn("missing", diagnosis["reason"])

    def test_summary_and_case_outputs_are_written_with_mocked_workflow(self):
        atoms = ["C", "C", "O"]
        reference_coordinates = [[0, 0, 0], [1, 0, 0], [0, 1, 0]]
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            reference = tmpdir / "ref.xyz"
            write_xyz(reference, atoms, reference_coordinates)
            spec_path = tmpdir / "benchmark.json"
            write_benchmark(spec_path, reference)
            output = tmpdir / "results"

            def fake_search(*args, **kwargs):
                return self.make_result(
                    Path(kwargs["root_directory"]) / "conformers",
                    generated=[("gen0", atoms, reference_coordinates, 0.0)],
                    selected=[("conf_0000", atoms, reference_coordinates, 0.0)],
                )

            with mock.patch.object(benchmark.conformer_workflow, "conformer_search", side_effect=fake_search):
                result = benchmark.run_conformer_benchmark(
                    spec_path,
                    output=output,
                    num_conformers=5,
                    num_seeds=1,
                    top_n=2,
                )

            summary_csv = output / "benchmark_summary.csv"
            summary_json = output / "benchmark_summary.json"
            diagnosis_json = output / "cases" / "case1" / "diagnosis.json"
            self.assertEqual(result["cases"][0]["failure_class"], "found")
            self.assertTrue(summary_csv.is_file())
            self.assertTrue(summary_json.is_file())
            self.assertTrue(diagnosis_json.is_file())
            with summary_csv.open(encoding="utf-8", newline="") as fp:
                row = next(csv.DictReader(fp))
            self.assertIn("lowest_rdkit_energy_native", row)
            self.assertIn("reference_selection_reason", row)

    def test_cli_help_smoke(self):
        from pyar.scripts import conformer_benchmark

        with self.assertRaises(SystemExit) as context:
            conformer_benchmark.argument_parse(["--help"])

        self.assertEqual(context.exception.code, 0)

    def test_cli_invokes_runner(self):
        from pyar.scripts import conformer_benchmark

        with tempfile.TemporaryDirectory() as tmpdir:
            spec_path = Path(tmpdir) / "benchmark.json"
            write_benchmark(spec_path, Path(tmpdir) / "ref.xyz")
            with mock.patch.object(
                conformer_benchmark,
                "run_conformer_benchmark",
                return_value={"name": "one-case", "cases": []},
            ) as runner:
                result = conformer_benchmark.main([str(spec_path), "--output", str(Path(tmpdir) / "out")])

        self.assertIsNone(result)
        self.assertEqual(runner.call_args.kwargs["num_conformers"], 100)
        self.assertIsNone(runner.call_args.kwargs["qc_params"])

    def test_rdkit_smoke_when_available(self):
        try:
            import rdkit  # noqa: F401
        except ImportError:
            self.skipTest("RDKit is not available")

        atoms = ["C", "C", "C", "C"]
        reference_coordinates = [[0, 0, 0], [1.5, 0, 0], [3.0, 0, 0], [4.5, 0, 0]]
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            reference = tmpdir / "butane_reference.xyz"
            write_xyz(reference, atoms, reference_coordinates)
            spec_path = tmpdir / "benchmark.json"
            spec_path.write_text(
                json.dumps(
                    {
                        "name": "rdkit-smoke",
                        "method": "reference",
                        "cases": [
                            {
                                "id": "butane",
                                "input": "CCCC",
                                "input_format": "smiles",
                                "reference_xyz": str(reference),
                            }
                        ],
                    }
                ),
                encoding="utf-8",
            )

            result = benchmark.run_conformer_benchmark(
                spec_path,
                output=tmpdir / "results",
                num_conformers=5,
                num_seeds=1,
                top_n=2,
                conformer_options={"torsion_kicks": False},
            )

        self.assertEqual(result["name"], "rdkit-smoke")
        self.assertEqual(len(result["cases"]), 1)
