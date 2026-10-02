import csv
import json
import shutil
from pathlib import Path

import pytest

from pyar.benchmarks import ts_optimizer as benchmark
from pyar.scripts.ts_optimizer_benchmark import argument_parse


DATA = Path(__file__).parent / "data" / "neb"


def _manifest(tmp_path, *, cases=None, settings=None, source=None):
    raw = {
        "name": "fixture-ts-benchmark",
        "source": source or {
            "name": "Fixture dataset", "version": "1", "license": "CC0-1.0",
            "doi": "10.0000/fixture",
        },
        "qc_model": {"software": "xtb", "xtb_model": "gfn2", "nprocs": 1},
        "settings": settings or {"ts_fmax": 0.03, "ts_max_cycles": 40},
        "cases": cases or [{
            "id": "case_0001_medium", "difficulty": "medium",
            "ts_guess": str(DATA / "guess.xyz"),
            "reference_ts": str(DATA / "guess.xyz"),
            "reactant": str(DATA / "hcn.xyz"), "product": str(DATA / "hnc.xyz"),
            "charge": 1, "multiplicity": 2, "seed": 17,
            "perturbation": {"type": "deterministic_displacement", "seed": 17},
            "reference_energy_hartree": -19.0,
        }],
    }
    path = tmp_path / "benchmark.json"
    path.write_text(json.dumps(raw), encoding="utf-8")
    return path


def _fake_run_neb(monkeypatch, *, fail_optimizer=None, nonstationary=False):
    calls = []

    def fake_run_neb(**kwargs):
        output = Path(kwargs["output"])
        stage = kwargs["stage"]
        optimizer = kwargs.get("ts_optimizer", "current")
        if stage == "ts":
            optimizer = kwargs["ts_optimizer"]
        calls.append((stage, dict(kwargs)))
        if stage == "ts" and optimizer == fail_optimizer:
            raise RuntimeError("fixture optimizer failure")

        result = {}
        if stage == "ts":
            shutil.copy2(kwargs["ts_geometry"], output / "ts_optimized.xyz")
            result = {
                "ts_optimization_converged": True,
                "ts_energy_hartree": -20.0,
                "optimizer_steps": 8,
                "backend_energy_gradient_evaluations": 9 if optimizer == "sella" else 7,
                "ts_optimization_wall_seconds": 10.0 if optimizer == "sella" else 8.0,
                "wall_seconds": 11.0,
            }
        elif stage == "frequency":
            result = {
                "stationary": not nonstationary,
                "first_order_saddle_confirmed": not nonstationary,
                "imaginary_frequency_count": 1 if not nonstationary else 0,
                "backend_energy_gradient_evaluations": 3,
                "wall_seconds": 5.0,
                "hessian_source": "finite_difference_cartesian",
                "hessian_evaluations": 11,
                "hessian_wall_seconds": 3.5,
            }
        elif stage == "relax":
            result = {"backend_energy_gradient_evaluations": 2, "wall_seconds": 7.0}
        elif stage == "irc":
            result = {"irc_converged": True, "backend_energy_gradient_evaluations": 4,
                      "wall_seconds": 8.0}
        elif stage == "endpoints":
            connected = "sella" not in str(output)
            result = {"reactant_product_connection_confirmed": connected,
                      "backend_energy_gradient_evaluations": 5, "wall_seconds": 9.0,
                      "forward_frequency": {
                          "hessian_source": "finite_difference_cartesian",
                          "hessian_evaluations": 13, "hessian_wall_seconds": 4.0,
                      },
                      "backward_frequency": {
                          "hessian_source": "finite_difference_cartesian",
                          "hessian_evaluations": 17, "hessian_wall_seconds": 6.0,
                      }}
        (output / f"{stage}_summary.json").write_text(
            json.dumps({"stage": stage, "status": "complete", **result}), encoding="utf-8",
        )
        return {"stage": stage, **result}

    monkeypatch.setattr(benchmark, "run_neb", fake_run_neb)
    monkeypatch.setattr(benchmark, "_backend_version", lambda qc: {
        "executable": "xtb", "path": "/fixture/xtb", "version": "fixture-version",
    })
    monkeypatch.setattr(benchmark, "_git_revision", lambda: "fixture-commit")
    return calls


def test_manifest_resolves_paths_and_requires_explicit_common_protocol(tmp_path):
    manifest = _manifest(tmp_path)
    spec = benchmark.load_ts_optimizer_benchmark(manifest)
    assert spec.cases[0].id == "case_0001_medium"
    assert Path(spec.cases[0].ts_guess).is_absolute()
    assert spec.qc_model["xtb_model"] == "gfn2"
    assert spec.settings["ts_fmax"] == 0.03
    assert spec.settings["ts_max_cycles"] == 40

    with pytest.raises(benchmark.TSOptimizerBenchmarkError, match="explicitly define xtb_model"):
        _manifest(tmp_path, cases=[{
            "id": "a", "ts_guess": str(DATA / "guess.xyz"),
            "reference_ts": str(DATA / "guess.xyz"),
        }], settings={"ts_fmax": 0.03, "ts_max_cycles": 40})
        raw = json.loads(manifest.read_text())
        raw["qc_model"].pop("xtb_model")
        manifest.write_text(json.dumps(raw))
        benchmark.load_ts_optimizer_benchmark(manifest)


@pytest.mark.parametrize("mutate, message", [
    (lambda x: x["cases"].append(dict(x["cases"][0])), "duplicate case id"),
    (lambda x: x["cases"][0].update(ts_guess="missing.xyz"), "does not exist"),
    (lambda x: x["settings"].update(ts_fmax=0), "ts_fmax must be positive"),
    (lambda x: x["settings"].update(ts_max_cycles=0), "ts_max_cycles must be a positive integer"),
    (lambda x: x["cases"][0].update(multiplicity=0), "multiplicity must be a positive integer"),
    (lambda x: x["cases"][0].update(product=None), "both reactant and product"),
    (lambda x: x["settings"].update(silent_fallback=True), "unsupported settings key"),
])
def test_manifest_rejects_invalid_inputs(tmp_path, mutate, message):
    path = _manifest(tmp_path)
    raw = json.loads(path.read_text())
    mutate(raw)
    path.write_text(json.dumps(raw))
    with pytest.raises(benchmark.TSOptimizerBenchmarkError, match=message):
        benchmark.load_ts_optimizer_benchmark(path)


def test_paired_runs_keep_same_geometry_and_separate_validation_costs(
    tmp_path, monkeypatch,
):
    calls = _fake_run_neb(monkeypatch)
    manifest = _manifest(tmp_path)
    geometric = benchmark.run_ts_optimizer_case(
        manifest, case_id="case_0001_medium", optimizer="geometric", output=tmp_path / "run",
    )
    sella = benchmark.run_ts_optimizer_case(
        manifest, case_id="case_0001_medium", optimizer="sella", output=tmp_path / "run",
    )

    assert [stage for stage, _ in calls] == [
        "ts", "frequency", "relax", "irc", "endpoints",
        "ts", "frequency", "relax", "irc", "endpoints",
    ]
    for stage, kwargs in (calls[0], calls[5]):
        assert kwargs["software"] == "xtb"
        assert kwargs["xtb_model"] == "gfn2"
        assert kwargs["charge"] == 1 and kwargs["multiplicity"] == 2
        assert kwargs["ts_fmax"] == 0.03 and kwargs["ts_max_cycles"] == 40
    assert geometric["outcome"] == "reaction_connected_success"
    assert sella["outcome"] == "first_order_saddle_wrong_connection"
    assert geometric["status"] == sella["status"] == "complete"
    assert geometric["reference_rmsd_angstrom"] == pytest.approx(0.0)
    assert geometric["reference_energy_difference_hartree"] == pytest.approx(-1.0)
    assert geometric["ts_backend_energy_gradient_evaluations"] == 7
    assert geometric["validation_backend_energy_gradient_evaluations"] == 14
    assert geometric["total_backend_energy_gradient_evaluations"] == 21
    assert geometric["validation_wall_seconds"] == 29.0
    assert geometric["validation_hessian_evaluations"] == 41
    assert geometric["validation_hessian_wall_seconds"] == 13.5
    assert geometric["hessian_sources"] == ["finite_difference_cartesian"]
    for optimizer in ("geometric", "sella"):
        directory = tmp_path / "run" / "cases" / "case_0001_medium" / optimizer
        assert (directory / "input_ts_guess.xyz").is_file()
        assert (directory / "input_reference_ts.xyz").is_file()
        assert (directory / "ts_summary.json").is_file()

    summary = benchmark.collect_ts_optimizer_benchmark(tmp_path / "run")
    assert summary["paired_cases"] == 1
    assert (tmp_path / "run" / "benchmark_manifest.json").read_bytes() == manifest.read_bytes()
    assert summary["optimizers"]["geometric"]["success_fraction_of_completed"] == 1.0
    assert summary["optimizers"]["sella"]["success_fraction_of_completed"] == 0.0
    with (tmp_path / "run" / "paired.csv").open(newline="") as stream:
        paired = next(csv.DictReader(stream))
    assert paired["pair_status"] == "paired"
    assert paired["discordant_success"] == "True"
    assert paired["total_evaluation_difference_geometric_minus_sella"] == "-2"


def test_collector_reports_missing_pair_without_counting_it_as_failure(tmp_path, monkeypatch):
    _fake_run_neb(monkeypatch)
    manifest = _manifest(tmp_path)
    benchmark.run_ts_optimizer_case(
        manifest, case_id="case_0001_medium", optimizer="geometric", output=tmp_path / "run",
    )
    summary = benchmark.collect_ts_optimizer_benchmark(tmp_path / "run")
    assert summary["expected_runs"] == 2
    assert summary["incomplete_runs"] == 1
    assert summary["optimizers"]["sella"]["completed_runs"] == 0
    assert summary["optimizers"]["sella"]["outcomes"] == {"incomplete_run": 1}
    with (tmp_path / "run" / "runs.csv").open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    assert [row["outcome"] for row in rows] == ["reaction_connected_success", "incomplete_run"]


def test_collector_refuses_pairs_with_different_physical_or_input_identity(tmp_path, monkeypatch):
    _fake_run_neb(monkeypatch)
    manifest = _manifest(tmp_path)
    for optimizer in benchmark.OPTIMIZERS:
        benchmark.run_ts_optimizer_case(
            manifest, case_id="case_0001_medium", optimizer=optimizer, output=tmp_path / "run",
        )
    result_path = tmp_path / "run/cases/case_0001_medium/sella/benchmark_result.json"
    result = json.loads(result_path.read_text())
    result["qc_model"]["xtb_model"] = "gxtb"
    result_path.write_text(json.dumps(result))
    summary = benchmark.collect_ts_optimizer_benchmark(tmp_path / "run")
    pair = summary["paired"][0]
    assert pair["pair_status"] == "incomplete_or_mismatched"
    assert pair["pair_identity_matches"] is False
    assert pair["total_evaluation_difference_geometric_minus_sella"] is None


def test_collector_refuses_pair_from_different_execution_environment(tmp_path, monkeypatch):
    _fake_run_neb(monkeypatch)
    manifest = _manifest(tmp_path)
    for optimizer in benchmark.OPTIMIZERS:
        benchmark.run_ts_optimizer_case(
            manifest, case_id="case_0001_medium", optimizer=optimizer, output=tmp_path / "run",
        )
    result_path = tmp_path / "run/cases/case_0001_medium/sella/benchmark_result.json"
    result = json.loads(result_path.read_text())
    result["runtime"]["platform"] = "different-worker"
    result_path.write_text(json.dumps(result))
    summary = benchmark.collect_ts_optimizer_benchmark(tmp_path / "run")
    assert summary["paired"][0]["pair_identity_matches"] is False
    assert summary["paired"][0]["pair_status"] == "incomplete_or_mismatched"


def test_collector_refuses_result_stored_under_wrong_optimizer(tmp_path, monkeypatch):
    _fake_run_neb(monkeypatch)
    manifest = _manifest(tmp_path)
    for optimizer in benchmark.OPTIMIZERS:
        benchmark.run_ts_optimizer_case(
            manifest, case_id="case_0001_medium", optimizer=optimizer, output=tmp_path / "run",
        )
    result_path = tmp_path / "run/cases/case_0001_medium/sella/benchmark_result.json"
    result = json.loads(result_path.read_text())
    result["optimizer"] = "geometric"
    result_path.write_text(json.dumps(result))
    summary = benchmark.collect_ts_optimizer_benchmark(tmp_path / "run")
    assert summary["paired"][0]["pair_identity_matches"] is False
    assert summary["paired"][0]["pair_status"] == "incomplete_or_mismatched"


def test_no_endpoint_data_stops_after_independent_frequency_validation(tmp_path, monkeypatch):
    calls = _fake_run_neb(monkeypatch)
    manifest = _manifest(tmp_path)
    raw = json.loads(manifest.read_text())
    raw["cases"][0].pop("reactant")
    raw["cases"][0].pop("product")
    manifest.write_text(json.dumps(raw))
    result = benchmark.run_ts_optimizer_case(
        manifest, case_id="case_0001_medium", optimizer="geometric", output=tmp_path / "run",
    )
    assert result["outcome"] == "validated_first_order_saddle"
    assert [stage for stage, _ in calls] == ["ts", "frequency"]


def test_non_xtb_backend_passes_method_basis_without_xtb_model(tmp_path, monkeypatch):
    calls = _fake_run_neb(monkeypatch)
    manifest = _manifest(tmp_path)
    raw = json.loads(manifest.read_text())
    raw["qc_model"] = {"software": "orca", "method": "PBE0", "basis": "def2-SVP"}
    manifest.write_text(json.dumps(raw))
    result = benchmark.run_ts_optimizer_case(
        manifest, case_id="case_0001_medium", optimizer="geometric", output=tmp_path / "run",
    )
    ts_call = calls[0][1]
    assert ts_call["software"] == "orca"
    assert ts_call["method"] == "PBE0"
    assert ts_call["basis"] == "def2-SVP"
    assert "xtb_model" not in ts_call
    assert "xtb_model" not in result["qc_model"]


def test_run_rejects_a_pair_already_claimed_by_another_job(tmp_path, monkeypatch):
    _fake_run_neb(monkeypatch)
    manifest = _manifest(tmp_path)
    job_dir = tmp_path / "run/cases/case_0001_medium/geometric"
    job_dir.mkdir(parents=True)
    (job_dir / ".benchmark-job.lock").write_text("claimed")
    with pytest.raises(benchmark.TSOptimizerBenchmarkError, match="already running"):
        benchmark.run_ts_optimizer_case(
            manifest, case_id="case_0001_medium", optimizer="geometric", output=tmp_path / "run",
        )


def test_run_preserves_partial_artifacts_after_interrupted_job(tmp_path, monkeypatch):
    _fake_run_neb(monkeypatch)
    manifest = _manifest(tmp_path)
    job_dir = tmp_path / "run/cases/case_0001_medium/geometric"
    job_dir.mkdir(parents=True)
    partial = job_dir / "ts_path.xyz"
    partial.write_text("partial data must not be overwritten")
    with pytest.raises(benchmark.TSOptimizerBenchmarkError, match="partial or uncollected artifacts"):
        benchmark.run_ts_optimizer_case(
            manifest, case_id="case_0001_medium", optimizer="geometric", output=tmp_path / "run",
        )
    assert partial.read_text() == "partial data must not be overwritten"


def test_failed_optimizer_is_recorded_without_trying_other_optimizer(tmp_path, monkeypatch):
    calls = _fake_run_neb(monkeypatch, fail_optimizer="sella")
    manifest = _manifest(tmp_path)
    result = benchmark.run_ts_optimizer_case(
        manifest, case_id="case_0001_medium", optimizer="sella", output=tmp_path / "run",
    )
    assert result["status"] == "failed"
    assert result["outcome"] == "optimizer_exception"
    assert result["failure_stage"] == "ts"
    assert "fixture optimizer failure" in result["error"]
    assert [stage for stage, _ in calls] == ["ts"]


def test_converged_nonstationary_result_is_distinct_from_optimizer_failure(
    tmp_path, monkeypatch,
):
    _fake_run_neb(monkeypatch, nonstationary=True)
    result = benchmark.run_ts_optimizer_case(
        _manifest(tmp_path), case_id="case_0001_medium", optimizer="geometric",
        output=tmp_path / "run",
    )
    assert result["status"] == "complete"
    assert result["outcome"] == "converged_not_stationary"
    assert result["stages"] == ["ts", "frequency"]


def test_run_preserves_input_pair_hash_and_rejects_overwrite(tmp_path, monkeypatch):
    _fake_run_neb(monkeypatch)
    manifest = _manifest(tmp_path)
    first = benchmark.run_ts_optimizer_case(
        manifest, case_id="case_0001_medium", optimizer="geometric", output=tmp_path / "run",
    )
    assert first["input_sha256"] == first["input_hashes"]["ts_guess"]
    with pytest.raises(benchmark.TSOptimizerBenchmarkError, match="fresh output directory"):
        benchmark.run_ts_optimizer_case(
            manifest, case_id="case_0001_medium", optimizer="geometric", output=tmp_path / "run",
        )
    with pytest.raises(benchmark.TSOptimizerBenchmarkError, match="no optimizer fallback"):
        benchmark.run_ts_optimizer_case(
            manifest, case_id="case_0001_medium", optimizer="native", output=tmp_path / "other",
        )


def test_run_refuses_geometry_changed_after_pair_manifest_creation(tmp_path, monkeypatch):
    _fake_run_neb(monkeypatch)
    guess = tmp_path / "guess.xyz"
    shutil.copy2(DATA / "guess.xyz", guess)
    manifest = _manifest(tmp_path)
    raw = json.loads(manifest.read_text())
    raw["cases"][0]["ts_guess"] = str(guess)
    manifest.write_text(json.dumps(raw))
    benchmark.run_ts_optimizer_case(
        manifest, case_id="case_0001_medium", optimizer="geometric", output=tmp_path / "run",
    )
    shutil.copy2(DATA / "hcn.xyz", guess)
    with pytest.raises(benchmark.TSOptimizerBenchmarkError, match="input files changed"):
        benchmark.run_ts_optimizer_case(
            manifest, case_id="case_0001_medium", optimizer="sella", output=tmp_path / "run",
        )


def test_run_rejects_modified_copied_manifest(tmp_path, monkeypatch):
    _fake_run_neb(monkeypatch)
    manifest = _manifest(tmp_path)
    benchmark.run_ts_optimizer_case(
        manifest, case_id="case_0001_medium", optimizer="geometric", output=tmp_path / "run",
    )
    source_copy = tmp_path / "run" / "benchmark_manifest.json"
    source_copy.write_text("{}")
    with pytest.raises(benchmark.TSOptimizerBenchmarkError, match="differs from the original manifest"):
        benchmark.run_ts_optimizer_case(
            manifest, case_id="case_0001_medium", optimizer="sella", output=tmp_path / "run",
        )


def test_cli_has_independent_run_and_collect_commands():
    args = argument_parse([
        "run", "benchmark.json", "--case", "case1", "--optimizer", "geometric",
        "--output", "benchmark_run",
    ])
    assert (args.command, args.case, args.optimizer, args.output) == (
        "run", "case1", "geometric", "benchmark_run",
    )
    collect = argument_parse(["collect", "benchmark_run"])
    assert (collect.command, collect.output) == ("collect", "benchmark_run")
