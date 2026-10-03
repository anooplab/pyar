import json
import shutil

import pytest

from benchmarks.ts_optimizers.analysis.analyze_orca_optts_comparison import (
    analyze,
    render,
)


def _result(case_id, method, outcome, input_hash="same"):
    result = {
        "case_id": case_id,
        "optimizer": method,
        "outcome": outcome,
        "input_sha256": input_hash,
        "difficulty": case_id.rsplit("_", 1)[-1],
        "ts_optimizer_steps": 10,
        "ts_wall_seconds": 1.0,
    }
    if method == "orca_optts":
        result["ts_orca_gradient_events"] = 9
        result["ts_result"] = {"optimizer_converged": True}
    else:
        result["ts_backend_energy_gradient_evaluations"] = 11
    return result


def _write(directory, relative_path, result):
    path = directory / relative_path
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result), encoding="utf-8")


def _fixture_runs(tmp_path, *, mismatched_hash=False):
    orca = tmp_path / "orca"
    paired = tmp_path / "paired"
    outcomes = {
        "MR_1_easy": {"orca_optts": "reaction_connected_success", "geometric": "failed", "sella": "failed"},
        "MR_1_hard": {"orca_optts": "failed", "geometric": "reaction_connected_success", "sella": "failed"},
        "MR_2_easy": {"orca_optts": "failed", "geometric": "failed", "sella": "failed"},
        "MR_2_hard": {"orca_optts": "failed", "geometric": "failed", "sella": "failed"},
    }
    for case_id, methods in outcomes.items():
        for method, outcome in methods.items():
            digest = "different" if mismatched_hash and case_id == "MR_1_easy" and method == "sella" else "same"
            base = orca if method == "orca_optts" else paired
            relative = (
                f"cases/{case_id}/orca_optts/benchmark_result.json"
                if method == "orca_optts"
                else f"cases/{case_id}/{method}/benchmark_result.json"
            )
            _write(base, relative, _result(case_id, method, outcome, digest))
    return orca, paired


def test_analysis_requires_exactly_matched_input_hashes(tmp_path):
    orca, paired = _fixture_runs(tmp_path, mismatched_hash=True)
    with pytest.raises(ValueError, match="input geometry hash mismatch"):
        analyze(orca, paired)


def test_report_distinguishes_per_start_and_reaction_cluster_results(tmp_path):
    orca, paired = _fixture_runs(tmp_path)
    report = render(analyze(orca, paired))
    assert "Paired starts: 4; reaction clusters: 2" in report
    assert "All three methods used the exact same starting geometry" in report
    assert "ORCA only / PyAR only" in report
    assert "do not establish that the methods are equivalent" in report
    assert "Direct paired optimizer comparison" in report


def test_analysis_accepts_optimizer_results_from_separate_run_directories(tmp_path):
    orca, paired = _fixture_runs(tmp_path)
    geometric = tmp_path / "geometric-only"
    sella = tmp_path / "sella-only"
    for method, destination in (("geometric", geometric), ("sella", sella)):
        for source in paired.glob(f"cases/*/{method}/benchmark_result.json"):
            target = destination / source.relative_to(paired)
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source, target)

    results = analyze(
        orca, paired, geometric_directory=geometric, sella_directory=sella,
    )
    assert len(results["geometric"]) == len(results["sella"]) == 4


def test_report_can_compare_internal_and_cartesian_sella_runs(tmp_path):
    orca, paired = _fixture_runs(tmp_path)
    report = render(analyze(orca, paired, sella_cartesian_directory=paired))

    assert "Sella coordinate-system comparison" in report
    assert "full outcome label changed for 0/4 starts" in report


def test_report_describes_matched_gradient_and_internal_coordinate_protocol(tmp_path):
    orca, paired = _fixture_runs(tmp_path)
    for path in orca.glob("cases/*/orca_optts/benchmark_result.json"):
        result = json.loads(path.read_text(encoding="utf-8"))
        result["settings"] = {
            "reference_pyAR_ts_fmax_applied_to_orca": True,
            "reference_pyAR_ts_fmax_ev_per_angstrom": 0.02,
            "orca_optts_max_cycles": 200,
            "orca_coordinate_system": "redundant_internal",
        }
        result["ts_result"] = {"optimizer_converged": True}
        path.write_text(json.dumps(result), encoding="utf-8")
    for path in paired.glob("cases/*/sella/benchmark_result.json"):
        result = json.loads(path.read_text(encoding="utf-8"))
        result["settings"] = {
            "sella_internal_coordinates": True,
            "ts_fmax": 0.02,
            "ts_max_cycles": 200,
        }
        path.write_text(json.dumps(result), encoding="utf-8")
    for path in paired.glob("cases/*/geometric/benchmark_result.json"):
        result = json.loads(path.read_text(encoding="utf-8"))
        result["settings"] = {"ts_fmax": 0.02, "ts_max_cycles": 200}
        path.write_text(json.dumps(result), encoding="utf-8")

    report = render(analyze(orca, paired))

    assert "ORCA redundant internal coordinates" in report
    assert "geomeTRIC delocalized internals" in report
    assert "Sella internal coordinates" in report
    assert "full convergence definitions are not identical" in report
    assert "All three methods agreed" in report
    assert "Provisional default and recovery recommendation" in report


def test_report_does_not_invent_pilot_measurements(tmp_path):
    orca, paired = _fixture_runs(tmp_path)
    report = render(analyze(orca, paired))
    assert 'Provenance incomplete' in report
    assert '9e-12' not in report
    assert '15-reaction' not in report
    assert 'resamples the 15' not in report


@pytest.mark.parametrize('field', ['Hamiltonian', 'charge', 'imaginary_frequency_threshold', 'reference_ts'])
def test_analysis_rejects_incompatible_physics_and_validation(tmp_path, field):
    orca, paired = _fixture_runs(tmp_path)
    for method in ('geometric', 'sella'):
        path = paired / f'cases/MR_1_easy/{method}/benchmark_result.json'
        result = json.loads(path.read_text())
        mismatch = method == 'sella'
        result['qc_model'] = {'xtb_model': 'gxtb' if mismatch and field == 'Hamiltonian' else 'gfn2'}
        result['case_qc_settings'] = {'charge': int(mismatch and field == 'charge')}
        result['settings'] = {'imaginary_frequency_threshold': 30 if mismatch and field == 'imaginary_frequency_threshold' else 20}
        result['input_hashes'] = {'reference_ts': 'different' if mismatch and field == 'reference_ts' else 'same'}
        path.write_text(json.dumps(result))
    with pytest.raises(ValueError, match='mismatch'):
        analyze(orca, paired)
