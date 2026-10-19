import hashlib
import json
from pathlib import Path

import pytest

from tools.scientific_validation import EvidenceError, evaluate, main
from pyar.backend_capabilities import BACKEND_CAPABILITIES

pytestmark = [pytest.mark.unit, pytest.mark.regression]


@pytest.fixture
def manifest(tmp_path):
    reference = tmp_path / "reference.dat"
    reference.write_text("independent reference", encoding="utf-8")
    digest = hashlib.sha256(reference.read_bytes()).hexdigest()
    input_path = tmp_path / "input.xyz"
    input_path.write_text("1\nH\nH 0 0 0\n", encoding="utf-8")
    input_digest = hashlib.sha256(input_path.read_bytes()).hexdigest()
    log_path = tmp_path / "job" / "output.log"
    log_path.parent.mkdir()
    log_path.write_text("recorded raw backend output\n", encoding="utf-8")
    return {
        "schema_version": 1,
        "benchmark_id": "test-case",
        "cases": [
            {
                "id": "h2-single-point",
                "tier": "tier2",
                "system": "H2",
                "charge": 0,
                "multiplicity": 1,
                "method": "GFN2-xTB",
                "backend": "xtb",
                "input": "input.xyz",
                "input_sha256": input_digest,
                "reference": {
                    "source": "independent test source",
                    "version": "1",
                    "provenance": "fixture created independently of the evaluated result",
                    "path": "reference.dat",
                    "sha256": digest,
                },
                "acceptance": {
                    "energy": {"expected": -1.0, "tolerance": 0.01, "unit": "hartree"}
                },
            }
        ],
    }


def result(**updates):
    row = {
        "id": "h2-single-point",
        "status": "PASS",
        "charge": 0,
        "multiplicity": 1,
        "method": "GFN2-xTB",
        "backend": "xtb",
        "backend_version": "test 1.0",
        "input_sha256": hashlib.sha256(b"1\nH\nH 0 0 0\n").hexdigest(),
        "metrics": {"energy": -1.002},
        "logs": ["job/output.log"],
    }
    row.update(updates)
    return {"schema_version": 1, "candidate_sha": "b" * 40, "cases": [row]}


def test_manifest_backed_complete_result_passes(manifest, tmp_path):
    report = evaluate(manifest, result(), tmp_path, expected_candidate_sha="b" * 40)
    assert report["verdict"] == "PASS"
    assert report["counts"]["PASS"] == 1
    assert report["cases"][0]["logs"] == [
        {
            "path": "job/output.log",
            "sha256": hashlib.sha256(b"recorded raw backend output\n").hexdigest(),
        }
    ]


@pytest.mark.parametrize(
    "updates",
    [
        {"metrics": {"energy": -1.2}},
        {"method": "GFN1-xTB"},
        {"input_sha256": "c" * 64},
        {"backend_version": None},
        {"logs": []},
    ],
)
def test_invalid_scientific_evidence_fails_closed(manifest, tmp_path, updates):
    report = evaluate(
        manifest, result(**updates), tmp_path, expected_candidate_sha="b" * 40
    )
    assert report["verdict"] == "FAIL"
    assert report["counts"]["FAIL"] == 1


def test_missing_case_is_not_tested_and_blocks_verdict(manifest, tmp_path):
    report = evaluate(
        manifest,
        {"schema_version": 1, "candidate_sha": "b" * 40, "cases": []},
        tmp_path,
        expected_candidate_sha="b" * 40,
    )
    assert report["cases"][0]["status"] == "NOT_TESTED"
    assert report["verdict"] == "FAIL"


def test_not_applicable_without_authorized_exception_is_blocked(manifest, tmp_path):
    report = evaluate(
        manifest,
        result(status="NOT_APPLICABLE"),
        tmp_path,
        expected_candidate_sha="b" * 40,
    )
    assert report["cases"][0]["status"] == "BLOCKED"
    assert report["verdict"] == "FAIL"


def test_approved_not_applicable_case_is_reported_as_exception(manifest, tmp_path):
    report = evaluate(
        manifest,
        result(
            status="NOT_APPLICABLE",
            authorized_exception={
                "reason": "route excluded from this release",
                "authorized_by": "scientific owner",
                "approval_reference": "decision-123",
            },
        ),
        tmp_path,
        expected_candidate_sha="b" * 40,
    )
    assert report["cases"][0]["status"] == "NOT_APPLICABLE"
    assert report["verdict"] == "PASS"


def test_changed_reference_hash_is_rejected(manifest, tmp_path):
    (tmp_path / "reference.dat").write_text("changed", encoding="utf-8")
    with pytest.raises(EvidenceError, match="hash differs"):
        evaluate(manifest, result(), tmp_path, expected_candidate_sha="b" * 40)


def test_candidate_sha_must_match_approved_head(manifest, tmp_path):
    with pytest.raises(EvidenceError, match="differs from expected candidate"):
        evaluate(manifest, result(), tmp_path, expected_candidate_sha="c" * 40)


def test_missing_raw_log_cannot_pass(manifest, tmp_path):
    (tmp_path / "job/output.log").unlink()
    report = evaluate(manifest, result(), tmp_path, expected_candidate_sha="b" * 40)
    assert report["verdict"] == "FAIL"
    assert "raw log artifact is missing" in report["cases"][0]["failures"][0]


def test_candidate_sha_must_be_hex(manifest, tmp_path):
    invalid = result()
    invalid["candidate_sha"] = "z" * 40
    with pytest.raises(EvidenceError, match="full 40-character Git SHA"):
        evaluate(manifest, invalid, tmp_path, expected_candidate_sha="b" * 40)


def test_backend_matrix_covers_current_capability_registry():
    matrix_path = (
        Path(__file__).parents[1]
        / "benchmarks/scientific_validation/backend_matrix.json"
    )
    matrix = json.loads(matrix_path.read_text(encoding="utf-8"))
    # Some older integration tests register synthetic backends globally and
    # leave them behind; exclude those known test-only names from this inventory.
    test_only_names = {"fake", "missing_provider"}
    assert {row["backend"] for row in matrix["entries"]} == (
        set(BACKEND_CAPABILITIES) - test_only_names
    )
    for row in matrix["entries"]:
        capability = BACKEND_CAPABILITIES[row["backend"]]
        assert row["required_executables"] == sorted(capability.required_executables)
        assert row["required_python_modules"] == sorted(
            capability.required_python_modules
        )
        assert row["real_calculation"]["status"] in {"NOT_TESTED", "PASS", "FAIL"}


def test_cli_writes_machine_and_human_reports_for_missing_case(
    manifest, tmp_path, monkeypatch
):
    manifest_path = tmp_path / "cases.json"
    results_path = tmp_path / "results.json"
    json_output = tmp_path / "report.json"
    markdown_output = tmp_path / "report.md"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    results_path.write_text(
        json.dumps({"schema_version": 1, "candidate_sha": "b" * 40, "cases": []}),
        encoding="utf-8",
    )
    monkeypatch.setattr(
        "sys.argv",
        [
            "scientific_validation.py",
            str(manifest_path),
            str(results_path),
            "--candidate-sha",
            "b" * 40,
            "--json",
            str(json_output),
            "--markdown",
            str(markdown_output),
        ],
    )
    assert main() == 1
    assert '"verdict": "FAIL"' in json_output.read_text(encoding="utf-8")
    assert "NOT_TESTED" in markdown_output.read_text(encoding="utf-8")
