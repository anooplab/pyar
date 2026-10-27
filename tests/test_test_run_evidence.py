import json
from pathlib import Path
import subprocess
import xml.etree.ElementTree as ET

import pytest

from tools import run_tests_with_evidence as evidence


def _write_junit(path: Path, *, declared: int, testcase_count: int) -> None:
    suite = ET.Element(
        "testsuite", tests=str(declared), failures="0", errors="0", skipped="0"
    )
    for index in range(testcase_count):
        ET.SubElement(suite, "testcase", classname="test_example", name=f"test_{index}")
    ET.ElementTree(suite).write(path, encoding="utf-8", xml_declaration=True)


def test_junit_summary_preserves_declared_and_emitted_counts(tmp_path):
    report = tmp_path / "junit.xml"
    _write_junit(report, declared=7, testcase_count=4)

    assert evidence._junit_summary(report) == {
        "declared_tests": 7,
        "testcase_elements": 4,
        "declared_matches_testcases": False,
        "failures": 0,
        "errors": 0,
        "skipped": 0,
        "suite_timestamps": [None],
    }


def test_evidence_run_binds_reports_and_does_not_rewrite_count_mismatch(
    tmp_path, monkeypatch
):
    report_payload = b'<testsuite tests="7" failures="0" errors="0" skipped="0"><testcase/><testcase/><testcase/><testcase/></testsuite>'

    def fake_subprocess(command, *, cwd, env, check, **kwargs):
        assert command[1:3] == ["-m", "pytest"]
        junit = cwd / next(
            arg.split("=", 1)[1] for arg in command if arg.startswith("--junitxml=")
        )
        coverage = cwd / next(
            arg.split("--cov-report=json:", 1)[1]
            for arg in command
            if arg.startswith("--cov-report=json:")
        )
        junit.write_bytes(report_payload)
        coverage.write_text(
            json.dumps(
                {
                    "meta": {
                        "timestamp": evidence._utc_now().isoformat(),
                        "version": "7.16.1",
                        "branch_coverage": True,
                    }
                }
            ),
            encoding="utf-8",
        )
        return subprocess.CompletedProcess(
            command, 0, stdout="3 subtests passed\n", stderr=""
        )

    monkeypatch.setattr(evidence.subprocess, "run", fake_subprocess)
    monkeypatch.setattr(evidence, "_git_sha", lambda root: "source-sha")
    run_root = tmp_path / "runs"

    assert evidence.run(tmp_path, run_root, ["-q"]) == 0

    run_dir = next(path for path in run_root.iterdir() if path.is_dir())
    manifest = json.loads((run_dir / "manifest.json").read_text())
    assert manifest["source_sha"] == "source-sha"
    assert manifest["junit"]["declared_tests"] == 7
    assert manifest["junit"]["testcase_elements"] == 4
    assert manifest["junit"]["declared_matches_testcases"] is False
    assert manifest["junit"]["subtest_count_reconciles"] is True
    assert manifest["pytest_summary"] == {"passed": 3, "failed": 0, "skipped": 0}
    assert manifest["artifacts"]["junit"]["sha256"] == evidence._sha256(
        run_dir / "junit.xml"
    )
    assert "per-subtest XML records are absent" in manifest["limitations"][0]
    assert (run_dir / "junit.xml").read_bytes() == report_payload


def test_failed_pytest_exit_code_and_missing_reports_are_recorded(
    tmp_path, monkeypatch
):
    monkeypatch.setattr(
        evidence.subprocess,
        "run",
        lambda *args, **kwargs: subprocess.CompletedProcess(
            args[0], 2, stdout="", stderr=""
        ),
    )
    monkeypatch.setattr(evidence, "_git_sha", lambda root: "source-sha")

    assert evidence.run(tmp_path, tmp_path / "runs", []) == 2
    run_dir = next((tmp_path / "runs").iterdir())
    manifest = json.loads((run_dir / "manifest.json").read_text())
    assert manifest["pytest_exit_code"] == 2
    assert manifest["runner_exit_code"] == 2
    assert set(manifest["artifacts"]) == {"pytest_stdout", "pytest_stderr"}
    assert set(manifest["limitations"]) == {
        "junit report was not produced",
        "coverage report was not produced",
    }


def test_malformed_report_is_recorded_as_incomplete_evidence(tmp_path, monkeypatch):
    def fake_subprocess(command, **kwargs):
        junit = kwargs["cwd"] / next(
            arg.split("=", 1)[1] for arg in command if arg.startswith("--junitxml=")
        )
        coverage = kwargs["cwd"] / next(
            arg.split("--cov-report=json:", 1)[1]
            for arg in command
            if arg.startswith("--cov-report=json:")
        )
        junit.write_text("<testsuite", encoding="utf-8")
        coverage.write_text("{}", encoding="utf-8")
        return subprocess.CompletedProcess(command, 0, stdout="", stderr="")

    monkeypatch.setattr(evidence.subprocess, "run", fake_subprocess)
    monkeypatch.setattr(evidence, "_git_sha", lambda root: "source-sha")

    assert evidence.run(tmp_path, tmp_path / "runs", []) == 2
    run_dir = next((tmp_path / "runs").iterdir())
    manifest = json.loads((run_dir / "manifest.json").read_text())
    assert manifest["evidence_complete"] is False
    assert manifest["junit"]["valid"] is False
    assert manifest["coverage"]["valid"] is False
    assert sum(
        message.startswith(("JUnit report is invalid:", "coverage report is invalid:"))
        for message in manifest["evidence_errors"]
    ) == 2


@pytest.mark.parametrize("bad_xml", ["", "<testsuite", "<not-junit/>"])
def test_junit_summary_rejects_invalid_or_missing_suite(tmp_path, bad_xml):
    report = tmp_path / "junit.xml"
    report.write_text(bad_xml, encoding="utf-8")
    with pytest.raises((ET.ParseError, ValueError)):
        evidence._junit_summary(report)
