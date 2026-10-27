"""Run pytest into an isolated, source-bound evidence directory."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
from importlib.metadata import version
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import uuid
import xml.etree.ElementTree as ET


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _utc_now() -> datetime:
    return datetime.now(timezone.utc)


def _parse_timestamp(value: str) -> datetime:
    parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    if parsed.tzinfo is None:
        parsed = parsed.astimezone()
    return parsed.astimezone(timezone.utc)


def _junit_summary(path: Path) -> dict[str, object]:
    root = ET.parse(path).getroot()
    suites = [root] if root.tag == "testsuite" else list(root.findall("testsuite"))
    if not suites:
        raise ValueError("JUnit XML contains no top-level testsuite")
    cases = list(root.iter("testcase"))
    declared = sum(int(suite.attrib.get("tests", "0")) for suite in suites)
    return {
        "declared_tests": declared,
        "testcase_elements": len(cases),
        "declared_matches_testcases": declared == len(cases),
        "failures": sum(int(suite.attrib.get("failures", "0")) for suite in suites),
        "errors": sum(int(suite.attrib.get("errors", "0")) for suite in suites),
        "skipped": sum(int(suite.attrib.get("skipped", "0")) for suite in suites),
        "suite_timestamps": [suite.attrib.get("timestamp") for suite in suites],
    }


def _git_sha(root: Path) -> str:
    result = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=root, text=True, capture_output=True
    )
    if result.returncode:
        raise RuntimeError(f"could not identify source commit: {result.stderr.strip()}")
    return result.stdout.strip()


def _source_tree_sha256(root: Path) -> str:
    files = sorted(
        path
        for path in (root / "pyar").rglob("*.py")
        if "AIMNet2" not in path.relative_to(root).parts
        and "__pycache__" not in path.parts
    )
    digest = hashlib.sha256()
    for path in files:
        relative = path.relative_to(root).as_posix()
        file_digest = _sha256(path)
        digest.update(relative.encode("utf-8"))
        digest.update(b"\0")
        digest.update(file_digest.encode("ascii"))
        digest.update(b"\n")
    return digest.hexdigest()


def _subtest_summary(output: str) -> dict[str, int] | None:
    result = {"passed": 0, "failed": 0, "skipped": 0}
    found = False
    for outcome in result:
        match = re.search(rf"\b(\d+) subtests? {outcome}\b", output)
        if match:
            result[outcome] = int(match.group(1))
            found = True
    return result if found else None


def run(root: Path, output_root: Path, pytest_args: list[str]) -> int:
    run_id = uuid.uuid4().hex
    run_dir = output_root / run_id
    run_dir.mkdir(parents=True, exist_ok=False)
    junit = run_dir / "junit.xml"
    coverage = run_dir / "coverage.json"
    coverage_data = run_dir / ".coverage"
    source_tree_before = _source_tree_sha256(root)
    started = _utc_now()
    junit_arg = os.path.relpath(junit, root)
    coverage_arg = os.path.relpath(coverage, root)
    coverage_data_arg = os.path.relpath(coverage_data, root)
    command = [
        sys.executable,
        "-m",
        "pytest",
        "--cov=pyar",
        f"--cov-report=json:{coverage_arg}",
        f"--junitxml={junit_arg}",
        *pytest_args,
    ]
    environment = os.environ.copy()
    environment["COVERAGE_FILE"] = str(root / coverage_data_arg)
    environment["PYAR_TEST_RUN_ID"] = run_id
    process = subprocess.run(
        command, cwd=root, env=environment, check=False, capture_output=True, text=True
    )
    finished = _utc_now()
    stdout = process.stdout or ""
    stderr = process.stderr or ""
    (run_dir / "pytest.stdout.txt").write_text(stdout, encoding="utf-8")
    (run_dir / "pytest.stderr.txt").write_text(stderr, encoding="utf-8")
    sys.stdout.write(stdout)
    sys.stderr.write(stderr)

    manifest: dict[str, object] = {
        "schema": 1,
        "run_id": run_id,
        "source_sha": _git_sha(root),
        "source_tree_sha256_before": source_tree_before,
        "source_tree_sha256_after": _source_tree_sha256(root),
        "command": command,
        "cwd": ".",
        "python_version": sys.version.split()[0],
        "pytest_version": version("pytest"),
        "evidence_runner_sha256": _sha256(Path(__file__).resolve()),
        "coverage_data_path": coverage_data_arg,
        "started_utc": started.isoformat(),
        "finished_utc": finished.isoformat(),
        "pytest_exit_code": process.returncode,
        "runner_exit_code": process.returncode,
        "artifacts": {},
        "junit": None,
        "evidence_errors": [],
        "pytest_summary": _subtest_summary(stdout + "\n" + stderr),
        "coverage": None,
        "limitations": [],
    }
    if manifest["source_tree_sha256_before"] != manifest["source_tree_sha256_after"]:
        message = "first-party production source tree changed while pytest was running"
        manifest["limitations"].append(message)
        manifest["evidence_errors"].append(message)
    artifacts: dict[str, dict[str, object]] = {}
    for label, path in (
        ("junit", junit),
        ("coverage", coverage),
        ("pytest_stdout", run_dir / "pytest.stdout.txt"),
        ("pytest_stderr", run_dir / "pytest.stderr.txt"),
    ):
        if path.is_file():
            stat = path.stat()
            artifacts[label] = {
                "path": path.name,
                "sha256": _sha256(path),
                "size_bytes": stat.st_size,
                "mtime_utc": datetime.fromtimestamp(
                    stat.st_mtime, timezone.utc
                ).isoformat(),
            }
        else:
            if label in {"junit", "coverage"}:
                message = f"{label} report was not produced"
                manifest["limitations"].append(message)
                manifest["evidence_errors"].append(message)
    manifest["artifacts"] = artifacts

    if junit.is_file():
        try:
            manifest["junit"] = _junit_summary(junit)
            manifest["junit"]["valid"] = True
            if manifest["pytest_summary"] is not None:
                subtests = sum(manifest["pytest_summary"].values())
                mismatch = (
                    manifest["junit"]["declared_tests"]
                    - manifest["junit"]["testcase_elements"]
                )
                manifest["junit"]["subtest_count_reconciles"] = mismatch == subtests
                if mismatch and mismatch != subtests:
                    manifest["evidence_errors"].append(
                        "JUnit count difference is not reconciled by pytest subtest summary"
                    )
            if not manifest["junit"]["declared_matches_testcases"]:
                manifest["limitations"].append(
                    "JUnit groups unittest subtest outcomes under parent testcases; "
                    "per-subtest XML records are absent and the declared count must "
                    "not be rewritten for parsers that require one testcase per subtest"
                )
        except (ET.ParseError, ValueError) as error:
            message = f"JUnit report is invalid: {error}"
            manifest["junit"] = {"valid": False, "error": str(error)}
            manifest["limitations"].append(message)
            manifest["evidence_errors"].append(message)
    if coverage.is_file():
        try:
            data = json.loads(coverage.read_text(encoding="utf-8"))
            coverage_time = _parse_timestamp(data["meta"]["timestamp"])
            manifest["coverage"] = {
                "valid": True,
                "timestamp": data["meta"]["timestamp"],
                "coverage_version": data["meta"]["version"],
                "branch_coverage": data["meta"]["branch_coverage"],
            }
            if not (started <= coverage_time <= finished):
                message = "coverage timestamp is outside the pytest process interval"
                manifest["limitations"].append(message)
                manifest["evidence_errors"].append(message)
        except (json.JSONDecodeError, KeyError, TypeError, ValueError) as error:
            message = f"coverage report is invalid: {error}"
            manifest["coverage"] = {"valid": False, "error": str(error)}
            manifest["limitations"].append(message)
            manifest["evidence_errors"].append(message)
    for label in ("junit", "coverage"):
        artifact = artifacts.get(label)
        if artifact is not None:
            mtime = _parse_timestamp(str(artifact["mtime_utc"]))
            if not (started <= mtime <= finished):
                manifest["limitations"].append(
                    f"{label} artifact modification time is outside the pytest process interval"
                )
                manifest["evidence_errors"].append(
                    f"{label} artifact modification time is outside the pytest process interval"
                )
    manifest["evidence_complete"] = not manifest["evidence_errors"]
    if process.returncode == 0 and not manifest["evidence_complete"]:
        manifest["runner_exit_code"] = 2
    (run_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"Evidence manifest: {run_dir / 'manifest.json'}")
    return manifest["runner_exit_code"]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--output-root", type=Path, default=Path("reports/test-runs"))
    parser.add_argument("pytest_args", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    root = args.root.resolve()
    output_root = args.output_root
    if not output_root.is_absolute():
        output_root = root / output_root
    pytest_args = args.pytest_args
    if pytest_args[:1] == ["--"]:
        pytest_args = pytest_args[1:]
    return run(root, output_root, pytest_args)


if __name__ == "__main__":
    raise SystemExit(main())
