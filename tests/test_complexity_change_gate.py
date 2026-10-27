import json
import subprocess
from pathlib import Path

from tools.check_complexity_changes import analyze_complexity_change, main


def _git(repo, *args):
    return subprocess.run(
        ["git", *args], cwd=repo, check=True, capture_output=True, text=True
    ).stdout.strip()


def _repository(tmp_path, source):
    repo = tmp_path / "repo"
    (repo / "pyar").mkdir(parents=True)
    (repo / "pyar" / "__init__.py").write_text("", encoding="utf-8")
    (repo / "pyar" / "sample.py").write_text(source, encoding="utf-8")
    _git(repo, "init", "-q")
    _git(repo, "config", "user.name", "PyAR test")
    _git(repo, "config", "user.email", "pyar-test@example.invalid")
    _git(repo, "add", "pyar")
    _git(repo, "commit", "-qm", "baseline")
    return repo, _git(repo, "rev-parse", "HEAD")


def _commit_source(repo, source):
    (repo / "pyar" / "sample.py").write_text(source, encoding="utf-8")
    _git(repo, "add", "pyar/sample.py")
    _git(repo, "commit", "-qm", "candidate")
    return _git(repo, "rev-parse", "HEAD")


def test_no_change_has_zero_decision_delta(tmp_path):
    repo, sha = _repository(tmp_path, "def work(value):\n    return value\n")

    report = analyze_complexity_change(repo, sha, sha, enforce=True)

    assert report["changed_files"] == []
    assert report["decision_delta"] == 0
    assert report["changed_functions"] == {"modified": [], "added": [], "removed": []}
    assert report["status"] == "pass"
    assert "unavailable" in report["coverage_status"]


def test_new_decision_requires_explicit_pr_reason_and_owner(tmp_path):
    repo, base = _repository(tmp_path, "def work(value):\n    return value\n")
    head = _commit_source(
        repo,
        "def work(value):\n    if value:\n        return 1\n    return 0\n",
    )
    event = tmp_path / "event.json"
    event.write_text(json.dumps({"pull_request": {"body": "No exception."}}))

    report = analyze_complexity_change(repo, base, head, event_json=event, enforce=True)

    assert report["decision_delta"] == 1
    assert len(report["changed_functions"]["modified"]) == 1
    assert report["status"] == "fail"


def test_cli_returns_failure_for_unexplained_enforced_growth(tmp_path):
    repo, base = _repository(tmp_path, "def work(value):\n    return value\n")
    head = _commit_source(
        repo,
        "def work(value):\n    if value:\n        return 1\n    return 0\n",
    )
    output = tmp_path / "reports" / "complexity-change.json"

    exit_code = main(
        [
            "--repo",
            str(repo),
            "--base",
            base,
            "--head",
            head,
            "--output",
            str(output),
            "--enforce",
        ]
    )

    assert exit_code == 1
    assert json.loads(output.read_text())["status"] == "fail"


def test_decision_preserving_helper_extraction_is_not_metric_gamed(tmp_path):
    repo, base = _repository(
        tmp_path,
        "def work(value):\n    if value:\n        return 1\n    return 0\n",
    )
    head = _commit_source(
        repo,
        "def work(value):\n    return _choose(value)\n\n"
        "def _choose(value):\n    if value:\n        return 1\n    return 0\n",
    )

    report = analyze_complexity_change(repo, base, head, enforce=True)

    assert report["base_decision_total"] == report["head_decision_total"] == 1
    assert report["decision_delta"] == 0
    assert len(report["changed_functions"]["modified"]) == 1
    assert len(report["changed_functions"]["added"]) == 1
    assert report["status"] == "pass"


def test_positive_delta_passes_with_explicit_reason_and_owner(tmp_path):
    repo, base = _repository(tmp_path, "def work(value):\n    return value\n")
    head = _commit_source(
        repo,
        "def work(value):\n    if value:\n        return 1\n    return 0\n",
    )
    event = tmp_path / "event.json"
    event.write_text(
        json.dumps(
            {
                "pull_request": {
                    "body": "Complexity-Increase-Reason: adds explicit input policy\n"
                    "Complexity-Increase-Owner: @maintainer"
                }
            }
        )
    )

    report = analyze_complexity_change(repo, base, head, event_json=event, enforce=True)

    assert report["decision_delta"] == 1
    assert report["status"] == "pass_with_review_note"
    assert report["exception"] == {
        "reason": "adds explicit input policy",
        "owner": "@maintainer",
    }


def test_removing_decisions_reduces_budget_without_exception(tmp_path):
    repo, base = _repository(
        tmp_path,
        "def work(value):\n    if value:\n        return 1\n    return 0\n",
    )
    head = _commit_source(repo, "def work(value):\n    return int(bool(value))\n")

    report = analyze_complexity_change(repo, base, head, enforce=True)

    assert report["decision_delta"] == -1
    assert report["status"] == "pass"


def test_changed_function_source_is_reported_even_when_cc_is_unchanged(tmp_path):
    repo, base = _repository(tmp_path, "def work(value):\n    return value + 1\n")
    head = _commit_source(repo, "def work(value):\n    return value + 2\n")

    report = analyze_complexity_change(repo, base, head)

    assert report["decision_delta"] == 0
    assert len(report["changed_functions"]["modified"]) == 1
    before = report["changed_functions"]["modified"][0]["before"]
    after = report["changed_functions"]["modified"][0]["after"]
    assert before["source_sha256"] != after["source_sha256"]
    assert before["cc_radon"] == after["cc_radon"] == 1


def test_renamed_function_is_reported_as_removed_and_added(tmp_path):
    repo, base = _repository(tmp_path, "def old_name(value):\n    return value\n")
    head = _commit_source(repo, "def new_name(value):\n    return value\n")

    report = analyze_complexity_change(repo, base, head)

    assert report["decision_delta"] == 0
    assert [
        row["qualified_name"] for row in report["changed_functions"]["removed"]
    ] == ["old_name"]
    assert [row["qualified_name"] for row in report["changed_functions"]["added"]] == [
        "new_name"
    ]
