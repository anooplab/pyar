"""Fail-closed checks for promoting approved PyAR release artifacts."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import urllib.error
import urllib.request


class PromotionError(ValueError):
    """Raised when promotion evidence does not match the approved release."""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_run(
    record: dict,
    *,
    run_id: str,
    workflow: str,
    sha: str,
    required_jobs: list[str],
) -> None:
    if str(record.get("id")) != str(run_id):
        raise PromotionError(f"workflow run id differs from approved run {run_id}")
    if record.get("path") != workflow:
        raise PromotionError(f"run {run_id} is not from {workflow}")
    if record.get("head_sha") != sha:
        raise PromotionError(f"run {run_id} is for a different candidate SHA")
    if record.get("status") != "completed" or record.get("conclusion") != "success":
        raise PromotionError(f"run {run_id} has not completed successfully")
    jobs = record.get("jobs")
    if not isinstance(jobs, list) or not jobs:
        raise PromotionError(f"run {run_id} has no recorded jobs")
    failed = [
        job.get("name", "unnamed") for job in jobs if job.get("conclusion") != "success"
    ]
    if failed:
        raise PromotionError(f"run {run_id} has non-success jobs: {', '.join(failed)}")
    observed_jobs = {job.get("name") for job in jobs}
    if observed_jobs != set(required_jobs):
        raise PromotionError(
            f"run {run_id} job inventory differs from the approved list"
        )
    if workflow.endswith("/ci.yml") and "R1 required checks" not in observed_jobs:
        raise PromotionError("CI run is missing the R1 aggregate required-check job")
    if (
        workflow.endswith("/scientific-tier2.yml")
        and "Real xTB GFN2 scan smoke" not in observed_jobs
    ):
        raise PromotionError("Tier 2 run is missing the real-backend smoke job")


def validate_release(
    approval: dict,
    dist_dir: Path,
    *,
    version: str,
    candidate_sha: str,
    ci_run: dict,
    tier2_run: dict,
) -> dict[str, str]:
    if not re.fullmatch(r"[0-9]+\.[0-9]+\.[0-9]+", version):
        raise PromotionError("release version must be a three-part numeric version")
    if approval.get("schema_version") != 1:
        raise PromotionError("unsupported release approval schema")
    if approval.get("status") != "APPROVED":
        raise PromotionError("release approval is not APPROVED")
    if approval.get("version") != version:
        raise PromotionError("release version differs from approval")
    if approval.get("candidate_sha") != candidate_sha:
        raise PromotionError("candidate SHA differs from approval")
    if not re.fullmatch(r"[0-9a-f]{40}", candidate_sha):
        raise PromotionError("candidate SHA must be a full lowercase Git SHA")
    if approval.get("target_registry") != "https://upload.pypi.org/legacy/":
        raise PromotionError("target registry is not the approved PyPI registry")
    if not approval.get("approved_by") or not approval.get("approved_at_utc"):
        raise PromotionError("owner approval identity or timestamp is missing")
    if approval.get("scientific_scope_status") != "PASS_SCOPED":
        raise PromotionError("approved scientific scope is not marked PASS_SCOPED")
    if not approval.get("scientific_exclusions"):
        raise PromotionError("scientific exclusions are not recorded")

    run_ids = approval.get("workflow_runs", {})
    required_jobs = approval.get("required_jobs", {})
    validate_run(
        ci_run,
        run_id=str(run_ids.get("ci")),
        workflow=".github/workflows/ci.yml",
        sha=candidate_sha,
        required_jobs=required_jobs.get("ci", []),
    )
    validate_run(
        tier2_run,
        run_id=str(run_ids.get("tier2")),
        workflow=".github/workflows/scientific-tier2.yml",
        sha=candidate_sha,
        required_jobs=required_jobs.get("tier2", []),
    )

    expected = approval.get("artifacts")
    if not isinstance(expected, dict) or len(expected) != 2:
        raise PromotionError("approval must contain exactly wheel and sdist hashes")
    actual_paths = sorted(
        path
        for path in dist_dir.iterdir()
        if path.is_file()
        and (path.name.endswith(".whl") or path.name.endswith(".tar.gz"))
    )
    if len(actual_paths) != 2 or {p.name for p in actual_paths} != set(expected):
        raise PromotionError("downloaded artifact set differs from approval")
    actual = {path.name: sha256(path) for path in actual_paths}
    if actual != expected:
        raise PromotionError("downloaded artifact hashes differ from owner approval")
    expected_names = {
        f"pyar_chem-{version}-py3-none-any.whl",
        f"pyar_chem-{version}.tar.gz",
    }
    if set(expected) != expected_names:
        raise PromotionError("approved distribution names do not match the version")
    return actual


def assert_unpublished(version: str) -> None:
    url = f"https://pypi.org/pypi/pyar-chem/{version}/json"
    try:
        with urllib.request.urlopen(url, timeout=15):
            raise PromotionError(f"PyPI version {version} already exists")
    except urllib.error.HTTPError as exc:
        if exc.code != 404:
            raise PromotionError(
                f"cannot verify PyPI version availability: HTTP {exc.code}"
            ) from exc
    except urllib.error.URLError as exc:
        raise PromotionError(
            f"cannot reach PyPI to check version availability: {exc}"
        ) from exc


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("approval", type=Path)
    parser.add_argument("dist", type=Path)
    parser.add_argument("--version", required=True)
    parser.add_argument("--candidate-sha", required=True)
    parser.add_argument("--ci-run", type=Path, required=True)
    parser.add_argument("--tier2-run", type=Path, required=True)
    parser.add_argument("--assert-unpublished", action="store_true")
    args = parser.parse_args()
    approval = json.loads(args.approval.read_text(encoding="utf-8"))
    ci_run = json.loads(args.ci_run.read_text(encoding="utf-8"))
    tier2_run = json.loads(args.tier2_run.read_text(encoding="utf-8"))
    hashes = validate_release(
        approval,
        args.dist,
        version=args.version,
        candidate_sha=args.candidate_sha,
        ci_run=ci_run,
        tier2_run=tier2_run,
    )
    if args.assert_unpublished:
        assert_unpublished(args.version)
    print(
        json.dumps(
            {
                "status": "PASS",
                "version": args.version,
                "candidate_sha": args.candidate_sha,
                "artifacts": hashes,
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
