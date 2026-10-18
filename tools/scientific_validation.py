#!/usr/bin/env python3
"""Validate scientific benchmark evidence against an immutable case manifest.

This tool evaluates recorded evidence; it does not run chemistry or infer that
an unavailable backend has passed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

ALLOWED_STATUS = {"PASS", "FAIL", "NOT_TESTED", "NOT_APPLICABLE", "BLOCKED"}
REQUIRED_CASE_FIELDS = {
    "id",
    "tier",
    "system",
    "charge",
    "multiplicity",
    "method",
    "backend",
    "input",
    "input_sha256",
    "reference",
    "acceptance",
}


class EvidenceError(ValueError):
    """Raised for malformed or unverifiable benchmark evidence."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_manifest(manifest: dict, base_dir: Path) -> list[dict]:
    if manifest.get("schema_version") != 1:
        raise EvidenceError("manifest schema_version must be 1")
    cases = manifest.get("cases")
    if not isinstance(cases, list) or not cases:
        raise EvidenceError("manifest must contain at least one case")
    seen = set()
    for case in cases:
        missing = REQUIRED_CASE_FIELDS - case.keys()
        if missing:
            raise EvidenceError(
                f"case is missing required fields: {', '.join(sorted(missing))}"
            )
        if case["id"] in seen:
            raise EvidenceError(f"duplicate case id: {case['id']}")
        seen.add(case["id"])
        if case["tier"] not in {"tier1", "tier2", "tier3"}:
            raise EvidenceError(f"{case['id']}: unsupported tier {case['tier']!r}")
        for key in ("charge", "multiplicity"):
            if not isinstance(case[key], int):
                raise EvidenceError(f"{case['id']}: {key} must be an integer")
        if case["multiplicity"] < 1:
            raise EvidenceError(f"{case['id']}: multiplicity must be positive")
        if (
            not isinstance(case["input_sha256"], str)
            or len(case["input_sha256"]) != 64
            or any(
                char not in "0123456789abcdef" for char in case["input_sha256"].lower()
            )
        ):
            raise EvidenceError(
                f"{case['id']}: input_sha256 must be a SHA-256 hex digest"
            )
        input_path = (base_dir / case["input"]).resolve()
        if base_dir.resolve() not in input_path.parents:
            raise EvidenceError(
                f"{case['id']}: input path must stay within manifest directory"
            )
        if not input_path.is_file() or sha256_file(input_path) != case["input_sha256"]:
            raise EvidenceError(
                f"{case['id']}: input file is missing or its hash differs"
            )
        reference = case["reference"]
        if not isinstance(reference, dict) or not all(
            reference.get(k) for k in ("source", "version", "provenance", "sha256")
        ):
            raise EvidenceError(
                f"{case['id']}: reference source, version, provenance, and sha256 are required"
            )
        if len(reference["sha256"]) != 64 or any(
            char not in "0123456789abcdef" for char in reference["sha256"].lower()
        ):
            raise EvidenceError(
                f"{case['id']}: reference sha256 must be a SHA-256 hex digest"
            )
        ref_path = reference.get("path")
        if ref_path:
            path = (base_dir / ref_path).resolve()
            if base_dir.resolve() not in path.parents:
                raise EvidenceError(
                    f"{case['id']}: reference path must stay within manifest directory"
                )
            if not path.is_file() or sha256_file(path) != reference["sha256"]:
                raise EvidenceError(
                    f"{case['id']}: reference file is missing or its hash differs"
                )
        acceptance = case["acceptance"]
        if not isinstance(acceptance, dict) or not acceptance:
            raise EvidenceError(
                f"{case['id']}: at least one acceptance metric is required"
            )
        for metric, rule in acceptance.items():
            if not isinstance(rule, dict) or not isinstance(
                rule.get("expected"), (int, float)
            ):
                raise EvidenceError(
                    f"{case['id']}.{metric}: numeric expected value is required"
                )
            tolerance = rule.get("tolerance")
            if (
                not isinstance(tolerance, (int, float))
                or tolerance < 0
                or not math.isfinite(tolerance)
            ):
                raise EvidenceError(
                    f"{case['id']}.{metric}: finite nonnegative tolerance is required"
                )
            if not rule.get("unit"):
                raise EvidenceError(f"{case['id']}.{metric}: unit is required")
    return cases


def evaluate(manifest: dict, results: dict, base_dir: Path) -> dict:
    cases = validate_manifest(manifest, base_dir)
    if results.get("schema_version") != 1:
        raise EvidenceError("results schema_version must be 1")
    candidate = results.get("candidate_sha")
    if not isinstance(candidate, str) or len(candidate) != 40:
        raise EvidenceError("results candidate_sha must be a full 40-character Git SHA")
    observed = results.get("cases", [])
    if not isinstance(observed, list):
        raise EvidenceError("results cases must be a list")
    by_id = {}
    for item in observed:
        case_id = item.get("id")
        if case_id in by_id:
            raise EvidenceError(f"duplicate result case id: {case_id}")
        by_id[case_id] = item
    unknown = by_id.keys() - {case["id"] for case in cases}
    if unknown:
        raise EvidenceError(
            f"results contain unmanifested cases: {', '.join(sorted(unknown))}"
        )

    evaluated = []
    for case in cases:
        result = by_id.get(case["id"])
        if result is None:
            evaluated.append(
                {
                    "id": case["id"],
                    "tier": case["tier"],
                    "system": case["system"],
                    "charge": case["charge"],
                    "multiplicity": case["multiplicity"],
                    "method": case["method"],
                    "backend": case["backend"],
                    "backend_version": None,
                    "reference": case["reference"],
                    "expected": case["acceptance"],
                    "observed": {},
                    "logs": [],
                    "status": "NOT_TESTED",
                    "failures": ["no result recorded"],
                }
            )
            continue
        status = result.get("status")
        if status not in ALLOWED_STATUS:
            raise EvidenceError(f"{case['id']}: invalid result status {status!r}")
        failures = []
        if status == "NOT_APPLICABLE":
            exception = result.get("authorized_exception", {})
            if not all(
                exception.get(key)
                for key in ("reason", "authorized_by", "approval_reference")
            ):
                status = "BLOCKED"
                failures.append("NOT_APPLICABLE lacks an authorized exception record")
            else:
                failures.append(
                    "authorized exception: "
                    + exception["reason"]
                    + f" (approved by {exception['authorized_by']}; {exception['approval_reference']})"
                )
        if status == "FAIL":
            failures.append(
                str(result.get("failure", "runner reported scientific failure"))
            )
        elif status in {"NOT_TESTED", "BLOCKED"} and result.get("failure"):
            failures.append(str(result["failure"]))
        if status == "PASS":
            for field in (
                "charge",
                "multiplicity",
                "method",
                "backend",
                "input_sha256",
            ):
                if result.get(field) != case.get(field):
                    failures.append(f"{field} differs from frozen manifest")
            if not result.get("backend_version"):
                failures.append("backend_version is missing")
            metrics = result.get("metrics", {})
            for name, rule in case["acceptance"].items():
                observed_value = metrics.get(name)
                if not isinstance(observed_value, (int, float)) or not math.isfinite(
                    observed_value
                ):
                    failures.append(f"{name}: observed finite numeric value is missing")
                    continue
                deviation = abs(observed_value - rule["expected"])
                if deviation > rule["tolerance"]:
                    failures.append(
                        f"{name}: observed {observed_value} {rule['unit']}, expected "
                        f"{rule['expected']} +/- {rule['tolerance']} {rule['unit']}"
                    )
            if not result.get("logs"):
                failures.append("raw log/artifact references are missing")
            if failures:
                status = "FAIL"
        evaluated.append(
            {
                "id": case["id"],
                "tier": case["tier"],
                "system": case["system"],
                "charge": case["charge"],
                "multiplicity": case["multiplicity"],
                "method": case["method"],
                "backend": case["backend"],
                "backend_version": result.get("backend_version"),
                "reference": case["reference"],
                "expected": case["acceptance"],
                "observed": result.get("metrics", {}),
                "logs": result.get("logs", []),
                "status": status,
                "failures": failures,
            }
        )

    counts = {
        status: sum(row["status"] == status for row in evaluated)
        for status in sorted(ALLOWED_STATUS)
    }
    verdict = (
        "PASS"
        if counts["FAIL"] == counts["BLOCKED"] == counts["NOT_TESTED"] == 0
        else "FAIL"
    )
    manifest_digest = hashlib.sha256(
        json.dumps(manifest, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    return {
        "schema_version": 1,
        "benchmark_id": manifest.get("benchmark_id", "unnamed"),
        "candidate_sha": candidate,
        "manifest_sha256": manifest_digest,
        "verdict": verdict,
        "counts": counts,
        "cases": evaluated,
    }


def render_markdown(report: dict) -> str:
    lines = [
        f"# Scientific validation: {report['benchmark_id']}",
        "",
        f"Candidate: `{report['candidate_sha']}`  ",
        f"Verdict: **{report['verdict']}**  ",
        "Counts: "
        + ", ".join(f"{key}={value}" for key, value in report["counts"].items()),
        "",
        "| Case | Tier | System / method / backend | Status | Expected vs observed | Logs | Failure details |",
        "| --- | --- | --- | --- | --- | --- | --- |",
    ]
    for case in report["cases"]:
        identity = f"{case['system']}; {case['method']}; {case['backend']} {case['backend_version'] or ''}"
        values = json.dumps(
            {"expected": case["expected"], "observed": case["observed"]}, sort_keys=True
        )
        lines.append(
            f"| {case['id']} | {case['tier']} | {identity} | {case['status']} | "
            f"`{values}` | {', '.join(case['logs'])} | {'; '.join(case['failures'])} |"
        )
    return "\n".join(lines) + "\n"


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("results", type=Path)
    parser.add_argument("--json", dest="json_output", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
        results = json.loads(args.results.read_text(encoding="utf-8"))
        report = evaluate(manifest, results, args.manifest.parent)
    except (OSError, json.JSONDecodeError, EvidenceError) as exc:
        parser.error(str(exc))
    args.json_output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    args.markdown.write_text(render_markdown(report), encoding="utf-8")
    return 0 if report["verdict"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
