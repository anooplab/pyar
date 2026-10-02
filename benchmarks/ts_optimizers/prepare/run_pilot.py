"""Run a shuffled, paired TS-optimizer pilot without overwriting job data."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import random
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[3]


def _sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def argument_parse(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", help="prepared pyar-benchmark-ts JSON manifest")
    parser.add_argument("--output", required=True, help="benchmark run directory")
    parser.add_argument("--schedule-seed", type=int, default=20261003)
    return parser.parse_args(argv)


def run_pilot(manifest_path, output, *, schedule_seed=20261003):
    manifest_path = Path(manifest_path).expanduser().resolve()
    output = Path(output).expanduser().resolve()
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    cases = manifest.get("cases")
    if not isinstance(cases, list) or not cases:
        raise ValueError("manifest has no cases")
    tasks = [
        {"case_id": case["id"], "optimizer": optimizer}
        for case in cases for optimizer in ("geometric", "sella")
    ]
    random.Random(schedule_seed).shuffle(tasks)
    output.mkdir(parents=True, exist_ok=True)
    schedule_path = output / "pilot_execution.json"
    plan = {
        "schema_version": 1,
        "manifest_sha256": _sha256(manifest_path),
        "schedule_seed": schedule_seed,
        "randomized_jobs": tasks,
        "jax_cache_policy": "fresh temporary JAX compilation cache for each Sella job",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    run_manifest_path = output / "run_manifest.json"
    if run_manifest_path.is_file():
        run_manifest = json.loads(run_manifest_path.read_text(encoding="utf-8"))
        if run_manifest.get("manifest_sha256") != plan["manifest_sha256"]:
            raise ValueError("run_manifest.json belongs to a different benchmark manifest")
    if schedule_path.exists():
        previous = json.loads(schedule_path.read_text(encoding="utf-8"))
        for key in ("manifest_sha256", "schedule_seed", "randomized_jobs", "jax_cache_policy"):
            if previous.get(key) != plan[key]:
                raise ValueError(
                    "pilot_execution.json belongs to a different benchmark or schedule"
                )
        plan = previous
    else:
        temporary = schedule_path.with_suffix(".tmp")
        temporary.write_text(json.dumps(plan, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        os.replace(temporary, schedule_path)

    environment = os.environ.copy()
    current_pythonpath = environment.get("PYTHONPATH")
    environment["PYTHONPATH"] = (
        str(REPOSITORY_ROOT) if not current_pythonpath
        else str(REPOSITORY_ROOT) + os.pathsep + current_pythonpath
    )
    for index, task in enumerate(plan["randomized_jobs"], start=1):
        result_path = (
            output / "cases" / task["case_id"] / task["optimizer"] / "benchmark_result.json"
        )
        if result_path.is_file():
            print(
                f"[{index}/{len(tasks)}] already recorded: "
                f"{task['case_id']} {task['optimizer']}", flush=True,
            )
            continue
        command = [
            sys.executable, "-m", "pyar.scripts.ts_optimizer_benchmark", "run",
            str(manifest_path), "--case", task["case_id"], "--optimizer", task["optimizer"],
            "--output", str(output),
        ]
        if task["optimizer"] == "sella":
            with tempfile.TemporaryDirectory(prefix="pyar-sella-jax-") as cache_dir:
                job_environment = dict(environment, JAX_COMPILATION_CACHE_DIR=cache_dir)
                completed = subprocess.run(
                    command, cwd=REPOSITORY_ROOT, env=job_environment,
                    capture_output=True, text=True, check=False,
                )
        else:
            completed = subprocess.run(
                command, cwd=REPOSITORY_ROOT, env=environment,
                capture_output=True, text=True, check=False,
            )
        print(f"[{index}/{len(tasks)}] {task['case_id']} {task['optimizer']} ", end="", flush=True)
        if completed.returncode:
            print(
                "runner error: "
                f"{completed.stderr.strip() or completed.stdout.strip()}", flush=True,
            )
            continue
        try:
            result = json.loads(completed.stdout)
            print(f"{result['outcome']}", flush=True)
        except (KeyError, json.JSONDecodeError):
            print("finished without a parseable result line", flush=True)


def main(argv=None):
    args = argument_parse(argv)
    run_pilot(args.manifest, args.output, schedule_seed=args.schedule_seed)
    return None


if __name__ == "__main__":
    main()
