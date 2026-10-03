import json
import random
import tempfile
from pathlib import Path
from types import SimpleNamespace

import pytest

from benchmarks.ts_optimizers.prepare import run_pilot


def test_pilot_runner_shuffles_pairs_and_uses_fresh_sella_cache(tmp_path, monkeypatch):
    manifest = tmp_path / "pilot.json"
    manifest.write_text(json.dumps({"cases": [{"id": "r1_easy"}, {"id": "r2_hard"}]}))
    output = tmp_path / "run"
    calls = []

    def fake_run(command, *, cwd, env, capture_output, text, check):
        case_id = command[command.index("--case") + 1]
        optimizer = command[command.index("--optimizer") + 1]
        output_dir = Path(command[command.index("--output") + 1])
        calls.append((case_id, optimizer, env.get("JAX_COMPILATION_CACHE_DIR")))
        result = output_dir / "cases" / case_id / optimizer / "benchmark_result.json"
        result.parent.mkdir(parents=True, exist_ok=True)
        result.write_text(json.dumps({"case_id": case_id, "optimizer": optimizer}))
        return SimpleNamespace(returncode=0, stdout=json.dumps({"outcome": "complete"}), stderr="")

    monkeypatch.setattr(run_pilot.subprocess, "run", fake_run)
    run_pilot.run_pilot(manifest, output, schedule_seed=33)
    execution = json.loads((output / "pilot_execution.json").read_text())
    expected = [
        {"case_id": case_id, "optimizer": optimizer}
        for case_id in ("r1_easy", "r2_hard") for optimizer in ("geometric", "sella")
    ]
    random.Random(33).shuffle(expected)
    assert execution["randomized_jobs"] == expected
    assert len(calls) == 4
    sella_cache_dirs = [
        Path(cache_dir) for _, optimizer, cache_dir in calls if optimizer == "sella"
    ]
    assert len(sella_cache_dirs) == 2
    assert all(not cache_dir.exists() for cache_dir in sella_cache_dirs)
    assert all(cache_dir.parent == Path(tempfile.gettempdir()) for cache_dir in sella_cache_dirs)

    run_pilot.run_pilot(manifest, output, schedule_seed=33)
    assert len(calls) == 4
    with pytest.raises(ValueError, match="different benchmark or schedule"):
        run_pilot.run_pilot(manifest, output, schedule_seed=34)


def test_pilot_runner_rejects_mismatched_existing_run_manifest(tmp_path):
    manifest = tmp_path / "pilot.json"
    manifest.write_text(json.dumps({"cases": [{"id": "r1_easy"}]}))
    output = tmp_path / "run"
    output.mkdir()
    (output / "run_manifest.json").write_text(json.dumps({"manifest_sha256": "wrong"}))
    with pytest.raises(ValueError, match="different benchmark manifest"):
        run_pilot.run_pilot(manifest, output)
