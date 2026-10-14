"""Filesystem integrity and process-boundary contracts, without calculators."""
import json
import os
from pathlib import Path
import subprocess
import sys
from uuid import uuid4

import pytest

from pyar.persistence import (ArtifactRef, ArtifactError, RunLayout, RunLock, LockError,
                              detect_version, LayoutError, write_json, read_json,
                              canonical_json_bytes)


def authorities(root, workflow="grow", run_id=None):
    run_id = run_id or str(uuid4())
    write_json(root / "layout.json", {"format": "pyar-run-layout", "layout_version": 2,
                                      "workflow": workflow, "run_id": run_id})
    if workflow != "conformer":
        write_json(root / "request.json", {"request_schema_version": 1,
                                           "workflow": workflow, "run_id": run_id})
    if workflow != "scan-bond":
        write_json(root / "state.json", {"workflow_state_schema_version": {"grow": 2, "aggregate": 2, "microsolvation": 2,
                                                                            "reaction": 3, "conformer": 4}[workflow],
                                         "workflow": workflow, "run_id": run_id})
    return run_id


def test_paths_and_detection_are_read_only(tmp_path):
    root = tmp_path / "grow"
    layout = RunLayout(root, "grow")
    assert layout.step(1000) == root / "steps/step_1000"
    assert layout.results == root / "results"
    assert detect_version(root, "grow") == 0
    assert not root.exists()
    with pytest.raises(Exception):
        layout.root = tmp_path
    for index in (-1, True, 1.5):
        with pytest.raises(ValueError):
            layout.step(index)


@pytest.mark.parametrize("path", ["", "/etc/passwd", "../outside", "a/../b", "a//b", "a/./b",
                                  "C:/data", "a\\b", "./x", "a/", "\x00"])
def test_unsafe_artifact_paths_rejected(tmp_path, path):
    with pytest.raises(ValueError):
        RunLayout(tmp_path, "grow").path(path)
    with pytest.raises(ArtifactError):
        ArtifactRef(str(uuid4()), "result", path, "0" * 64)


def test_artifacts_keep_identity_verify_bytes_and_reject_escape(tmp_path):
    root = tmp_path / "run"
    root.mkdir()
    (root / "structure.xyz").write_bytes(b"published geometry")
    ref = ArtifactRef.capture(root, "structure.xyz", "selected")
    loaded = ArtifactRef.from_dict(ref.to_dict())
    assert loaded == ref
    assert loaded.verify(root) == root / "structure.xyz"
    (root / "structure.xyz").write_bytes(b"changed geometry")
    with pytest.raises(ArtifactError, match="checksum"):
        loaded.verify(root)
    (root / "structure.xyz").unlink()
    with pytest.raises(ArtifactError, match="Cannot read"):
        loaded.verify(root)
    (root / "escape").symlink_to(tmp_path, target_is_directory=True)
    with pytest.raises(ArtifactError, match="escapes"):
        ArtifactRef.capture(root, "escape/outside", "result")
    with pytest.raises(ValueError, match="escapes"):
        RunLayout(root, "grow").path("escape/outside")


def test_json_determinism_and_failed_publication_preserve_state(tmp_path, monkeypatch):
    assert canonical_json_bytes({"z": "é", "a": 1}) == b'{"a":1,"z":"\xc3\xa9"}\n'
    path = tmp_path / "state.json"
    write_json(path, {"old": True})
    before = path.read_bytes()
    for value in ({"energy": float("nan")}, {1: "coerced"}):
        with pytest.raises(ValueError):
            write_json(path, value)
        assert path.read_bytes() == before
    with pytest.raises(FileExistsError):
        write_json(path, {"new": True}, overwrite=False)
    assert path.read_bytes() == before
    def fail(*args):
        raise OSError("simulated fsync crash")
    monkeypatch.setattr(os, "fsync", fail)
    with pytest.raises(OSError):
        write_json(path, {"new": True})
    assert path.read_bytes() == before
    assert not list(tmp_path.glob(".pyar-json-*"))


def test_strict_json_and_immutable_publication(tmp_path):
    path = tmp_path / "commit.json"
    write_json(path, {"stage": "complete"}, overwrite=False)
    assert read_json(path) == {"stage": "complete"}
    for raw in ('{"x":1,"x":2}', '{"x":NaN}', '{"x":Infinity}', '{"x":1e999}'):
        path.write_text(raw)
        with pytest.raises(ValueError):
            read_json(path)
    link = tmp_path / "symlink.json"
    link.symlink_to(path)
    with pytest.raises(ValueError, match="symlink"):
        write_json(link, {})


@pytest.mark.parametrize("workflow", ["grow", "aggregate", "microsolvation", "reaction", "scan-bond", "conformer"])
def test_version_detection_requires_explicit_identity(tmp_path, workflow):
    authorities(tmp_path, workflow)
    before = {p.name: p.read_bytes() for p in tmp_path.iterdir()}
    assert detect_version(tmp_path, workflow) == 2
    assert before == {p.name: p.read_bytes() for p in tmp_path.iterdir()}
    marker = read_json(tmp_path / "layout.json")
    for change in ({"layout_version": 99}, {"layout_version": True}, {"run_id": "wrong"},
                   {"workflow": "not-a-workflow"}):
        write_json(tmp_path / "layout.json", dict(marker, **change))
        with pytest.raises(LayoutError):
            detect_version(tmp_path, workflow)
    write_json(tmp_path / "layout.json", marker)
    authority = tmp_path / ("state.json" if workflow == "conformer" else "request.json")
    write_json(authority, {"version": 1})
    with pytest.raises(LayoutError, match="authority"):
        detect_version(tmp_path, workflow)


def test_legacy_detection_delegates_and_never_migrates(tmp_path):
    path = tmp_path / "request.json"
    path.write_text('{"request_signature":"legacy scan"}')
    before = path.read_bytes()
    with pytest.raises(LayoutError, match="unrecognized"):
        detect_version(tmp_path, "scan-bond")
    assert detect_version(tmp_path, "scan-bond", validate_v1=lambda root:
                          read_json(root / "request.json")["request_signature"] == "legacy scan") == 1
    assert path.read_bytes() == before and not (tmp_path / "layout.json").exists()
    path.unlink()
    write_json(tmp_path / "layout.json", {"format": "pyar-run-layout", "layout_version": 2,
                                         "workflow": "grow", "run_id": str(uuid4())})
    with pytest.raises(LayoutError, match="incomplete"):
        detect_version(tmp_path, "grow")


def child(root, code):
    return subprocess.run([sys.executable, "-c", "from pyar.persistence import RunLock,LockError; " + code,
                           str(root)], capture_output=True, text=True, timeout=15)


@pytest.mark.skipif(os.name != "posix", reason="local POSIX flock contract")
def test_lock_excludes_another_process_and_preserves_inode(tmp_path):
    root = tmp_path / "fresh"
    guard = RunLock(root)
    with guard:
        inode = guard.path.stat().st_ino
        blocked = child(root, "import sys; RunLock(sys.argv[1],recover_stale=True).acquire()")
        assert blocked.returncode != 0 and "active writer" in blocked.stderr
        assert not root.exists()
    with RunLock(root):
        assert guard.path.stat().st_ino == inode
    assert read_json(guard.path)["active"] is False


@pytest.mark.skipif(os.name != "posix", reason="local POSIX flock contract")
def test_process_death_requires_explicit_recovery(tmp_path):
    root = tmp_path / "fresh"
    crashed = child(root, "import sys,os; g=RunLock(sys.argv[1]).acquire(); os._exit(0)")
    assert crashed.returncode == 0, crashed.stderr
    path = RunLock(root).path
    before = path.read_bytes()
    with pytest.raises(LockError, match="Stale"):
        RunLock(root).acquire()
    assert path.read_bytes() == before
    with RunLock(root, recover_stale=True):
        assert "recovery" in read_json(path)
    with RunLock(root):
        pass


def test_corrupt_and_symlink_lock_fail_closed(tmp_path):
    root = tmp_path / "run"
    guard = RunLock(root)
    guard.path.write_bytes(b"interrupted metadata")
    with pytest.raises(LockError, match="corrupt"):
        guard.acquire()
    with RunLock(root, recover_stale=True):
        assert read_json(guard.path)["recovery"]["previous_metadata_sha256"]
    guard.path.unlink()
    target = tmp_path / "target"
    target.write_text("untouched")
    guard.path.symlink_to(target)
    with pytest.raises(LockError):
        guard.acquire()
    assert target.read_text() == "untouched"


def test_invalid_serialization_leaves_no_directories(tmp_path):
    with pytest.raises(ValueError):
        write_json(tmp_path / "missing/commit.json", {"bad": float("inf")})
    assert not (tmp_path / "missing").exists()


def test_publication_failures_on_both_sides_of_replace(tmp_path, monkeypatch):
    from pyar.persistence import json_store
    path = tmp_path / "state.json"
    write_json(path, {"old": 1})
    def fail(*args):
        raise OSError("publication interrupted")
    with monkeypatch.context() as patcher:
        patcher.setattr(os, "replace", fail)
        with pytest.raises(OSError):
            write_json(path, {"new": 2})
    assert read_json(path) == {"old": 1}
    with monkeypatch.context() as patcher:
        patcher.setattr(json_store, "fsync_directory", fail)
        with pytest.raises(OSError):
            write_json(path, {"new": 2})
    # A post-publication durability error is surfaced, not rolled back silently.
    assert read_json(path) == {"new": 2}
    assert not list(tmp_path.glob(".pyar-json-*"))


def test_unknown_authority_schema_and_symlinks_rejected(tmp_path):
    authorities(tmp_path)
    state = read_json(tmp_path / "state.json")
    write_json(tmp_path / "state.json", dict(state, workflow_state_schema_version=99))
    with pytest.raises(LayoutError):
        detect_version(tmp_path, "grow")
    (tmp_path / "state.json").unlink()
    (tmp_path / "state.json").symlink_to(tmp_path / "request.json")
    with pytest.raises(LayoutError):
        detect_version(tmp_path, "grow")


@pytest.mark.skipif(os.name != "posix", reason="local POSIX flock contract")
def test_lock_rejects_special_files_and_bad_recovery_option(tmp_path):
    with pytest.raises(LockError):
        RunLock(tmp_path / "run", recover_stale="yes")
    lock = RunLock(tmp_path / "run")
    os.mkfifo(lock.path)
    with pytest.raises(LockError, match="regular file"):
        lock.acquire()


def test_immutable_record_race_has_exactly_one_winner(tmp_path):
    from concurrent.futures import ThreadPoolExecutor
    from threading import Barrier
    barrier = Barrier(2)
    path = tmp_path / "commit.json"
    def publish(value):
        barrier.wait(timeout=5)
        try:
            write_json(path, {"winner": value}, overwrite=False)
            return value
        except FileExistsError:
            return None
    with ThreadPoolExecutor(max_workers=2) as pool:
        results = list(pool.map(publish, [1, 2]))
    winner, = [item for item in results if item is not None]
    assert read_json(path) == {"winner": winner}
    assert not list(tmp_path.glob(".pyar-json-*"))


def test_new_directory_entries_are_flushed(tmp_path, monkeypatch):
    from pyar.persistence import json_store
    flushed = []
    original = json_store.fsync_directory
    def record(path):
        flushed.append(Path(path))
        original(path)
    monkeypatch.setattr(json_store, "fsync_directory", record)
    write_json(tmp_path / "new/nested/state.json", {"complete": True})
    assert flushed == [tmp_path, tmp_path / "new", tmp_path / "new/nested"]


@pytest.mark.skipif(os.name != "posix", reason="local POSIX flock contract")
def test_exception_releases_writer_and_keeps_recovery_evidence(tmp_path):
    guard = RunLock(tmp_path / "run")
    guard.path.write_text("crashed metadata")
    with pytest.raises(RuntimeError, match="workflow failure"):
        with RunLock(guard.root, recover_stale=True):
            raise RuntimeError("workflow failure")
    recovered = read_json(guard.path)["recovery"]
    with RunLock(guard.root):
        pass
    assert read_json(guard.path)["recovery"] == recovered
