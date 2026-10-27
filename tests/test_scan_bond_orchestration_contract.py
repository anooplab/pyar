"""Observable scan-bond outcomes across preparation and stage boundaries."""

import importlib
import csv
import hashlib
import io
import json
from pathlib import Path

import numpy as np
import pytest

from pyar.backends import write_xyz
from pyar.backends.orca_scan import OrcaBondScanResult
from pyar.core.molecule import Molecule


workflow = importlib.import_module("pyar.workflows.scan_bond")


@pytest.fixture
def fragment_inputs(tmp_path):
    left, right = tmp_path / "a.xyz", tmp_path / "b.xyz"
    write_xyz(["H"], [[0.0, 0.0, 0.0]], left)
    write_xyz(["H"], [[1.2, 0.0, 0.0]], right)
    return left, right


def _orientation():
    molecule = Molecule(
        ["H", "H"], [[0.0, 0.0, 0.0], [1.2, 0.0, 0.0]], name="orientation_000"
    )
    molecule.fragments = [[0], [1]]
    return molecule


def _settings():
    return {"software": "orca", "method": "BP86", "basis": "def2-SVP"}


@pytest.mark.parametrize(
    "settings,orientations,scan_step,scan_points,message",
    [
        ({"software": "invalid"}, 0, -1, 2, "unsupported: invalid"),
        (_settings(), 0, -1, 2, "orientation count must be a positive integer"),
        (_settings(), 1, -1, 2, "scan step must be finite and positive"),
        (_settings(), 1, 0.1, 2, "cannot be used together"),
    ],
)
def test_rejected_route_keeps_first_error_and_does_not_allocate(
    tmp_path, fragment_inputs, settings, orientations, scan_step, scan_points, message
):
    output = tmp_path / "scan"
    with pytest.raises(ValueError, match=message):
        workflow.run_scan_bond(
            *fragment_inputs,
            (0, 0),
            orientations,
            settings,
            output,
            scan_step=scan_step,
            scan_points=scan_points,
            through="scan",
        )
    assert not output.exists()


def test_failed_placement_remains_one_failed_orientation_without_backend(
    tmp_path, fragment_inputs, monkeypatch
):
    monkeypatch.setattr(
        workflow, "_make_orientations", lambda *args: [(0, None, "placement failed")]
    )
    monkeypatch.setattr(
        workflow, "separated_reactant_identity", lambda *args: {"inchi": "reactants"}
    )
    monkeypatch.setattr(
        workflow,
        "run_orca_bond_scan",
        lambda *args: pytest.fail("placement failure started a backend"),
    )
    settings = _settings()
    output = tmp_path / "scan"
    response = workflow.run_scan_bond(
        *fragment_inputs, (0, 0), 1, settings, output, through="scan"
    )

    assert settings == _settings()
    assert response["status"] == "failed"
    assert response["scan_status"] == "scan_incomplete"
    assert response["results"] == json.loads((output / "summary.json").read_text())
    assert response["results"][0]["status"] == "scan_failed"
    assert response["results"][0]["scan_status"] == "orientation_generation_failed"
    assert response["results"][0]["error"] == "placement failed"
    assert (
        json.loads((output / "orientation_000/result.json").read_text())
        == response["results"][0]
    )
    assert not (output / "orientation_000/scan_commit.json").exists()


def test_profile_failure_preserves_scan_evidence_but_never_commits(
    tmp_path, fragment_inputs, monkeypatch
):
    monkeypatch.setattr(
        workflow, "_make_orientations", lambda *args: [(0, _orientation(), None)]
    )
    monkeypatch.setattr(
        workflow, "separated_reactant_identity", lambda *args: {"inchi": "reactants"}
    )

    def scan(molecule, request, directory, settings):
        directory.mkdir(parents=True, exist_ok=True)
        distances = [request.start_distance_angstrom, request.end_distance_angstrom]
        frames = []
        for distance in distances:
            frame = np.array([[0.0, 0.0, 0.0], [distance, 0.0, 0.0]])
            frames.append((molecule.atoms_list, frame, "frame"))
        profile = [
            {"target_distance_angstrom": distance, "energy_hartree": energy}
            for distance, energy in zip(distances, [-1.0, -0.9])
        ]
        paths = [
            directory / name
            for name in (
                "scan_trajectory.xyz",
                "final_scan.xyz",
                "scan.inp",
                "scan.out",
            )
        ]
        for path in paths:
            path.write_text("native evidence\n")
        return OrcaBondScanResult(
            True,
            paths[0],
            paths[1],
            frames[-1][1],
            paths[2],
            paths[3],
            "success",
            frames,
            profile,
        )

    monkeypatch.setattr(workflow, "run_orca_bond_scan", scan)
    monkeypatch.setattr(
        workflow,
        "_write_scan_profile",
        lambda *args, **kwargs: (_ for _ in ()).throw(OSError("profile write failed")),
    )
    output = tmp_path / "scan"
    response = workflow.run_scan_bond(
        *fragment_inputs,
        (0, 0),
        1,
        _settings(),
        output,
        scan_end=0.6,
        scan_points=2,
        through="scan",
    )
    result = response["results"][0]

    assert result["status"] == "scan_success_analysis_failed"
    assert result["scan_status"] == "success"
    assert result["error"] == "profile write failed"
    assert response["status"] == "failed"
    assert not (output / "orientation_000/scan_commit.json").exists()
    assert (
        output / "orientation_000/attempts/attempt_000/execution/scan/scan.out"
    ).read_text() == "native evidence\n"


def test_completed_artifact_hashes_stream_without_whole_file_reads(tmp_path, monkeypatch):
    paths = workflow._orientation_paths(tmp_path, 0, None)
    paths.directory.mkdir()
    for path in paths.completed_artifacts:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"evidence" if path.name != "scan_trajectory.xyz" else b"x" * 3_000_000)
    expected = {
        str(path.relative_to(paths.directory)): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in paths.completed_artifacts
    }

    def reject_whole_file_read(self):
        raise AssertionError(f"whole-file read attempted for {self}")

    monkeypatch.setattr(Path, "read_bytes", reject_whole_file_read)
    result = {"orientation": 0, "status": "scan_success_relax_success"}
    workflow._save_orientation_result(
        paths.directory, result,
        qc_params={"total_multiplicity": 1, "total_multiplicity_source": "explicit"},
        intramolecular=False, completed_artifacts=paths.completed_artifacts, layout=None,
    )
    assert result["completed_artifact_hashes"] == expected
    recovered = workflow._recover_orientation(paths, 0, {}, None)
    assert recovered.complete
    assert recovered.previous_result["completed_artifact_hashes"] == expected


def test_interrupted_derived_report_keeps_previous_complete_bytes(tmp_path, monkeypatch):
    target = tmp_path / "result.json"
    previous = b'{\n  "orientation": 0,\n  "status": "scan_failed"\n}'
    target.write_bytes(previous)
    real_open = Path.open

    class InterruptedTemporary:
        def __init__(self, stream):
            self.stream = stream
            self.name = stream.name

        def __enter__(self):
            self.stream.__enter__()
            return self

        def __exit__(self, *args):
            return self.stream.__exit__(*args)

        def fileno(self):
            return self.stream.fileno()

        def write(self, payload):
            self.stream.write(payload[:12])
            raise KeyboardInterrupt("injected partial temporary write")

    def interrupt_temporary_open(path, *args, **kwargs):
        stream = real_open(path, *args, **kwargs)
        return InterruptedTemporary(stream) if path.name.startswith(".pyar-report-") else stream

    monkeypatch.setattr(Path, "open", interrupt_temporary_open)
    with pytest.raises(KeyboardInterrupt, match="partial temporary write"):
        workflow._publish_derived_report(target, '{"new": true}')
    assert target.read_bytes() == previous
    assert json.loads(target.read_text())["status"] == "scan_failed"
    assert not list(tmp_path.glob(".pyar-report-*"))


def test_derived_report_fsync_failure_after_replace_remains_parseable(tmp_path, monkeypatch):
    from pyar.persistence import json_store

    target = tmp_path / "result.json"
    target.write_text('{"old": true}')
    real_fsync = json_store.os.fsync
    calls = 0

    def fail_after_replace(descriptor):
        nonlocal calls
        calls += 1
        if calls == 2:
            raise OSError("injected directory fsync failure")
        return real_fsync(descriptor)

    monkeypatch.setattr(json_store.os, "fsync", fail_after_replace)
    with pytest.raises(OSError, match="directory fsync failure"):
        workflow._publish_derived_report(target, '{"new": true}')
    assert target.read_bytes() == b'{"new": true}'
    assert json.loads(target.read_text()) == {"new": True}
    assert not list(tmp_path.glob(".pyar-report-*"))


def test_derived_report_bytes_match_previous_text_and_csv_writers(tmp_path):
    value = {"name": "café", "status": "complete"}
    json_text = json.dumps(value, indent=2, sort_keys=True, default=str)
    old_json, new_json = tmp_path / "old.json", tmp_path / "new.json"
    old_json.write_text(json_text)
    workflow._publish_derived_report(new_json, json_text)
    assert new_json.read_bytes() == old_json.read_bytes()
    assert new_json.stat().st_mode & 0o777 == old_json.stat().st_mode & 0o777
    new_json.chmod(0o640)
    workflow._publish_derived_report(new_json, json_text)
    assert new_json.stat().st_mode & 0o777 == 0o640

    old_csv, new_csv = tmp_path / "old.csv", tmp_path / "new.csv"
    with old_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=sorted(value))
        writer.writeheader()
        writer.writerow(value)
    with io.StringIO(newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=sorted(value))
        writer.writeheader()
        writer.writerow(value)
        workflow._publish_derived_report(new_csv, handle.getvalue())
    assert new_csv.read_bytes() == old_csv.read_bytes()


def test_derived_report_rejects_symlink_without_touching_target(tmp_path):
    outside = tmp_path / "outside.json"
    outside.write_text('{"keep": true}')
    report = tmp_path / "result.json"
    report.symlink_to(outside)
    with pytest.raises(ValueError, match="Refusing report publication through a symlink"):
        workflow._publish_derived_report(report, '{"new": true}')
    assert outside.read_text() == '{"keep": true}'
