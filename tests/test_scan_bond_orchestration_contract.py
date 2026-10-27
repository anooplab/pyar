"""Observable scan-bond outcomes across preparation and stage boundaries."""

import importlib
import json

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
