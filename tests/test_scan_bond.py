import json
import importlib
from pathlib import Path

import numpy as np
import pytest

from pyar.backends import write_xyz
from pyar.backends.orca_scan import OrcaBondScanRequest, OrcaBondScanResult, parse_xyz_trajectory, run_orca_bond_scan
from pyar.core.molecule import Molecule
from pyar.workflows.scan_bond import absolute_target_indices, run_scan_bond, scan_point_count


def test_fragment_local_indices_are_merged_asymmetric():
    fragment_a = Molecule(["C"] * 5, np.zeros((5, 3)))
    assert absolute_target_indices(fragment_a, 3, 1) == (3, 6)
    assert absolute_target_indices(fragment_a, 0, 0) == (0, 5)


@pytest.mark.parametrize("bad", [-1, 5])
def test_fragment_a_index_is_validated(bad):
    fragment_a = Molecule(["C"] * 5, np.zeros((5, 3)))
    with pytest.raises(ValueError, match="fragment A"):
        absolute_target_indices(fragment_a, bad, 0)


@pytest.mark.parametrize("start,end,step,expected", [
    (2.0, 1.0, 0.1, 11),
    (1.0, 2.0, 0.3, 5),
    (1.0, 1.01, 0.1, 2),
])
def test_scan_point_count_includes_endpoints(start, end, step, expected):
    assert scan_point_count(start, end, step=step) == expected


def test_scan_point_count_validates_exclusive_controls():
    with pytest.raises(ValueError, match="cannot"):
        scan_point_count(2.0, 1.0, step=0.1, points=4)
    with pytest.raises(ValueError, match="differ"):
        scan_point_count(1.0, 1.0, step=0.1)
    with pytest.raises(ValueError, match="positive"):
        scan_point_count(2.0, 1.0, step=0.0)


def test_multixyz_parser_returns_final_frame_and_rejects_truncation(tmp_path):
    trajectory = tmp_path / "scan.allxyz"
    trajectory.write_text(
        "2\nframe 1\nH 0 0 0\nH 1 0 0\n"
        "2\nframe 2\nH 0 0 0\nH 1.2 0 0\n"
    )
    frames = parse_xyz_trajectory(trajectory, 2, ["H", "H"])
    assert len(frames) == 2
    np.testing.assert_allclose(frames[-1][1][1], [1.2, 0.0, 0.0])
    trajectory.write_text("2\ntruncated\nH 0 0 0\n")
    with pytest.raises(ValueError, match="Truncated"):
        parse_xyz_trajectory(trajectory, 2, ["H", "H"])


def test_orca_scan_input_and_recovery(tmp_path, monkeypatch):
    molecule = Molecule(["H", "H"], np.asarray([[0., 0., 0.], [2., 0., 0.]]))
    molecule.scftype = "uks"
    request = OrcaBondScanRequest(0, 1, 2.0, 1.0, 11)
    scan_dir = tmp_path / "scan"

    def fake_run(command, stdout_path=None, stderr_path=None):
        assert Path.cwd() == scan_dir
        Path(stdout_path).write_text("****ORCA TERMINATED NORMALLY****\n")
        Path("scan.allxyz").write_text(
            "2\nscan\nH 0 0 0\nH 2 0 0\n"
            "2\nscan\nH 0 0 0\nH 1 0 0\n"
        )
        return 0

    monkeypatch.setattr("pyar.backends.orca_scan.require_executable", lambda *args: "orca")
    monkeypatch.setattr("pyar.backends.orca_scan.run_command", fake_run)
    result = run_orca_bond_scan(
        molecule, request, scan_dir,
        {"method": "BP86", "basis": "def2-SVP", "nprocs": 1, "scf_cycles": 1000,
         "opt_cycles": 25, "opt_threshold": "tight"},
    )
    assert result.success
    assert "B 0 1 = 2.0000000000, 1.0000000000, 11" in result.input_path.read_text()
    assert "UKS" in result.input_path.read_text()
    assert "MaxIter 25" in result.input_path.read_text()
    assert result.final_geometry_path.exists()
    assert result.trajectory_path.exists()
    assert result.final_coordinates[1, 0] == 1.0


def test_workflow_relaxes_scan_frame_and_writes_summary(tmp_path, monkeypatch):
    a = tmp_path / "a.xyz"
    b = tmp_path / "b.xyz"
    write_xyz(["H"], [[0., 0., 0.]], a)
    write_xyz(["H"], [[1.8, 0., 0.]], b)
    scan_frame = np.asarray([[0., 0., 0.], [1.0, 0., 0.]])

    def fake_scan(molecule, request, directory, qc_params):
        directory.mkdir(parents=True, exist_ok=True)
        input_path = directory / "scan.inp"
        output_path = directory / "scan.out"
        trajectory = directory / "scan_trajectory.xyz"
        final_path = directory / "final_scan.xyz"
        input_path.write_text("scan")
        output_path.write_text("normal")
        write_xyz(molecule.atoms_list, scan_frame, final_path)
        trajectory.write_text("2\nframe\nH 0 0 0\nH 1 0 0\n")
        return OrcaBondScanResult(True, trajectory, final_path, scan_frame,
                                  input_path, output_path, "success")

    def fake_optimise(molecule, qc_params):
        molecule.coordinates = np.asarray([[0., 0., 0.], [0.7, 0., 0.]])
        Path("job_relaxed").mkdir()
        write_xyz(molecule.atoms_list, molecule.coordinates,
                  Path("job_relaxed") / "result_relaxed.xyz")
        return True

    scan_workflow = importlib.import_module("pyar.workflows.scan_bond")
    monkeypatch.setattr(scan_workflow, "run_orca_bond_scan", fake_scan)
    monkeypatch.setattr(scan_workflow, "optimise", fake_optimise)
    result = run_scan_bond(
        a, b, (0, 0), 1,
        {"software": "orca", "method": "BP86", "basis": "def2-SVP",
         "nprocs": 1, "scf_cycles": 1000},
        tmp_path / "scan_bond", scan_end=0.5, scan_points=2,
    )
    assert result["results"][0]["status"] == "scan_success_relax_success", result
    assert result["results"][0]["target_bond_present_after_relaxation"]
    assert (tmp_path / "scan_bond" / "summary.json").exists()
    assert (tmp_path / "scan_bond" / "orientation_000" / "result_relaxed.xyz").exists()
    assert (tmp_path / "scan_bond" / "orientation_000" / "result.json").exists()
