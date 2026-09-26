import json
import importlib
from pathlib import Path
from unittest.mock import Mock

import numpy as np
import pytest

from pyar.backends import write_xyz
from pyar.backends.orca_scan import (
    OrcaBondScanRequest,
    OrcaBondScanResult,
    _orca_keyword,
    parse_orca_scan_profile,
    parse_xyz_trajectory,
    load_orca_bond_scan_result,
    run_orca_bond_scan,
)
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
        "2\nframe 1\nH 0 0 0\nH 1 0 0\n>\n"
        "2\nframe 2\nH 0 0 0\nH 1.2 0 0\n"
    )
    frames = parse_xyz_trajectory(trajectory, 2, ["H", "H"])
    assert len(frames) == 2
    np.testing.assert_allclose(frames[-1][1][1], [1.2, 0.0, 0.0])
    trajectory.write_text("2\ntruncated\nH 0 0 0\n")
    with pytest.raises(ValueError, match="Truncated"):
        parse_xyz_trajectory(trajectory, 2, ["H", "H"])


def test_multixyz_parser_rejects_nonfinite_coordinates(tmp_path):
    trajectory = tmp_path / "scan.allxyz"
    trajectory.write_text("2\nframe\nH 0 0 0\nH nan 0 0\n")
    with pytest.raises(ValueError, match="Non-finite"):
        parse_xyz_trajectory(trajectory, 2, ["H", "H"])


def test_orca_scan_input_and_recovery(tmp_path, monkeypatch):
    molecule = Molecule(["H", "H"], np.asarray([[0., 0., 0.], [2., 0., 0.]]))
    molecule.scftype = "uks"
    request = OrcaBondScanRequest(0, 1, 2.0, 1.0, 2)
    scan_dir = tmp_path / "scan"

    def fake_run(command, stdout_path=None, stderr_path=None):
        assert Path.cwd() == scan_dir
        Path(stdout_path).write_text("****ORCA TERMINATED NORMALLY****\n")
        Path("scan.allxyz").write_text(
            "2\nscan\nH 0 0 0\nH 2 0 0\n"
            "2\nscan\nH 0 0 0\nH 1 0 0\n"
        )
        Path("scan.relaxscanact.dat").write_text("2.0 -1.0\n1.0 -0.5\n")
        return 0

    monkeypatch.setattr("pyar.backends.orca_scan.require_executable", lambda *args: "orca")
    monkeypatch.setattr("pyar.backends.orca_scan.run_command", fake_run)
    result = run_orca_bond_scan(
        molecule, request, scan_dir,
        {"method": "BP86", "basis": "def2-SVP", "nprocs": 1, "scf_cycles": 1000,
         "opt_cycles": 25, "opt_threshold": "tight"},
    )
    assert result.success
    assert "B 0 1 = 2.0000000000, 1.0000000000, 2" in result.input_path.read_text()
    assert "UKS" in result.input_path.read_text()
    assert "MaxIter 25" in result.input_path.read_text()
    assert result.final_geometry_path.exists()
    assert result.trajectory_path.exists()
    assert result.final_coordinates[1, 0] == 1.0
    assert result.profile == [
        {"target_distance_angstrom": 2.0, "energy_hartree": -1.0},
        {"target_distance_angstrom": 1.0, "energy_hartree": -0.5},
    ]


def test_parse_orca_scan_profile_rejects_invalid_rows(tmp_path):
    profile = tmp_path / "scan.relaxscanact.dat"
    profile.write_text("2.0 -1.0\ninvalid row\n")
    with pytest.raises(ValueError, match="row"):
        parse_orca_scan_profile(profile)


def test_recovered_scan_rejects_profile_off_requested_grid(tmp_path):
    molecule = Molecule(["H", "H"], np.asarray([[0., 0., 0.], [2., 0., 0.]]))
    request = OrcaBondScanRequest(0, 1, 2.0, 1.0, 2)
    scan_dir = tmp_path / "scan"
    scan_dir.mkdir()
    (scan_dir / "scan.out").write_text("****ORCA TERMINATED NORMALLY****")
    (scan_dir / "scan.allxyz").write_text(
        "2\nframe 1\nH 0 0 0\nH 2 0 0\n"
        "2\nframe 2\nH 0 0 0\nH 1 0 0\n"
    )
    (scan_dir / "scan.relaxscanact.dat").write_text("1.9 -1.0\n1.0 -0.5\n")
    result = load_orca_bond_scan_result(molecule, request, scan_dir)
    assert not result.success
    assert result.status == "scan_profile_grid_mismatch"


def test_scan_profile_writes_relative_energies_and_internal_maximum(tmp_path):
    from pyar.workflows.scan_bond import _write_scan_profile

    scan_dir = tmp_path / "orientation_000" / "scan"
    scan_dir.mkdir(parents=True)
    frames = [
        (["H", "H"], np.asarray([[0., 0., 0.], [2.0, 0., 0.]]), "first"),
        (["H", "H"], np.asarray([[0., 0., 0.], [1.5, 0., 0.]]), "maximum"),
        (["H", "H"], np.asarray([[0., 0., 0.], [1.0, 0., 0.]]), "last"),
    ]
    profile = [
        {"target_distance_angstrom": 2.0, "energy_hartree": -10.0},
        {"target_distance_angstrom": 1.5, "energy_hartree": -9.9},
        {"target_distance_angstrom": 1.0, "energy_hartree": -9.95},
    ]

    result = _write_scan_profile(scan_dir, frames, profile)
    assert result["maximum_scan_index"] == 2
    assert result["maximum_is_internal"]
    assert result["barrier_from_first_scan_point_kcal_mol"] == pytest.approx(
        0.1 * 627.509474
    )
    assert "not a confirmed transition state" in result["interpretation"]
    rows = (scan_dir / "scan_profile.csv").read_text().splitlines()
    assert len(rows) == 4
    metadata = json.loads(Path(result["ts_candidate_metadata"]).read_text())
    assert metadata["candidate_files"]["pre_maximum"].endswith("pre_maximum.xyz")
    assert metadata["candidate_files"]["post_maximum"].endswith("post_maximum.xyz")


def test_orca_keyword_maps_merged_uhf_to_unrestricted_dft():
    assert "UKS" in _orca_keyword({"method": "BP86", "basis": "def2-SVP"}, "uhf")


def test_orca_scan_rejects_incomplete_trajectory(tmp_path, monkeypatch):
    molecule = Molecule(["H", "H"], np.asarray([[0., 0., 0.], [2., 0., 0.]]))
    request = OrcaBondScanRequest(0, 1, 2.0, 1.0, 3)

    def fake_run(command, stdout_path=None, stderr_path=None):
        Path(stdout_path).write_text("****ORCA TERMINATED NORMALLY****\n")
        Path("scan.allxyz").write_text("2\nframe 1\nH 0 0 0\nH 2 0 0\n")
        return 0

    monkeypatch.setattr("pyar.backends.orca_scan.require_executable", lambda *args: "orca")
    monkeypatch.setattr("pyar.backends.orca_scan.run_command", fake_run)
    result = run_orca_bond_scan(
        molecule, request, tmp_path / "incomplete",
        {"method": "BP86", "basis": "def2-SVP", "nprocs": 1, "scf_cycles": 1000},
    )
    assert not result.success
    assert result.status == "trajectory_incomplete"


@pytest.mark.parametrize(
    "relaxed_identity,target_distance,identity_changed,target_contact,identity_status",
    [
        ({"inchi": "reactants", "smiles": "H.H"}, 0.7, False, True, "success"),
        ({"inchi": "product", "smiles": "H-H"}, 0.7, True, True, "success"),
        ({"inchi": "product", "smiles": "H-H"}, 3.0, True, False, "success"),
        (None, 0.7, None, True, "relaxed_identity_failed"),
    ],
)
def test_workflow_relaxes_scan_frame_and_writes_summary(
    tmp_path, monkeypatch, relaxed_identity, target_distance,
    identity_changed, target_contact, identity_status,
):
    a = tmp_path / "a.xyz"
    b = tmp_path / "b.xyz"
    write_xyz(["H"], [[0., 0., 0.]], a)
    write_xyz(["H"], [[1.8, 0., 0.]], b)
    scan_frame = np.asarray([[0., 0., 0.], [1.0, 0., 0.]])
    reactant_identity = {"inchi": "reactants", "smiles": "H.H"}

    def fake_scan(molecule, request, directory, qc_params):
        nonlocal scan_calls
        scan_calls += 1
        directory.mkdir(parents=True, exist_ok=True)
        input_path = directory / "scan.inp"
        output_path = directory / "scan.out"
        trajectory = directory / "scan_trajectory.xyz"
        final_path = directory / "final_scan.xyz"
        input_path.write_text("scan")
        output_path.write_text("normal")
        write_xyz(molecule.atoms_list, scan_frame, final_path)
        trajectory.write_text(
            "2\nframe 1\nH 0 0 0\nH 1.2 0 0\n"
            "2\nframe 2\nH 0 0 0\nH 1.0 0 0\n"
        )
        return OrcaBondScanResult(True, trajectory, final_path, scan_frame,
                                  input_path, output_path, "success",
                                  frames=[
                                      (["H", "H"], np.asarray([[0., 0., 0.], [1.2, 0., 0.]]), "frame 1"),
                                      (["H", "H"], scan_frame, "frame 2"),
                                  ],
                                  profile=[
                                      {"target_distance_angstrom": 1.2, "energy_hartree": -1.0},
                                      {"target_distance_angstrom": 1.0, "energy_hartree": -0.9},
                                  ])

    scan_calls = 0

    def fake_optimise(molecule, qc_params):
        molecule.coordinates = np.asarray([[0., 0., 0.], [target_distance, 0., 0.]])
        Path("job_relaxed").mkdir()
        write_xyz(molecule.atoms_list, molecule.coordinates,
                  Path("job_relaxed") / "result_relaxed.xyz")
        return True

    scan_workflow = importlib.import_module("pyar.workflows.scan_bond")
    monkeypatch.setattr(scan_workflow, "run_orca_bond_scan", fake_scan)
    monkeypatch.setattr(scan_workflow, "optimise", fake_optimise)
    monkeypatch.setattr(
        scan_workflow, "separated_reactant_identity", lambda *args: reactant_identity
    )
    if relaxed_identity is None:
        monkeypatch.setattr(
            scan_workflow, "molecule_identity_from_xyz",
            Mock(side_effect=ValueError("identity conversion failed")),
        )
    else:
        monkeypatch.setattr(
            scan_workflow, "molecule_identity_from_xyz", lambda *args: relaxed_identity
        )
    result = run_scan_bond(
        a, b, (0, 0), 1,
        {"software": "orca", "method": "BP86", "basis": "def2-SVP",
         "nprocs": 1, "scf_cycles": 1000},
        tmp_path / "scan_bond", scan_end=0.5, scan_points=2,
    )
    assert result["results"][0]["status"] == "scan_success_relax_success", result
    orientation_result = result["results"][0]
    assert orientation_result["target_bond_present_after_relaxation"] is target_contact
    assert orientation_result["identity_status"] == identity_status
    assert orientation_result["reactant_identity"] == reactant_identity
    if identity_status == "success":
        assert orientation_result["relaxed_identity"] == relaxed_identity
        assert orientation_result["product_identity_changed"] is identity_changed
    else:
        assert orientation_result["product_identity_changed"] is None
        assert "identity conversion failed" in orientation_result["identity_error"]
    assert (tmp_path / "scan_bond" / "summary.json").exists()
    assert (tmp_path / "scan_bond" / "orientation_000" / "result_relaxed.xyz").exists()
    assert (tmp_path / "scan_bond" / "orientation_000" / "result.json").exists()
    assert (tmp_path / "scan_bond" / "orientation_000" / "scan" / "scan_profile.csv").exists()
    assert (tmp_path / "scan_bond" / "orientation_000" / "ts_candidates" / "highest_scan_energy.xyz").exists()
    rerun = run_scan_bond(
        a, b, (0, 0), 1,
        {"software": "orca", "method": "BP86", "basis": "def2-SVP",
         "nprocs": 1, "scf_cycles": 1000},
        tmp_path / "scan_bond", scan_end=0.5, scan_points=2,
    )
    assert scan_calls == 1
    assert rerun["results"][0]["status"] == "scan_success_relax_success"


def test_restart_reuses_completed_scan_after_relaxation_failure(tmp_path, monkeypatch):
    a = tmp_path / "a.xyz"
    b = tmp_path / "b.xyz"
    write_xyz(["H"], [[0., 0., 0.]], a)
    write_xyz(["H"], [[1.8, 0., 0.]], b)
    scan_frame = np.asarray([[0., 0., 0.], [1.0, 0., 0.]])
    scan_calls = 0
    relax_calls = 0

    def fake_scan(molecule, request, directory, qc_params):
        nonlocal scan_calls
        scan_calls += 1
        directory.mkdir(parents=True, exist_ok=True)
        input_path = directory / "scan.inp"
        output_path = directory / "scan.out"
        trajectory = directory / "scan_trajectory.xyz"
        final_path = directory / "final_scan.xyz"
        input_path.write_text("scan")
        output_path.write_text("normal")
        write_xyz(molecule.atoms_list, scan_frame, final_path)
        trajectory.write_text("2\nframe\nH 0 0 0\nH 1 0 0\n")
        frames = [
            (["H", "H"], np.asarray([[0., 0., 0.], [1.5, 0., 0.]]), "frame 1"),
            (["H", "H"], scan_frame, "frame 2"),
        ]
        profile = [
            {"target_distance_angstrom": 1.5, "energy_hartree": -1.1},
            {"target_distance_angstrom": 1.0, "energy_hartree": -1.0},
        ]
        return OrcaBondScanResult(True, trajectory, final_path, scan_frame,
                                  input_path, output_path, "success",
                                  frames=frames, profile=profile)

    def fake_optimise(molecule, qc_params):
        nonlocal relax_calls
        relax_calls += 1
        if relax_calls == 1:
            return False
        Path("job_relaxed").mkdir()
        write_xyz(molecule.atoms_list, molecule.coordinates,
                  Path("job_relaxed") / "result_relaxed.xyz")
        return True

    scan_workflow = importlib.import_module("pyar.workflows.scan_bond")
    monkeypatch.setattr(scan_workflow, "run_orca_bond_scan", fake_scan)
    monkeypatch.setattr(scan_workflow, "optimise", fake_optimise)
    monkeypatch.setattr(scan_workflow, "separated_reactant_identity",
                        lambda *args: {"inchi": "reactants", "smiles": "H.H"})
    monkeypatch.setattr(scan_workflow, "molecule_identity_from_xyz",
                        lambda *args: {"inchi": "product", "smiles": "H-H"})
    params = {"software": "orca", "method": "BP86", "basis": "def2-SVP",
              "nprocs": 1, "scf_cycles": 1000}
    output = tmp_path / "scan_bond"
    first = run_scan_bond(a, b, (0, 0), 1, params, output, scan_points=2)
    assert first["results"][0]["status"] == "scan_success_relax_failed"
    completed_scan = OrcaBondScanResult(
        True, output / "orientation_000" / "scan" / "scan_trajectory.xyz",
        output / "orientation_000" / "scan" / "final_scan.xyz", scan_frame,
        output / "orientation_000" / "scan" / "scan.inp",
        output / "orientation_000" / "scan" / "scan.out", "success",
        frames=[
            (["H", "H"], np.asarray([[0., 0., 0.], [1.5, 0., 0.]]), "frame 1"),
            (["H", "H"], scan_frame, "frame 2"),
        ],
        profile=[
            {"target_distance_angstrom": 1.5, "energy_hartree": -1.1},
            {"target_distance_angstrom": 1.0, "energy_hartree": -1.0},
        ],
    )
    monkeypatch.setattr(scan_workflow, "load_orca_bond_scan_result",
                        lambda *args: completed_scan)
    second = run_scan_bond(a, b, (0, 0), 1, params, output, scan_points=2)
    assert second["results"][0]["status"] == "scan_success_relax_success"
    assert scan_calls == 1
    assert relax_calls == 2


def test_restart_rejects_changed_request(tmp_path, monkeypatch):
    a = tmp_path / "a.xyz"
    b = tmp_path / "b.xyz"
    write_xyz(["H"], [[0., 0., 0.]], a)
    write_xyz(["H"], [[1.8, 0., 0.]], b)
    workflow = importlib.import_module("pyar.workflows.scan_bond")
    monkeypatch.setattr(workflow, "separated_reactant_identity",
                        lambda *args: {"inchi": "reactants", "smiles": "H.H"})
    def fail_without_running_orca(molecule, request, directory, qc_params):
        directory.mkdir(parents=True, exist_ok=True)
        return OrcaBondScanResult(False, None, None, None,
                                  directory / "scan.inp", directory / "scan.out", "mock_failure")
    monkeypatch.setattr(workflow, "run_orca_bond_scan", fail_without_running_orca)
    params = {"software": "orca", "method": "BP86", "basis": "def2-SVP",
              "nprocs": 1, "scf_cycles": 1000}
    output = tmp_path / "scan_bond"
    run_scan_bond(a, b, (0, 0), 1, params, output, scan_end=0.5, scan_points=2)
    with pytest.raises(FileExistsError, match="different request"):
        run_scan_bond(a, b, (0, 0), 1, params, output, scan_end=0.6, scan_points=2)


def test_workflow_default_scan_end_is_0_8_times_covalent_radii(tmp_path, monkeypatch):
    a = tmp_path / "a.xyz"
    b = tmp_path / "b.xyz"
    write_xyz(["H"], [[0., 0., 0.]], a)
    write_xyz(["H"], [[1.8, 0., 0.]], b)
    observed = {}

    def fake_scan(molecule, request, directory, qc_params):
        observed["end_distance"] = request.end_distance_angstrom
        observed["n_points"] = request.n_points
        directory.mkdir(parents=True, exist_ok=True)
        input_path = directory / "scan.inp"
        output_path = directory / "scan.out"
        trajectory = directory / "scan_trajectory.xyz"
        final_path = directory / "final_scan.xyz"
        input_path.write_text("scan")
        output_path.write_text("normal")
        coordinates = np.asarray([[0., 0., 0.], [request.end_distance_angstrom, 0., 0.]])
        write_xyz(molecule.atoms_list, coordinates, final_path)
        fractions = np.linspace(0.0, 1.0, request.n_points)
        frames = []
        trajectory_text = []
        profile = []
        for frame_index, fraction in enumerate(fractions):
            distance = 1.5 * sum(Molecule.from_xyz(path).covalent_radius[0] for path in (a, b))
            distance = distance + fraction * (request.end_distance_angstrom - distance)
            coordinates_i = np.asarray([[0., 0., 0.], [distance, 0., 0.]])
            frames.append((["H", "H"], coordinates_i, f"frame {frame_index + 1}"))
            trajectory_text.append(f"2\nframe {frame_index + 1}\nH 0 0 0\nH {distance} 0 0\n")
            profile.append({"target_distance_angstrom": distance,
                            "energy_hartree": -1.0 + 0.001 * frame_index})
        coordinates = frames[-1][1]
        trajectory.write_text("".join(trajectory_text))
        return OrcaBondScanResult(True, trajectory, final_path, coordinates,
                                  input_path, output_path, "success", frames=frames,
                                  profile=profile)

    def fake_optimise(molecule, qc_params):
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
        tmp_path / "scan_bond",
    )
    radii_sum = sum(Molecule.from_xyz(path).covalent_radius[0] for path in (a, b))
    assert observed["end_distance"] == pytest.approx(0.8 * radii_sum)
    assert result["results"][0]["target_distance_scan_end_angstrom"] == pytest.approx(
        0.8 * radii_sum
    )
    request = json.loads((tmp_path / "scan_bond" / "request.json").read_text())
    assert request["default_scan_end_factor"] == pytest.approx(0.8)
    assert request["default_scan_end_angstrom"] == pytest.approx(0.8 * radii_sum)
