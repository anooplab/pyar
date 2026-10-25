"""Fixed-geometry contracts, failure boundaries and scientific qualification."""

import json
from pathlib import Path

import numpy as np
import pytest

from pyar import modern_cli
from pyar.scripts.singlepoint import build_parser
from pyar.singlepoint import request, native, service, vibrations, worker
from pyar.singlepoint.io import validate_result


@pytest.fixture
def geometries(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    paths = []
    for name, length in [("a", 0.74), ("b", 0.85)]:
        path = tmp_path / f"{name}.xyz"
        path.write_text(f"2\n{name}\nH 0 0 0\nH {length} 0 0\n")
        paths.append(str(path))
    return paths


@pytest.fixture
def harmonic(monkeypatch):
    monkeypatch.setattr(request, "preflight", lambda _: [])
    monkeypatch.setattr("pyar.scripts.singlepoint.preflight", lambda _: [])

    def evaluate(settings, molecule, coordinates, directory, gradient):
        return float(np.sum(coordinates**2) / 2), coordinates if gradient else None, {}

    monkeypatch.setitem(native.EVALUATORS, "xtb", evaluate)
    return evaluate


@pytest.mark.parametrize(
    "command", ["energy", "engrad", "engra", "gradient", "frequency"]
)
def test_batch_layout_and_json(geometries, harmonic, command, capsys):
    args = [command, *geometries, "--json"]
    if command == "frequency":
        with pytest.raises(SystemExit) as error:
            modern_cli.main(args)
        assert error.value.code == 1
    else:
        modern_cli.main(args)
    payload = json.loads(capsys.readouterr().out)
    result = payload["result"]
    assert result["input_count"] == 2
    assert result["completed_count"] == 2
    root = Path(result["run_directory"])
    assert (root / "summary.csv").is_file()
    assert (root / "pyar-run.toml").is_file()
    for name in ("a", "b"):
        job = root / f"job_{name}"
        assert (job / "input.xyz").is_file()
        data = json.loads((job / "result.json").read_text())
        assert data["calculation_status"] == "completed"
        assert (job / "backend/reference/evaluation.json").is_file()
        if command == "frequency":
            assert (job / "hessian.txt").is_file()
            assert data["qualification"] == "nonstationary"
            assert data["rigid_mode_count"] == 5


@pytest.mark.parametrize("command", ["energy", "engrad", "frequency"])
def test_check_and_dry_run_no_writes(geometries, harmonic, command, capsys):
    before = set(Path.cwd().rglob("*"))
    modern_cli.main([command, *geometries, "--dry-run", "--json"])
    report = json.loads(capsys.readouterr().out)
    assert report["status"] == "ready"
    assert report["plan"]["backend"] == "xtb"
    assert set(Path.cwd().rglob("*")) == before


@pytest.mark.parametrize(
    "options, message",
    [
        (["--backend", "bogus"], "Unsupported"),
        (["--method", "PBE"], "does not support"),
        (["--basis", "STO-3G"], "does not support"),
        (["--xtb-model", "bogus"], "xtb-model"),
        (["--nprocs", "0"], "positive"),
        (["--scf-cycles", "-1"], "positive"),
        (["--accuracy", "nan"], "positive"),
        (["--memory-mb", "20"], "only for ORCA"),
        (["--solvent-model", "gbsa"], "requires --solvent"),
        (["--solvent", "water", "--solvent-model", "smd"], "Choose a solvent"),
        (["--solvent", "water\n--opt"], "single backend keyword"),
        (["--dispersion", "D4"], "only for ORCA"),
        (["--backend", "orca", "--method", "PBE Opt"], "single backend keyword"),
        (["--backend", "orca", "--method", "Opt"], "workflow or execution"),
        (["--backend", "orca", "--basis", "NumFreq"], "workflow or execution"),
        (["--backend", "orca", "--method", "gxtb"], "external g-xTB"),
        (
            ["--backend", "orca", "--method", "r2scan-3c", "--basis", "STO-3G"],
            "own basis",
        ),
        (["--backend", "gaussian", "--dispersion", "D4"], "not D4"),
        (["--backend", "gaussian", "--accuracy", "1"], "only for xTB"),
        (["--backend", "ani", "--model", "bad"], "Choose --model"),
        (["--backend", "ani", "--scf-cycles", "10"], "ML potentials"),
        (["--backend", "ani", "--solvent", "water"], "does not support solvent"),
        (["--backend", "turbomole", "--nprocs", "2"], "nprocs 1"),
    ],
)
def test_invalid_request_before_execution(geometries, options, message):
    parser = build_parser("energy")
    with pytest.raises(ValueError, match=message):
        request.resolve(parser.parse_args([*geometries, *options]), "energy")
    assert not Path("energy").exists()


@pytest.mark.parametrize(
    "option", ["step", "max-gradient", "rms-gradient", "imaginary-threshold"]
)
def test_invalid_frequency_controls(geometries, option):
    args = build_parser("frequency").parse_args([*geometries, "--" + option, "0"])
    with pytest.raises(ValueError, match="positive"):
        request.resolve(args, "frequency")


def test_defaults_and_composite_basis(geometries):
    parser = build_parser()
    settings, molecules, root = request.resolve(parser.parse_args(geometries), "energy")
    assert settings["xtb_model"] == "gfn2"
    assert settings["model"] is None
    assert [m.multiplicity for m in molecules] == [1, 1]
    args = parser.parse_args(
        [*geometries, "--backend", "orca", "--method", "r2scan-3c"]
    )
    settings, _, _ = request.resolve(args, "energy")
    assert settings["basis"] is None
    args = parser.parse_args([*geometries, "--backend", "orca", "--solvent", "water"])
    assert request.resolve(args, "energy")[0]["solvent_model"] == "CPCM"


def test_duplicate_and_existing_output(geometries):
    parser = build_parser()
    with pytest.raises(ValueError, match="Duplicate job"):
        request.resolve(parser.parse_args([geometries[0]] * 2), "energy")
    Path("energy").mkdir()
    with pytest.raises(ValueError, match="already exists"):
        request.resolve(parser.parse_args(geometries), "energy")


def test_preflight_validates_all_inputs(geometries, harmonic):
    Path(geometries[1]).write_text("invalid\n")
    with pytest.raises(SystemExit):
        modern_cli.main(["energy", *geometries])
    assert not Path("energy").exists()


def test_partial_failure_continues_and_records(
    geometries, harmonic, monkeypatch, capsys
):
    def evaluate(settings, molecule, *args):
        if molecule.name == "a":
            raise RuntimeError("SCF failed")
        return harmonic(settings, molecule, *args)

    monkeypatch.setitem(native.EVALUATORS, "xtb", evaluate)
    with pytest.raises(SystemExit) as error:
        modern_cli.main(["engrad", *geometries, "--json"])
    assert error.value.code == 1
    report = json.loads(capsys.readouterr().out)["result"]
    assert report["status"] == "partial_failed"
    assert report["failed_count"] == 1
    assert report["results"][1]["calculation_status"] == "completed"
    assert Path("engrad/job_a/backend/reference").is_dir()


def test_config_roundtrip(geometries, harmonic):
    modern_cli.main(["engrad", *geometries, "--check", "--write-config", "engrad.toml"])
    modern_cli.main(["run", "engrad.toml"])
    assert Path("engrad/job_a/gradient.txt").is_file()


@pytest.mark.parametrize(
    "energy, gradient, required",
    [
        (float("nan"), None, False),
        (float("inf"), None, False),
        (1.0, None, True),
        (1.0, np.zeros((1, 3)), True),
        (1.0, np.full((2, 3), np.nan), True),
    ],
)
def test_result_validation(energy, gradient, required):
    with pytest.raises(ValueError):
        validate_result(energy, gradient, 2, required)


@pytest.mark.parametrize(
    "stationary, imaginary, expected",
    [
        (False, 0, "nonstationary"),
        (True, 0, "minimum_confirmed"),
        (True, 1, "first_order_saddle_candidate"),
        (True, 2, "higher_order_saddle_candidate"),
    ],
)
def test_qualification(stationary, imaginary, expected):
    assert vibrations.qualification(stationary, imaginary) == expected


@pytest.mark.scientific
@pytest.mark.tier1
def test_numerical_hessian_and_frequency_units():
    # Diatomic stretch E = 1/2 k (x2-x1-r0)^2. Analytic Hessian has one
    # stretch eigenvalue 2*k/m after removing 5 rigid motions.
    k = 0.4
    coordinates = np.array([[0.0, 0.0, 0.0], [1.4, 0.0, 0.0]])

    def evaluate(c, label):
        dx = c[1, 0] - c[0, 0] - 1.4
        return 0.5 * k * dx**2, np.array([[-k * dx, 0, 0], [k * dx, 0, 0]])

    hessian, asymmetry = vibrations.numerical_hessian(evaluate, coordinates, 0.005)
    assert hessian[0, 3] == pytest.approx(-k)
    assert asymmetry < 1e-12
    analysis = vibrations.analyze_hessian(hessian, coordinates, [1.0, 1.0], 20.0)
    assert analysis["frequencies_cm-1"] == pytest.approx(
        [np.sqrt(2 * k) * 5140.48714], rel=1e-7
    )
    assert analysis["rigid_mode_count"] == 5
    assert analysis["imaginary_frequency_count"] == 0
    unstable = vibrations.analyze_hessian(-hessian, coordinates, [1.0, 1.0], 20.0)
    assert unstable["imaginary_frequency_count"] == 1


@pytest.mark.parametrize(
    "coordinates, expected",
    [
        ([[0.0, 0.0, 0.0]], 3),
        ([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]], 5),
        ([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], 6),
    ],
)
def test_rigid_modes(coordinates, expected):
    basis, rank = vibrations.vibrational_basis(
        np.array(coordinates), np.ones(len(coordinates))
    )
    assert rank == expected
    assert basis.shape == (3 * len(coordinates), 3 * len(coordinates) - rank)


def test_invalid_mass():
    with pytest.raises(ValueError, match="masses"):
        vibrations.analyze_hessian(np.eye(3), np.zeros((1, 3)), [0.0], 20.0)


def test_missing_dependencies(monkeypatch):
    monkeypatch.setattr(request.shutil, "which", lambda _: None)
    with pytest.raises(ValueError, match="Missing executable 'orca'"):
        request.preflight({"software": "orca"})
    monkeypatch.setattr(request, "module_available", lambda _: False)
    with pytest.raises(ValueError, match="Missing Python package"):
        request.preflight({"software": "ani"})


def test_inputs_honor_requested_theory(geometries):
    args = build_parser().parse_args(
        [
            *geometries,
            "--backend",
            "orca",
            "--method",
            "PBE0",
            "--basis",
            "def2-TZVP",
            "--solvent",
            "water",
        ]
    )
    settings, molecules, _ = request.resolve(args, "energy")
    text = native.orca_input(settings, molecules[0], np.zeros((2, 3)), True)
    assert "! PBE0 def2-TZVP TightSCF NoAutoStart EnGrad CPCM(water)" in text
    assert "D3" not in text and "RI " not in text and "Opt" not in text
    args = build_parser().parse_args(
        [
            *geometries,
            "--backend",
            "gaussian",
            "--solvent",
            "water",
            "--dispersion",
            "D3BJ",
        ]
    )
    settings, molecules, _ = request.resolve(args, "energy")
    text = native.gaussian_input(settings, molecules[0], np.zeros((2, 3)), True)
    assert "PBE1PBE/def2SVP" in text
    assert "Force" in text and "SCRF=(PCM,Solvent=water)" in text
    assert "EmpiricalDispersion=GD3BJ" in text


@pytest.mark.parametrize(
    "text", ["no completion", "ORCA TERMINATED NORMALLY\nSCF NOT CONVERGED"]
)
def test_completion_requires_more_than_exit_zero(text):
    with pytest.raises(RuntimeError):
        native.require_completion(text, "ORCA TERMINATED NORMALLY")


def test_subprocess_evidence_on_failure(tmp_path):
    import sys

    with pytest.raises(RuntimeError, match="status 3"):
        native.run_program(
            [sys.executable, "-c", 'import sys; print("bad"); sys.exit(3)'],
            tmp_path,
            {"nprocs": 1},
        )
    assert "bad" in (tmp_path / "backend.out").read_text()
    assert json.loads((tmp_path / "backend_exit.json").read_text())["returncode"] == 3


def test_turbomole_parameter_translation(monkeypatch):
    import ase.calculators.turbomole

    captured = {}
    monkeypatch.setattr(
        ase.calculators.turbomole, "Turbomole", lambda **kw: captured.update(kw)
    )
    worker.turbomole_calculator(
        {
            "charge": -1,
            "multiplicity": 2,
            "settings": {"method": "PBE0", "basis": "def2-TZVP", "scf_cycles": 123},
        }
    )
    assert captured["density functional"] == "pbe0"
    assert captured["total charge"] == -1
    assert captured["multiplicity"] == 2
    assert captured["scf iterations"] == 123


def test_small_negative_modes_are_retained():
    coordinates = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    hessian = np.zeros((6, 6))
    hessian[0, 0] = hessian[3, 3] = -1e-7
    hessian[0, 3] = hessian[3, 0] = 1e-7
    result = vibrations.analyze_hessian(hessian, coordinates, [1.0, 1.0], 20.0)
    assert result["imaginary_frequency_count"] == 0
    assert len(result["small_negative_frequencies_cm-1"]) == 1


def test_frequency_at_harmonic_minimum(tmp_path, monkeypatch):
    molecule = __import__("pyar.core.molecule", fromlist=["Molecule"]).Molecule(
        ["H", "H"], [[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]], name="h2"
    )
    from pyar.singlepoint.io import BOHR_ANGSTROM

    coordinates = molecule.coordinates / BOHR_ANGSTROM
    distance = coordinates[1, 0]

    def evaluate(c, label):
        displacement = c[1, 0] - c[0, 0] - distance
        return 0.5 * displacement**2, np.array(
            [[-displacement, 0.0, 0.0], [displacement, 0.0, 0.0]]
        )

    result = service.frequency_result(
        evaluate,
        molecule,
        coordinates,
        np.zeros((2, 3)),
        tmp_path,
        {
            "max_gradient": 4.5e-4,
            "rms_gradient": 3e-4,
            "step": 0.005,
            "imaginary_threshold": 20.0,
        },
    )
    assert result["minimum_confirmed"]
    assert result["warning"] is None


def test_summaries_distinguish_failure_and_scientific_gate(tmp_path):
    assert service.batch_status(2, 0, 2) == "failed"
    assert service.batch_status(1, 0, 2) == "partial_failed"
    assert service.batch_status(0, 1, 2) == "scientific_gate_failed"
    assert service.batch_status(0, 0, 2) == "completed"


@pytest.mark.parametrize("name", ["a b", "..", "a/b"])
def test_invalid_names(name):
    with pytest.raises(ValueError, match="basename"):
        request.validate_name(name)


def test_profile_roundtrip(geometries, harmonic, monkeypatch):
    # A profile saves resolved theory, not output data or numerical results.
    from pyar.run_config import PROFILE_KIND, write_toml_document

    monkeypatch.setenv("XDG_CONFIG_HOME", str(Path.cwd() / "config"))
    profile = Path("config/pyar/profiles/sp.toml")
    write_toml_document(
        profile,
        kind=PROFILE_KIND,
        workflow="engrad",
        arguments=[],
        effective={"backend": "xtb", "xtb_model": "gfn2"},
    )
    modern_cli.main(["engrad", *geometries, "--profile", "sp", "--check"])
    assert not Path("engrad").exists()


def test_nonfinite_hessian_rejected():
    def evaluate(c, label):
        return 1.0, np.full((1, 3), np.nan)

    with pytest.raises(ValueError, match="Non-finite"):
        vibrations.numerical_hessian(evaluate, np.zeros((1, 3)), 0.005)


@pytest.mark.scientific
@pytest.mark.tier1
def test_nonlinear_frequency_rotation_and_translation_invariance():
    coordinates = np.array([[0.0, 0.0, 0.0], [0.957, 0.0, 0.0], [-0.24, 0.927, 0.0]])
    masses = np.array([16.0, 1.0, 1.0])
    basis, rank = vibrations.vibrational_basis(coordinates, masses)
    assert rank == 6
    weights = np.repeat(np.sqrt(masses), 3)
    hessian = (basis @ np.diag([0.1, 0.2, 0.3]) @ basis.T) * np.outer(weights, weights)
    original = vibrations.analyze_hessian(hessian, coordinates, masses, 20.0)
    rotation = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    transform = np.kron(np.eye(3), rotation)
    moved = vibrations.analyze_hessian(
        transform @ hessian @ transform.T,
        coordinates @ rotation.T + [4.0, -3.0, 2.0],
        masses,
        20.0,
    )
    assert original["frequencies_cm-1"] == pytest.approx(moved["frequencies_cm-1"])
