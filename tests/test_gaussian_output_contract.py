"""Synthetic Gaussian log boundaries; these tests do not run Gaussian."""

from pathlib import Path
import subprocess

import numpy as np
import pytest

from pyar.backend_errors import BackendExecutionError
from pyar.backends import gaussian
from pyar.core.molecule import Molecule


def log(energy="-1.25", second_z="0.74"):
    return (
        " SCF Done: E(RHF) = -1.00 A.U.\n"
        f" SCF Done: E(RHF) = {energy} A.U.\n"
        " Optimization completed\n Stationary point found.\n"
        " Standard orientation:\n ----\n Center Atomic Atomic Coordinates\n"
        " Number Number Type X Y Z\n ----\n"
        " 1 1 0 0.0 0.0 0.0\n"
        f" 2 1 0 0.0 0.0 {second_z}\n"
        " Normal termination of Gaussian\n"
    )


@pytest.fixture
def subject(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(gaussian, "require_executable", lambda *args: "synthetic-g16")
    molecule = Molecule(
        ["H", "H"],
        [[0.0, 0.0, 0.0], [0.0, 0.0, 0.74]],
        name="probe",
        charge=0,
        multiplicity=1,
    )
    return gaussian.Gaussian(
        molecule, dict(method="HF", basis="sto-3g", nprocs=1, scf_cycles=10)
    )


def process_output(monkeypatch, subject, payload, status=0):
    calls = []

    def execute(command, **kwargs):
        calls.append((command, kwargs))
        Path(kwargs["stdout_path"]).write_text(payload)
        return status

    monkeypatch.setattr(gaussian, "run_command", execute)
    return calls


def test_converged_native_log_preserves_energy_units_coordinates_and_command(
    subject, monkeypatch
):
    calls = process_output(monkeypatch, subject, log())
    assert subject.optimize() is True
    assert (
        subject.energy == -1.25
    )  # Synthetic parser witness, not a chemical reference.
    np.testing.assert_array_equal(
        subject.optimized_coordinates, [[0.0, 0.0, 0.0], [0.0, 0.0, 0.74]]
    )
    assert calls == [
        (
            ["synthetic-g16"],
            dict(
                stdout_path=subject.out_file,
                stderr_path=subject.out_file,
                stdin_path=subject.inp_file,
            ),
        )
    ]
    rows = Path(subject.result_xyz_file).read_text().splitlines()
    assert rows[:2] == ["  2", "probe:-1.25"]
    assert [row.split()[0] for row in rows[2:]] == ["H", "H"]
    np.testing.assert_array_equal(
        [[float(v) for v in row.split()[1:]] for row in rows[2:]],
        subject.optimized_coordinates,
    )
    assert "0 1\n" in Path(subject.inp_file).read_text()


@pytest.mark.parametrize(
    "kind",
    [
        "empty",
        "no-scf",
        "no-optimization",
        "no-normal",
        "trailing-text",
        "bad-energy",
        "truncated-coordinates",
    ],
)
def test_failed_or_malformed_log_cannot_replace_previous_result(
    subject, monkeypatch, kind
):
    payload = log()
    if kind == "empty":
        payload = ""
    elif kind == "no-scf":
        payload = (
            "\n".join(row for row in payload.splitlines() if "SCF Done" not in row)
            + "\n"
        )
    elif kind == "no-optimization":
        payload = payload.replace(
            "Optimization completed", "Optimization not completed"
        )
    elif kind == "no-normal":
        payload = payload.replace("Normal termination", "Error termination")
    elif kind == "trailing-text":
        payload += "unexpected trailing output\n"
    elif kind == "bad-energy":
        payload = log(energy="not-a-number")
    else:
        payload = payload.replace("2 1 0 0.0 0.0 0.74", "2 1 0 0.0")
    Path(subject.result_xyz_file).write_bytes(b"previous result witness\n")
    process_output(monkeypatch, subject, payload)
    with pytest.raises(BackendExecutionError):
        subject.optimize()
    assert Path(subject.result_xyz_file).read_bytes() == b"previous result witness\n"


@pytest.mark.parametrize(
    "payload",
    [
        log(energy="nan"),
        log(energy="inf"),
        log(second_z="nan"),
        log().replace("Standard orientation", "Unavailable orientation"),
        log().replace(" 1 1 0 0.0 0.0 0.0\n", "").replace(" 2 1 0 0.0 0.0 0.74\n", ""),
    ],
    ids=[
        "nan-energy",
        "infinite-energy",
        "nan-coordinate",
        "missing-orientation",
        "missing-coordinate-rows",
    ],
)
def test_invalid_result_is_rejected_before_publication_or_state_assignment(
    subject, monkeypatch, payload
):
    Path(subject.result_xyz_file).write_bytes(b"previous result witness\n")
    process_output(monkeypatch, subject, payload)
    with pytest.raises(BackendExecutionError):
        subject.optimize()
    assert Path(subject.result_xyz_file).read_bytes() == b"previous result witness\n"
    assert subject.energy == 0.0
    assert subject.optimized_coordinates == []


def test_nonzero_process_exit_cannot_adopt_valid_looking_stale_output(
    subject, monkeypatch
):
    Path(subject.result_xyz_file).write_bytes(b"previous result witness\n")
    process_output(monkeypatch, subject, log(), status=2)
    with pytest.raises(BackendExecutionError):
        subject.optimize()
    assert Path(subject.result_xyz_file).read_bytes() == b"previous result witness\n"


@pytest.mark.parametrize(
    "error",
    [
        FileNotFoundError("g16 missing"),
        PermissionError("g16 not executable"),
        subprocess.TimeoutExpired(["g16"], 1),
    ],
)
def test_process_errors_remain_typed_and_keep_the_original_cause(
    subject, monkeypatch, error
):
    def execute(*args, **kwargs):
        raise error

    monkeypatch.setattr(gaussian, "run_command", execute)
    with pytest.raises(BackendExecutionError) as caught:
        subject.optimize()
    assert caught.value.__cause__ is error
    assert not Path(subject.result_xyz_file).exists()
