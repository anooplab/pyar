"""Synthetic native-output boundaries, not a Turbomole calculation campaign."""

from pathlib import Path
import subprocess
from unittest.mock import Mock

import numpy as np
import pytest

from pyar.backend_errors import (
    BackendConvergenceError,
    BackendExecutionError,
    BackendGeometryError,
    BackendGradientError,
    BackendOutputError,
    BackendSCFError,
)
from pyar.backends import turbomole as tm
from pyar.constants import BOHR_ANGSTROM_LEGACY


@pytest.fixture
def subject(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    runner = tm.Turbomole.__new__(tm.Turbomole)
    runner.atoms_list = ["H", "He"]
    runner.start_coords = np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
    runner.number_of_atoms = 2
    runner.atoms_in_fragments = [[0], [1]]
    runner.job_name = "synthetic"
    runner.result_xyz_file = "result_synthetic.xyz"
    runner.basis, runner.method = "def2-SVP", "bp86"
    runner.charge, runner.multiplicity, runner.gamma = 1, 1, None
    runner.energy = runner.optimized_coordinates = None
    tm.make_coord(runner.atoms_list, runner.start_coords)
    Path("energy").write_text("$energy\n1 -1.25 0 0\n$end\n")
    return runner


def jobex_output(monkeypatch, *, flag=None, message="", failure=None):
    def execute(command, **kwargs):
        assert command == ["jobex", "-ri", "-c", "7"]
        kwargs["stdout"].write("synthetic jobex output\n")
        if flag:
            Path(flag).write_text(message)
        if failure:
            raise failure
        return 0

    monkeypatch.setattr(tm.subp, "check_call", execute)


@pytest.mark.parametrize(
    "message,error",
    [
        ("ERROR: Module", BackendExecutionError),
        ("ERROR in statpt step,", BackendGeometryError),
        ("ERROR in relax step,", BackendGeometryError),
        ("ERROR: your energy calculation did not converge !!,", BackendSCFError),
        ("ERROR in dscf step", BackendSCFError),
        ("ERROR in gradient step", BackendGradientError),
        ("unrecognized native failure", BackendExecutionError),
    ],
)
def test_jobex_classifies_native_failures_without_publishing(
    subject, monkeypatch, message, error
):
    jobex_output(monkeypatch, flag="GEO_OPT_FAILED", message=message)
    with pytest.raises(error):
        subject.run_turbomole_jobex(max_cycles=7)
    assert subject.energy is None
    assert not Path(subject.result_xyz_file).exists()
    assert Path("GEO_OPT_FAILED").read_text() == message


def test_convergence_requires_native_evidence_not_just_exit_zero(subject, monkeypatch):
    jobex_output(monkeypatch)
    with pytest.raises(BackendOutputError, match="no typed result"):
        subject.run_turbomole_jobex(max_cycles=7)
    assert not Path(subject.result_xyz_file).exists()


def test_converged_result_uses_hartree_energy_and_bohr_coordinates(
    subject, monkeypatch
):
    jobex_output(monkeypatch, flag="GEO_OPT_CONVERGED")
    assert subject.run_turbomole_jobex(max_cycles=7) is True
    assert (
        subject.energy == -1.25
    )  # Synthetic parser witness, not a chemical reference.
    np.testing.assert_array_equal(
        subject.optimized_coordinates, subject.start_coords * BOHR_ANGSTROM_LEGACY
    )
    assert Path(subject.result_xyz_file).is_file()


def test_cycle_limit_exposes_only_an_incomplete_frame(subject, monkeypatch):
    jobex_output(
        monkeypatch, flag="GEO_OPT_FAILED", message="OPTIMIZATION DID NOT CONVERGE"
    )
    with pytest.raises(BackendConvergenceError) as caught:
        subject.run_turbomole_jobex(max_cycles=7)
    assert caught.value.usable_geometry is True
    assert subject.energy == -1.25
    assert not Path(subject.result_xyz_file).exists()


@pytest.mark.parametrize(
    "failure",
    [FileNotFoundError("jobex missing"), subprocess.CalledProcessError(2, ["jobex"])],
)
def test_process_failure_cannot_reuse_a_stale_convergence_marker(
    subject, monkeypatch, failure
):
    Path("GEO_OPT_CONVERGED").write_text("old convergence witness")
    jobex_output(monkeypatch, failure=failure)
    with pytest.raises(BackendExecutionError):
        subject.run_turbomole_jobex(max_cycles=7)
    assert not Path(subject.result_xyz_file).exists()
    assert Path("GEO_OPT_CONVERGED").read_text() == "old convergence witness"


@pytest.mark.parametrize(
    "filename,payload",
    [
        ("energy", "truncated"),
        ("coord", "bad header\n"),
        ("coord", "$coord\ninvalid row\n"),
    ],
)
def test_malformed_native_output_raises_inside_typed_boundary(
    subject, monkeypatch, filename, payload
):
    Path(filename).write_text(payload)
    jobex_output(monkeypatch, flag="GEO_OPT_CONVERGED")
    with pytest.raises(BackendExecutionError):
        subject.run_turbomole_jobex(max_cycles=7)
    assert not Path(subject.result_xyz_file).exists()


@pytest.mark.parametrize("statuses", [(False, True), (False, False)])
def test_define_cartesian_fallback_preserves_electronic_state(
    subject, monkeypatch, statuses
):
    prepare = Mock(side_effect=statuses)
    monkeypatch.setattr(tm, "prepare_control", prepare)
    subject.run_turbomole_jobex = Mock(return_value=True)
    if statuses[-1]:
        assert subject.optimize() is True
        subject.run_turbomole_jobex.assert_called_once_with(max_cycles=200)
    else:
        with pytest.raises(BackendGeometryError):
            subject.optimize()
        subject.run_turbomole_jobex.assert_not_called()
    assert prepare.call_count == 2
    first, second = [call.kwargs for call in prepare.call_args_list]
    assert second == dict(first, coordinates="cartesian")
    assert (first["charge"], first["multiplicity"], first["func"], first["basis"]) == (
        1,
        1,
        "bp",
        "def2-SVP",
    )


@pytest.mark.parametrize(
    "step,error",
    [
        ("energy", BackendSCFError),
        ("gradient", BackendGradientError),
        ("update", BackendGeometryError),
    ],
)
def test_biased_loop_failure_is_not_success(subject, monkeypatch, step, error):
    subject.gamma = 100.0
    Path("control").write_text("$end\n")
    Path("gradient").write_text("synthetic gradient boundary\n")
    monkeypatch.setattr(tm, "prepare_control", lambda **kwargs: True)
    monkeypatch.setattr(
        tm, "calc_energy", lambda: (False, None) if step == "energy" else (True, -1.25)
    )
    monkeypatch.setattr(tm, "calc_gradients", lambda: step != "gradient")
    monkeypatch.setattr(
        tm.restraints, "isotropic", lambda *args: (0.0, np.zeros((2, 3)))
    )
    monkeypatch.setattr(
        tm, "rewrite_turbomole_energy_and_gradient_files", lambda *args: None
    )
    monkeypatch.setattr(tm, "update_coord", lambda: False)
    with pytest.raises(error):
        subject.optimize()
    assert not Path(subject.result_xyz_file).exists()


@pytest.mark.parametrize(
    "actual_results,cartesian,expected,calls",
    [
        ([True], True, True, 1),
        ([False, True], True, True, 2),
        ([False, False], True, True, 2),
        ([False, False], False, False, 2),
    ],
)
def test_coordinate_update_preserves_bounded_native_fallback_schedule(
    subject, monkeypatch, actual_results, cartesian, expected, calls
):
    module = Mock(return_value=False)
    monkeypatch.setattr(tm, "run_turbomole_module", module)
    monkeypatch.setattr(tm, "actual", Mock(side_effect=actual_results))
    fallback = Mock(return_value=cartesian)
    monkeypatch.setattr(tm, "statpt_cart", fallback)
    assert tm.update_coord() is expected
    assert module.call_args_list == [(("statpt",),)] * calls
    assert fallback.call_count == (1 if actual_results == [False, False] else 0)


def test_internal_coordinate_generation_without_coord_does_not_launch_define(
    tmp_path, monkeypatch
):
    monkeypatch.chdir(tmp_path)
    execute = Mock()
    monkeypatch.setattr(tm.subp, "check_call", execute)
    assert tm.generate_internal_coordinates() is False
    execute.assert_not_called()


@pytest.mark.parametrize("success_attempt", [1, 3, None])
@pytest.mark.parametrize("process_error", [False, True])
def test_internal_coordinate_generation_preserves_metric_budget(
    subject, monkeypatch, success_attempt, process_error
):
    Path("control").write_text("$redund_inp\n metric 3\n$end\n")
    Path("def_inp").write_text("synthetic define input\n")
    metric = Mock(wraps=tm.change_or_add_metfic_in_redund_inp)
    monkeypatch.setattr(tm, "change_or_add_metfic_in_redund_inp", metric)
    attempts = []

    def execute(command, **kwargs):
        assert command == ["define"]
        attempts.append(command)
        if len(attempts) == success_attempt:
            kwargs["stdout"].write("ATOMIC ATTRIBUTE DATA\n")
        elif process_error:
            raise subprocess.CalledProcessError(1, command)
        else:
            kwargs["stdout"].write("synthetic unsuccessful define output\n")

    monkeypatch.setattr(tm.subp, "check_call", execute)
    assert tm.generate_internal_coordinates() is (success_attempt is not None)
    expected = [3, 2, 1, -1, -2][:success_attempt]
    assert [call.args[0] for call in metric.call_args_list] == expected
    assert len(attempts) == len(expected)
    assert "$redund_inp" not in Path("control").read_text()
