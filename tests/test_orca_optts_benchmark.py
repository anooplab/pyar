from pyar.benchmarks.orca_optts import (
    build_orca_optts_input,
    parse_orca_optts_output,
)


def test_orca_optts_input_selects_gfn2_and_records_optimizer_policy():
    text = build_orca_optts_input(
        ["H", "H"], [[0, 0, 0], [0, 0, 0.8]], charge=0, multiplicity=1,
    )
    assert "! GFN2-xTB OptTS" in text
    assert "CoordSys redundant" in text
    assert "InHess XTB2" in text
    assert "Update Bofill" in text
    assert "MaxIter 200" in text


def test_orca_input_matches_manifest_gradient_threshold_in_atomic_units():
    text = build_orca_optts_input(
        ["H", "H"], [[0, 0, 0], [0, 0, 0.8]], charge=0, multiplicity=1,
        ts_fmax_ev_per_angstrom=0.02,
    )
    assert "TolMaxG 0.000388" in text
    assert "TolRMSG 0.000259" in text


def test_orca_parser_counts_optimizer_gradient_events_and_convergence():
    output = """
MaxIter 200
* GEOMETRY OPTIMIZATION CYCLE 1 *
Time for energy+gradient : 0.02 s
* GEOMETRY OPTIMIZATION CYCLE 2 *
Time for energy+gradient : 0.01 s
***        THE OPTIMIZATION HAS CONVERGED     ***
****ORCA TERMINATED NORMALLY****
"""
    parsed = parse_orca_optts_output(output)
    assert parsed["optimizer_converged"] is True
    assert parsed["failure_class"] is None
    assert parsed["optimizer_steps"] == 2
    assert parsed["orca_gradient_events"] == 2


def test_orca_parser_identifies_cycle_limit_only_from_failure_message():
    output = """
MaxIter 200
THE OPTIMIZATION DID NOT CONVERGE: maximum number of geometry optimization cycles reached
****ORCA TERMINATED NORMALLY****
"""
    parsed = parse_orca_optts_output(output)
    assert parsed["optimizer_converged"] is False
    assert parsed["failure_class"] == "max_iterations"


def test_orca_rejects_different_xtb_executable_for_validation(tmp_path, monkeypatch):
    import pytest
    from types import SimpleNamespace
    from pyar.benchmarks import orca_optts as runner
    from pyar.benchmarks.ts_optimizer import TSOptimizerBenchmarkError
    first, second = tmp_path / 'xtb-one', tmp_path / 'xtb-two'
    first.touch()
    second.touch()
    monkeypatch.setattr(runner.shutil, 'which', lambda _: str(first))
    spec = SimpleNamespace(qc_model={'xtb_model': 'gfn2'}, cases=[SimpleNamespace(id='case')])
    with pytest.raises(TSOptimizerBenchmarkError, match='must match xtb on PATH'):
        runner.run_orca_optts_case(spec, case_id='case', output=tmp_path / 'run',
                                  xtb_executable=second, orca_executable=first)
    assert not (tmp_path / 'run').exists()


def test_native_optts_runner_preserves_command_and_captured_process_output(tmp_path, monkeypatch):
    from types import SimpleNamespace
    from pyar.benchmarks import orca_optts as runner

    observed = {}

    def fake_run(command, **kwargs):
        observed.update(command=command, kwargs=kwargs)
        return SimpleNamespace(stdout="stdout", stderr="stderr", returncode=7)

    monkeypatch.setattr(runner.subprocess, "run", fake_run)
    output, returncode, elapsed = runner._execute_native_optts(
        "/opt/orca", tmp_path / "job.inp", tmp_path, {"XTBEXE": "/opt/xtb"},
    )

    assert observed["command"] == ["/opt/orca", str(tmp_path / "job.inp")]
    assert observed["kwargs"] == {
        "cwd": tmp_path, "env": {"XTBEXE": "/opt/xtb"},
        "capture_output": True, "text": True, "check": False,
    }
    assert output == "stdout\nstderr"
    assert returncode == 7
    assert elapsed >= 0


def test_native_optts_runner_maps_oserror_to_existing_failure_evidence(tmp_path, monkeypatch):
    from pyar.benchmarks import orca_optts as runner

    def fail_to_start(*_args, **_kwargs):
        raise OSError("executable disappeared")

    monkeypatch.setattr(runner.subprocess, "run", fail_to_start)
    output, returncode, elapsed = runner._execute_native_optts(
        "/opt/orca", tmp_path / "job.inp", tmp_path, {},
    )

    assert output == "OSError: executable disappeared"
    assert returncode == 127
    assert elapsed >= 0
