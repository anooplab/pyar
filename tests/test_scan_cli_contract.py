"""Scan adapter validation order and execution-publication boundaries."""

from pathlib import Path

import pytest

from pyar.backend_errors import BackendExecutionError
from pyar.persistence import LockError
from pyar.scripts import scan_bond
import pyar.scan_request as shared
import pyar.workflows.scan_bond as workflow


@pytest.fixture
def inputs(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    for name in ("a.xyz", "b.xyz"):
        Path(name).write_text("2\nfixture\nH 0 0 0\nH 0.74 0 0\n")
    return ["a.xyz", "b.xyz", "--atoms", "0", "0"]


@pytest.mark.parametrize(
    "intra, options, message",
    [
        (True, ["--orientations", "2"], "Intramolecular scans preserve"),
        (
            False,
            ["--xtb-model", "gxtb", "--orientations", "0"],
            "A backend is required.",
        ),
        (
            False,
            ["--backend", "orca", "--xtb-model", "gfn2", "--nprocs", "0"],
            "--xtb-model requires --backend xtb",
        ),
        (
            False,
            ["--backend", "orca", "--nprocs", "0"],
            "nprocs, opt-cycles and scf-cycles must be positive",
        ),
        (
            False,
            ["--backend", "orca", "--through", "all", "--images", "4"],
            "--basis is required for ORCA DFT methods",
        ),
        (
            False,
            [
                "--backend",
                "xtb",
                "--characterize-unique-products",
                "--interpolation",
                "geodesic",
            ],
            "--characterize-unique-products requires --through neb",
        ),
        (
            True,
            ["--backend", "xtb", "--through", "neb", "--interpolation", "geodesic"],
            "currently supports linear or idpp",
        ),
        (
            False,
            ["--backend", "xtb", "--through", "all", "--images", "4"],
            "images must be an odd integer",
        ),
    ],
)
def test_option_error_precedes_missing_input_and_never_executes(
    inputs, monkeypatch, capsys, intra, options, message
):
    Path("a.xyz").unlink()
    monkeypatch.setattr(
        scan_bond,
        "run_scan_bond_request",
        lambda *a, **k: pytest.fail("executed invalid request"),
    )
    base = ["a.xyz", "--atoms", "0", "1"] if intra else inputs
    with pytest.raises(SystemExit) as caught:
        scan_bond.main([*base, *options], prog="pyar scan-bond")
    assert caught.value.code == 2
    assert message in capsys.readouterr().err


@pytest.fixture
def lifecycle(monkeypatch):
    events = []
    original = workflow.scan_request_fields

    def fields(*args, **kwargs):
        events.append("fields")
        return original(*args, **kwargs)

    monkeypatch.setattr(workflow, "scan_request_fields", fields)
    monkeypatch.setattr(
        workflow, "validate_scan_restart", lambda *a, **k: events.append("restart")
    )
    monkeypatch.setattr(
        shared, "preflight_scan", lambda *a: events.append("preflight") or ["xtb"]
    )
    monkeypatch.setattr(
        scan_bond, "prepared_run", lambda *a, **k: events.append("prepared")
    )
    monkeypatch.setattr(scan_bond, "started_run", lambda: events.append("started"))
    monkeypatch.setattr(
        scan_bond, "finished_run", lambda result: events.append(("finished", result))
    )
    return events


@pytest.mark.parametrize("check", [False, True])
@pytest.mark.parametrize("recover", [False, True])
def test_restart_precedes_preflight_and_check_never_starts(
    inputs, lifecycle, monkeypatch, check, recover
):
    calls = []
    result = {"status": "complete", "output_dir": "scan_bond"}

    def engine(request, **kwargs):
        lifecycle.append("execute")
        calls.append((request, kwargs))
        return result

    monkeypatch.setattr(scan_bond, "run_scan_bond_request", engine)
    flags = (["--check"] if check else []) + (
        ["--recover-stale-lock"] if recover else []
    )
    returned = scan_bond.main([*inputs, "--backend", "xtb", *flags])
    expected = ["fields", "restart", "preflight", "prepared"]
    if check:
        assert returned is None and calls == []
    else:
        expected += ["started", "execute", ("finished", result)]
        assert returned is result
        assert calls[0][1] == ({"recover_stale": True} if recover else {})
        assert calls[0][0]["qc_params"]["xtb_model"] == "gfn2"
    assert lifecycle == expected


@pytest.mark.parametrize(
    "failure",
    [
        ValueError("changed checkpoint"),
        ImportError("missing API"),
        OSError("unreadable checkpoint"),
        RuntimeError("capability unavailable"),
    ],
)
def test_restart_failure_is_wrapped_before_preflight_or_preparation(
    inputs, lifecycle, monkeypatch, capsys, failure
):
    def reject(*args):
        lifecycle.append("reject")
        raise failure

    monkeypatch.setattr(workflow, "validate_scan_restart", reject)
    with pytest.raises(SystemExit) as caught:
        scan_bond.main([*inputs, "--backend", "xtb"])
    assert caught.value.code == 2
    assert lifecycle == ["fields", "reject"]
    error = capsys.readouterr().err
    assert (
        f"Preflight failed for: scan-bond\n{failure}\nNo calculations were started."
        in error
    )


@pytest.mark.parametrize(
    "failure",
    [
        ValueError("backend failed"),
        BackendExecutionError("native backend failed"),
        KeyboardInterrupt(),
        SystemExit(9),
    ],
)
def test_execution_errors_keep_their_identity_and_do_not_publish_finished(
    inputs, lifecycle, monkeypatch, failure
):
    def reject(*args, **kwargs):
        raise failure

    monkeypatch.setattr(scan_bond, "run_scan_bond_request", reject)
    with pytest.raises(type(failure)) as caught:
        scan_bond.main([*inputs, "--backend", "xtb"])
    assert caught.value is failure
    assert lifecycle == ["fields", "restart", "preflight", "prepared", "started"]


def test_lock_error_is_cli_error_after_started(inputs, lifecycle, monkeypatch, capsys):
    def locked(*args, **kwargs):
        raise LockError("active writer")

    monkeypatch.setattr(scan_bond, "run_scan_bond_request", locked)
    with pytest.raises(SystemExit) as caught:
        scan_bond.main([*inputs, "--backend", "xtb"])
    assert caught.value.code == 2
    assert "active writer" in capsys.readouterr().err
    assert lifecycle[-1] == "started"


@pytest.mark.parametrize(
    "status", ["failed", "partial_failed", "scientific_gate_failed"]
)
def test_noncomplete_result_is_published_before_nonzero_exit(
    inputs, lifecycle, monkeypatch, capsys, status
):
    result = {"status": status, "output_dir": "scan_bond"}
    monkeypatch.setattr(scan_bond, "run_scan_bond_request", lambda *a, **k: result)
    with pytest.raises(SystemExit) as caught:
        scan_bond.main([*inputs, "--backend", "xtb"])
    assert caught.value.code == 1
    assert lifecycle[-1] == ("finished", result)
    streams = capsys.readouterr()
    assert f"scan-bond: {status} through scan; summaries: scan_bond" in streams.out
    assert (
        streams.err
        == "One or more orientations failed; inspect summary.json for the failed stage and gate.\n"
    )


def test_unexpected_preflight_exception_is_not_reclassified(
    inputs, lifecycle, monkeypatch
):
    failure = TypeError("unexpected preflight contract failure")

    def reject(*args):
        raise failure

    monkeypatch.setattr(shared, "preflight_scan", reject)
    with pytest.raises(TypeError) as caught:
        scan_bond.main([*inputs, "--backend", "xtb"])
    assert caught.value is failure
    assert lifecycle == ["fields", "restart"]
