"""Read-only CLI settings and report boundaries, independent of subprocess tracing."""

from copy import deepcopy
import json
import os
from unittest.mock import Mock

import pytest

from pyar.run_config import RUN_RECORD_KIND, write_toml_document
from pyar.scripts import cli_meta, path_stage


def path_arguments(*options):
    args = path_stage.build_parser("neb").parse_args(
        ["a.xyz", "b.xyz", "--ts-guess", "guess.xyz", *options]
    )
    args.command = "neb"
    return args


@pytest.mark.parametrize(
    "options,match",
    [
        (["--backend", "orca", "--xtb-model", "gfn2"], "--xtb-model requires"),
        (["--backend", "psi4"], "Cartesian gradients"),
        (["--backend", "gaussian"], "explicit --method"),
        (["--backend", "gaussian", "--method", "BP86"], "explicit --basis"),
        (["--backend", "orca", "--method", "BP86"], "DFT methods require --basis"),
        (
            ["--backend", "orca", "--method", "GFN2-xTB", "--basis", "def2-SVP"],
            "defines its own basis",
        ),
        (
            ["--backend", "orca", "--method", "r2scan-3c", "--basis", "def2-SVP"],
            "defines its own basis",
        ),
        (
            ["--backend", "orca", "--method", "r2scan-3c"],
            "require --geometry-optimizer native",
        ),
        (["--backend", "xtb", "--nprocs", "0"], "must be positive"),
    ],
)
def test_invalid_path_settings_fail_read_only_before_execution(
    options, match, tmp_path, monkeypatch
):
    monkeypatch.chdir(tmp_path)
    engine = Mock(side_effect=AssertionError("must not launch"))
    monkeypatch.setattr(path_stage, "run_neb", engine)
    args = path_arguments(*options)
    before = deepcopy(vars(args))
    with pytest.raises(ValueError, match=match):
        path_stage._resolve_settings(args)
    assert vars(args) == before
    assert list(tmp_path.iterdir()) == []
    engine.assert_not_called()


@pytest.mark.parametrize(
    "options,backend,method,basis,model",
    [
        (["--backend", "xtb"], "xtb", None, None, "gfn2"),
        (["--backend", "xtb", "--xtb-model", "gxtb"], "xtb", None, None, "gxtb"),
        (
            ["--backend", "orca16", "--method", "BP86", "--basis", "def2-SVP"],
            "orca",
            "BP86",
            "def2-SVP",
            None,
        ),
        (
            ["--backend", "g16", "--method", "B3LYP", "--basis", "6-31G"],
            "gaussian",
            "B3LYP",
            "6-31G",
            None,
        ),
        (["--backend", "orca", "--method", "GFN2-xTB"], "orca", "GFN2-xTB", None, None),
    ],
)
def test_path_settings_preserve_owner_method_and_optimizer(
    options, backend, method, basis, model
):
    args = path_arguments(*options)
    before = deepcopy(vars(args))
    settings = path_stage._resolve_settings(args)
    assert settings["software"] == backend
    assert (settings["method"], settings["basis"]) == (method, basis)
    assert settings.get("xtb_model") == model
    assert settings["geometry_optimizer"] == "geometric"
    assert settings["opt_target"] == "minimum"
    assert settings["gamma"] is None
    assert vars(args) == before


@pytest.mark.parametrize("as_json", [False, True])
def test_inspect_existing_json_is_read_only_and_reports_resolved_request(
    tmp_path, capsys, as_json
):
    payload = {
        "workflow": "neb",
        "status": "complete",
        "metadata": {"status": "ignored fallback", "created_utc": "recorded time"},
        "effective": {"software": "xtb"},
        "state": {"stage": "neb"},
    }
    record = tmp_path / "state.json"
    record.write_text(json.dumps(payload))
    before = record.read_bytes()
    cli_meta.inspect_run([str(tmp_path), *(["--json"] if as_json else [])])
    output = capsys.readouterr()
    assert output.err == ""
    if as_json:
        assert json.loads(output.out) == payload
    else:
        assert "Status: complete" in output.out
        assert "ignored fallback" not in output.out
        assert "Created Utc: recorded time" in output.out
        assert "Resolved request:\n  software: xtb" in output.out
        assert "Workflow state:\n  stage: neb" in output.out
    assert record.read_bytes() == before
    assert list(tmp_path.iterdir()) == [record]


def test_inspect_prefers_latest_immutable_invocation_record(tmp_path, capsys):
    for index in (1, 2):
        record = write_toml_document(
            tmp_path / f"pyar-run-{index}.toml",
            kind=RUN_RECORD_KIND,
            workflow="neb",
            arguments=[],
            metadata={"status": f"attempt-{index}"},
        )
        os.utime(record, ns=(index * 1_000_000_000, index * 1_000_000_000))
    (tmp_path / "state.json").write_text('{"status":"mutable state"}')
    before = {p.name: p.read_bytes() for p in tmp_path.iterdir()}
    cli_meta.inspect_run([str(tmp_path), "--json"])
    assert json.loads(capsys.readouterr().out)["metadata"]["status"] == "attempt-2"
    assert {p.name: p.read_bytes() for p in tmp_path.iterdir()} == before


@pytest.mark.parametrize(
    "kind", ["empty-directory", "missing-file", "invalid-json", "invalid-toml"]
)
def test_inspect_invalid_record_has_parser_error_and_no_writes(tmp_path, capsys, kind):
    path = tmp_path
    if kind != "empty-directory":
        path = tmp_path / ("broken.toml" if kind == "invalid-toml" else "broken.json")
        if kind != "missing-file":
            path.write_text("not a valid report")
    before = {p.name: p.read_bytes() for p in tmp_path.iterdir()}
    with pytest.raises(SystemExit) as caught:
        cli_meta.inspect_run([str(path)])
    assert caught.value.code == 2
    output = capsys.readouterr()
    assert output.out == ""
    assert "pyar inspect: error:" in output.err
    assert {p.name: p.read_bytes() for p in tmp_path.iterdir()} == before


@pytest.mark.parametrize("as_json", [False, True])
def test_backend_inspection_reports_missing_dependencies_without_importing_them(
    monkeypatch, capsys, as_json
):
    monkeypatch.setattr(cli_meta.shutil, "which", lambda name: None)
    monkeypatch.setattr(cli_meta, "_module_available", lambda name: False)
    cli_meta.backends(["aiqm1_mlatom", *(["--json"] if as_json else [])])
    output = capsys.readouterr()
    assert output.err == ""
    if as_json:
        record = json.loads(output.out)["backends"]["aiqm1_mlatom"]
        assert record["availability"] == {
            "executables": {"g16": None},
            "python_modules": {"mlatom": False},
        }
        assert record["optional_extra"] == "ml"
    else:
        assert "g16: not found on PATH" in output.out
        assert "mlatom: not installed" in output.out
        assert "pyar-chem[ml]" in output.out


def test_unknown_backend_is_an_actionable_parser_error(capsys):
    with pytest.raises(SystemExit) as caught:
        cli_meta.backends(["not-a-backend"])
    assert caught.value.code == 2
    assert "Unknown backend" in capsys.readouterr().err


@pytest.mark.parametrize("as_json", [False, True])
def test_backend_listing_needs_no_local_dependency_probes(monkeypatch, capsys, as_json):
    probe = Mock(side_effect=AssertionError("listing must not probe dependencies"))
    monkeypatch.setattr(cli_meta.shutil, "which", probe)
    monkeypatch.setattr(cli_meta, "_module_available", probe)
    cli_meta.backends(["--json"] if as_json else [])
    output = capsys.readouterr().out
    if as_json:
        records = json.loads(output)["backends"]
        assert records["xtb"]["energy_gradient"] is True
        assert records["psi4"]["energy_gradient"] is False
        assert all("availability" not in record for record in records.values())
    else:
        assert "xtb (semiempirical)" in output
        assert "psi4 (dft_qc)" in output
    probe.assert_not_called()


@pytest.mark.parametrize("error", [ImportError, ValueError, AttributeError])
def test_broken_optional_module_spec_is_reported_unavailable(monkeypatch, error):
    def find_spec(name):
        raise error("missing or malformed module spec")

    monkeypatch.setattr(cli_meta.importlib.util, "find_spec", find_spec)
    assert cli_meta._module_available("optional-backend") is False
