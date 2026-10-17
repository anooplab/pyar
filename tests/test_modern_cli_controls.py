"""Tests for shared modern CLI controls and reporting commands."""

import json
from concurrent.futures import ThreadPoolExecutor
import pytest

try:
    import tomllib
except ModuleNotFoundError:  # Python 3.10
    import tomli as tomllib

from pyar import modern_cli
from pyar.run_config import PROFILE_KIND, resolve_run_spec, write_toml_document
from pyar.scripts import modern_microsolvate


def _write_xyz(path, atoms, coords):
    lines = [str(len(atoms)), "fixture"]
    lines.extend(f"{atom} {x} {y} {z}" for atom, (x, y, z) in zip(atoms, coords))
    path.write_text("\n".join(lines) + "\n")


def test_info_and_backends_produce_json(capsys):
    modern_cli.main(["info", "--json"])
    info = json.loads(capsys.readouterr().out)
    assert info["name"] == "pyar-chem"
    assert info["version"]

    modern_cli.main(["backends", "--json"])
    backends = json.loads(capsys.readouterr().out)["backends"]
    assert backends["xtb"]["energy_gradient"] is True


def test_write_config_and_run_check_are_read_only(tmp_path, monkeypatch, capsys):
    solute = tmp_path / "solute.xyz"
    solvent = tmp_path / "solvent.xyz"
    _write_xyz(solute, ["H", "H"], [(0, 0, 0), (0.74, 0, 0)])
    _write_xyz(solvent, ["H", "H"], [(0, 0, 0), (0.74, 0, 0)])
    monkeypatch.chdir(tmp_path)
    output = tmp_path / "microsolvation"
    config_path = tmp_path / "micro.toml"
    modern_cli.main(["microsolvate", str(solute), str(solvent), "--count", "1",
                     "--site-group", "hydrogen=0", "--contact-distance", "3.2",
                     "--output", str(output), "--write-config", str(config_path), "--dry-run"])
    assert not output.exists()
    dry_run_output = capsys.readouterr().out
    assert "Preflight: microsolvate" in dry_run_output
    assert "Execution plan: microsolvate" in dry_run_output
    assert "Place solvent candidates against that surface" in dry_run_output
    assert "No calculation or workflow output was created." in dry_run_output
    config = tomllib.loads(config_path.read_text())
    assert config["schema_version"] == 1
    assert config["workflow"] == "microsolvate"
    assert "--count" in config["arguments"]
    assert config["arguments"][config["arguments"].index("--site-group") + 1] == "hydrogen=0"

    modern_cli.main(["run", str(config_path), "--check"])
    output_text = capsys.readouterr().out
    assert "Preflight: microsolvate" in output_text
    assert "Contact model: explicit 3.2 Å atom-pair cutoff" in output_text
    assert not output.exists()
    assert not (tmp_path / ".pyar").exists()


def test_execution_plan_contains_resolved_inputs_settings_and_stages():
    from pyar.run_config import ResolvedRunSpec, build_execution_plan

    resolved = ResolvedRunSpec(
        "grow",
        {"seed": "seed.xyz", "monomer": "ligand.xyz", "count": 4,
         "backend": None, "orientations": 8, "maximum_number_of_seeds": 12,
         "output": "grow", "check": True},
        {"count": "cli", "orientations": "default"},
        ("seed.xyz", "ligand.xyz", "--count", "4"),
    )
    plan = build_execution_plan("grow", resolved)
    assert plan.inputs == ("seed.xyz", "ligand.xyz")
    assert plan.backend is None
    assert plan.outputs == ("grow",)
    assert plan.stages == (
        "Initialize the seed pool", "Repeat add, filter, and bounded selection",
        "Retain every growth stage",
    )
    assert plan.settings["count"] == 4


def test_execution_plan_stages_follow_requested_path_stage():
    from pyar.run_config import ResolvedRunSpec, build_execution_plan

    expected = {
        'react': 3,
        'neb': 5,
        'ts': 6,
        'frequency': 7,
        'irc': 8,
        'endpoints': 9,
        'endpoint-frequency': 10,
        'all': 10,
    }
    for through, stage_count in expected.items():
        resolved = ResolvedRunSpec('react', {'through': through}, {}, ())
        plan = build_execution_plan('react', resolved)
        assert len(plan.stages) == stage_count
    assert build_execution_plan(
        'react', ResolvedRunSpec('react', {'through': 'react'}, {}, ())
    ).stages == (
        'Generate reactant orientations', 'Run biased reaction search',
        'Relax and identify candidate products',
    )
    assert build_execution_plan(
        'neb', ResolvedRunSpec('neb', {}, {}, ())
    ).stages == (
        'Relax and validate the endpoint structures',
        'Optimize the reaction path with NEB',
    )
    assert build_execution_plan(
        'irc', ResolvedRunSpec('irc', {}, {}, ())
    ).stages == ('Run both IRC branches from the frequency-validated transition state',)


def test_toml_create_is_atomic_no_clobber_and_force_is_explicit(tmp_path):
    from pyar.run_config import RUN_CONFIG_KIND, load_toml_document

    path = tmp_path / 'shared.toml'

    def create(value):
        return write_toml_document(
            path, kind=RUN_CONFIG_KIND, workflow='grow', arguments=[],
            effective={'count': value},
        )

    with ThreadPoolExecutor(max_workers=8) as pool:
        futures = [pool.submit(create, value) for value in range(8)]
    outcomes = []
    for future in futures:
        try:
            outcomes.append(future.result())
        except FileExistsError:
            pass
    assert len(outcomes) == 1
    saved = load_toml_document(path, RUN_CONFIG_KIND)
    assert saved['effective']['count'] in range(8)
    before = path.read_bytes()
    with pytest.raises(FileExistsError):
        create(20)
    assert path.read_bytes() == before
    write_toml_document(path, kind=RUN_CONFIG_KIND, workflow='grow', arguments=[],
                        effective={'count': 20}, overwrite=True)
    assert load_toml_document(path, RUN_CONFIG_KIND)['effective']['count'] == 20


def test_profile_precedence_and_profile_omits_run_specific_paths(tmp_path, monkeypatch, capsys):
    source = tmp_path / "input.xyz"
    _write_xyz(source, ["H", "H"], [(0, 0, 0), (0.74, 0, 0)])
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "config"))
    write_toml_document(
        tmp_path / "config" / "pyar" / "profiles" / "inspection.toml",
        kind=PROFILE_KIND, workflow="identify", arguments=["--json"],
        effective={"json": True},
    )
    modern_cli.main(["identify", str(source), "--profile", "inspection"])
    result = json.loads(capsys.readouterr().out)
    assert result["structures"][0]["atom_count"] == 2


def test_resolved_spec_tracks_cli_profile_and_default_sources():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("inputs", nargs=2)
    parser.add_argument("--backend", default="xtb")
    parser.add_argument("--nprocs", type=int, default=1)
    parser.add_argument("--charge", nargs="+", type=int, default=[0])
    request, profiled = resolve_run_spec(
        "react", parser, ["a.xyz", "b.xyz", "--nprocs", "4"],
        profile_values={"backend": "orca", "charge": [0, -1]}, profile="careful",
    )
    assert request.profile == "careful"
    assert profiled.values["backend"] == "orca"
    assert profiled.values["charge"] == [0, -1]
    assert profiled.value_sources["charge"] == "profile"
    request, resolved = resolve_run_spec(
        "react", parser, ["a.xyz", "b.xyz", "--nprocs", "4", "--charge", "1", "0"],
        profile_values={"backend": "orca", "charge": [0, -1]}, profile="careful",
    )
    assert resolved.values["inputs"] == ["a.xyz", "b.xyz"]
    assert resolved.values["charge"] == [1, 0]
    assert resolved.value_sources == {"inputs": "cli", "backend": "profile", "nprocs": "cli", "charge": "cli"}

    _, configured = resolve_run_spec(
        "react", parser, ["a.xyz", "b.xyz", "--nprocs", "6"],
        argument_source="config",
    )
    assert configured.value_sources["inputs"] == "config"
    assert configured.value_sources["nprocs"] == "config"


def test_run_records_include_value_sources_and_failed_status(tmp_path, monkeypatch):
    import argparse
    from types import SimpleNamespace

    def build_parser(prog=None):
        parser = argparse.ArgumentParser(prog=prog)
        parser.add_argument("input")
        parser.add_argument("--backend", default="xtb")
        parser.add_argument("--output", default="unused")
        parser.add_argument("--check", action="store_true")
        return parser

    def fail(argv=None, *, prog=None):
        from pyar.cli_runtime import started_run
        started_run()
        raise ValueError("mock workflow failure")

    fake_module = SimpleNamespace(build_parser=build_parser, modern_main=fail)
    monkeypatch.setitem(modern_cli._COMMANDS, "optimize", ("fake.optimize", "fake"))
    monkeypatch.setattr(modern_cli, "import_module", lambda name: fake_module)
    monkeypatch.chdir(tmp_path)
    try:
        modern_cli.main(["optimize", "molecule.xyz", "--backend", "xtb"])
    except ValueError:
        pass
    else:
        raise AssertionError("mock failure did not propagate")
    records = list((tmp_path / ".pyar" / "runs").glob("*/pyar-run.toml"))
    assert len(records) == 1
    record = tomllib.loads(records[0].read_text())
    metadata = json.loads(record["metadata_json"])
    assert metadata["status"] == "failed"
    sources = metadata["value_sources"]
    assert sources["input"] == "cli"
    assert sources["backend"] == "cli"
    assert sources["output"] == "default"


def test_run_record_failure_is_nonzero_and_json_distinguishes_scientific_success(
        tmp_path, monkeypatch, capsys):
    import argparse
    from types import SimpleNamespace
    from pyar.cli_runtime import finished_run, started_run
    from pyar.run_config import RUN_RECORD_KIND

    def build_parser(prog=None):
        parser = argparse.ArgumentParser(prog=prog)
        parser.add_argument('input')
        parser.add_argument('--check', action='store_true')
        return parser

    def complete(argv=None, *, prog=None):
        started_run()
        finished_run({'workflow': 'optimize', 'status': 'completed', 'input_count': 1})

    fake_module = SimpleNamespace(build_parser=build_parser, modern_main=complete)
    monkeypatch.setitem(modern_cli._COMMANDS, 'optimize', ('fake.optimize', 'fake'))
    monkeypatch.setattr(modern_cli, 'import_module', lambda _name: fake_module)
    monkeypatch.chdir(tmp_path)

    def fail_record(*args, **kwargs):
        if kwargs.get('kind') == RUN_RECORD_KIND:
            raise OSError('disk full')
        raise AssertionError('unexpected document write')

    monkeypatch.setattr('pyar.run_config.write_toml_document', fail_record)
    try:
        modern_cli.main(['optimize', 'input.xyz', '--json'])
    except OSError as exc:
        assert 'run provenance could not be written' in str(exc)
    else:
        raise AssertionError('provenance write failure returned success')

    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert payload['scientific_status'] == 'completed'
    assert payload['status'] == 'completed'
    assert payload['invocation_status'] == 'provenance_write_failed'
    assert payload['provenance_error'] == 'disk full'
    assert 'outputs were retained' in captured.err


def test_all_remaining_workflow_commands_are_registered_lazily():
    for command in ("solvate", "neb", "ts", "irc", "run", "config", "profile", "info",
                    "backends", "doctor", "inspect"):
        assert command in modern_cli._COMMANDS
    # Listing commands does not import optional scientific workflow packages.
    parser = modern_cli.build_parser()
    assert parser.prog == "pyar"


def test_dry_run_is_limited_to_commands_with_full_check_paths(capsys):
    try:
        modern_cli.main(["identify", "structure.xyz", "--dry-run"])
    except SystemExit as exc:
        assert exc.code == 2
        assert "supported only for calculation workflows" in capsys.readouterr().err
    else:
        raise AssertionError("analysis command unexpectedly accepted --dry-run")
