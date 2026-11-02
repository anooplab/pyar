"""Observable dispatch contracts, independent of dispatcher helper structure."""

import argparse
import json
import logging
from types import SimpleNamespace

import pytest

from pyar import cli
from pyar.cli_runtime import finished_run, observe_run, prepared_run, started_run


@pytest.fixture
def task(monkeypatch):
    def install(entry, *, parser=True):
        def build_parser(prog=None):
            result = argparse.ArgumentParser(prog=prog)
            result.add_argument("--count", type=int, default=1)
            result.add_argument("--check", action="store_true")
            return result

        module = SimpleNamespace(main=entry)
        if parser:
            module.build_parser = build_parser
        monkeypatch.setattr(cli, "import_module", lambda name: module)
        return module

    return install


def prepare(count):
    prepared_run(
        SimpleNamespace(count=count, check=False),
        request={"count": count},
        backend={},
        molecules=[],
        requirements=[],
        outputs=[],
    )


@pytest.mark.parametrize(
    "failure", [ValueError("backend failed"), SystemExit(7), KeyboardInterrupt()]
)
@pytest.mark.parametrize("record_fails", [False, True])
def test_execution_failure_identity_survives_provenance_attempt(
    task, monkeypatch, capsys, failure, record_fails
):
    events = []

    def entry(*args, **kwargs):
        started_run()
        finished_run({"status": "scientific_gate_failed"})
        events.append("execute")
        raise failure

    def record(*args, **kwargs):
        assert kwargs["failure"] is failure
        events.append("record")
        if record_fails:
            raise OSError("record disk full")

    task(entry)
    monkeypatch.setattr(cli, "_write_execution_record", record)
    with pytest.raises(type(failure)) as caught:
        cli.invoke("grow", ["--json"])
    assert caught.value is failure
    assert events == ["execute", "record"]
    payload = json.loads(capsys.readouterr().out)
    assert payload["status"] == payload["scientific_status"] == "scientific_gate_failed"
    assert payload["invocation_status"] == "failed"
    assert payload.get("provenance_error") == (
        "record disk full" if record_fails else None
    )


def test_repeated_preparation_freezes_first_config_but_records_latest_request(
    tmp_path, task, monkeypatch
):
    from pyar.run_config import load_toml_document

    config = tmp_path / "run.toml"
    events = []

    def entry(*args, **kwargs):
        prepare(2)
        events.append(
            load_toml_document(config, "pyar-run-config")["effective"]["count"]
        )
        prepare(3)
        started_run()
        return {"status": "complete"}

    def record(command, request, resolved, observation, **kwargs):
        events.append(
            (
                resolved.values["count"],
                resolved.value_sources["count"],
                observation.result,
            )
        )

    task(entry)
    monkeypatch.setattr(cli, "_write_execution_record", record)
    result = cli.invoke("grow", ["--write-config", str(config)])
    assert result == {"status": "complete"}
    assert events == [2, (3, "workflow-default", result)]
    assert load_toml_document(config, "pyar-run-config")["effective"]["count"] == 2


@pytest.mark.parametrize("started", [False, True])
def test_prepared_config_failure_records_only_if_execution_started(
    tmp_path, task, monkeypatch, capsys, started
):
    failure = OSError("config disk full")
    records = []

    def entry(*args, **kwargs):
        if started:
            started_run()
        prepare(2)
        pytest.fail("must not execute after config publication failed")

    def write(*args, **kwargs):
        raise failure

    task(entry)
    monkeypatch.setattr("pyar.run_config.write_toml_document", write)
    monkeypatch.setattr(
        cli, "_write_execution_record", lambda *a, **k: records.append(k["failure"])
    )
    with pytest.raises(OSError) as caught:
        cli.invoke("grow", ["--json", "--write-config", str(tmp_path / "run.toml")])
    assert caught.value is failure
    assert records == ([failure] if started else [])
    assert json.loads(capsys.readouterr().out)["invocation_status"] == "failed"


@pytest.mark.parametrize(
    "arguments", [["--quiet"], ["--json"], ["--dry-run", "--json"]]
)
def test_stream_routing_and_logger_restoration(task, monkeypatch, capsys, arguments):
    logger = logging.getLogger("pyar")
    previous = logger.level, list(logger.handlers)
    recorded = []

    def entry(argv, **kwargs):
        print("workflow progress")
        if "--check" not in argv:
            started_run()
        return {"status": "complete"}

    task(entry)
    monkeypatch.setattr(
        cli, "_write_execution_record", lambda *a, **k: recorded.append(True)
    )
    cli.invoke("grow", arguments)
    streams = capsys.readouterr()
    if "--quiet" in arguments:
        assert streams.out == streams.err == ""
    else:
        payload = json.loads(streams.out)
        assert payload["mode"] == ("dry-run" if "--dry-run" in arguments else "run")
        assert streams.err == (
            "" if "--dry-run" in arguments else "workflow progress\n"
        )
    assert recorded == ([] if "--dry-run" in arguments else [True])
    assert (logger.level, logger.handlers) == previous


def test_help_precedes_profile_lookup_and_existing_config_rejection(
    tmp_path, task, capsys
):
    existing = tmp_path / "saved.toml"
    existing.write_bytes(b"preserve")
    task(lambda *a, **k: pytest.fail("help dispatched task"))
    with pytest.raises(SystemExit) as caught:
        cli.invoke(
            "grow",
            ["--help", "--profile", "does-not-exist", "--write-config", str(existing)],
        )
    assert caught.value.code == 0
    assert "Common execution controls" in capsys.readouterr().out
    assert existing.read_bytes() == b"preserve"


def test_nested_observer_is_restored_and_preflight_has_no_run_record(task, monkeypatch):
    records = []
    task(lambda *a, **k: prepare(2))
    monkeypatch.setattr(
        cli, "_write_execution_record", lambda *a, **k: records.append(True)
    )
    with observe_run() as outer:
        prepare(7)
        cli.invoke("grow", ["--check"])
        assert outer.options["count"] == 7 and not outer.started
        prepare(9)
        assert outer.options["count"] == 9
    assert records == []


def test_unhandled_provenance_exception_propagates_before_json(
    task, monkeypatch, capsys
):
    failure = TypeError("unexpected writer failure")

    def entry(*args, **kwargs):
        started_run()
        return {"status": "complete"}

    def record(*args, **kwargs):
        raise failure

    task(entry)
    monkeypatch.setattr(cli, "_write_execution_record", record)
    with pytest.raises(TypeError) as caught:
        cli.invoke("grow", ["--json"])
    assert caught.value is failure
    assert capsys.readouterr().out == ""


def test_source_overrides_and_profile_merge_preserve_caller_data(
    task, monkeypatch, capsys
):
    profile = {"count": 2}
    overrides = {"count": "explicit-replay-override"}
    calls = []

    def entry(argv, **kwargs):
        calls.append(list(argv))
        started_run()
        return {"status": "complete"}

    task(entry)
    monkeypatch.setattr(cli, "_write_execution_record", lambda *a, **k: None)
    result = cli.invoke(
        "grow",
        ["--count", "4", "--json"],
        config_file="source.toml",
        profile_values=profile,
        source_overrides=overrides,
    )
    payload = json.loads(capsys.readouterr().out)
    assert result == {"status": "complete"}
    assert payload["resolved"]["values"]["count"] == 4
    assert payload["resolved"]["value_sources"]["count"] == "explicit-replay-override"
    assert calls == [["--count", "4"]]
    assert profile == {"count": 2} and overrides == {
        "count": "explicit-replay-override"
    }


def test_config_failure_after_return_retains_scientific_result(
    task, monkeypatch, tmp_path, capsys
):
    failure = OSError("config publication failed")
    recorded = []

    def entry(*args, **kwargs):
        started_run()
        return {"status": "complete", "energy": -1.0}

    def write(*args, **kwargs):
        raise failure

    def record(command, request, resolved, observation, **kwargs):
        recorded.append((observation.result, kwargs["failure"]))

    task(entry)
    monkeypatch.setattr("pyar.run_config.write_toml_document", write)
    monkeypatch.setattr(cli, "_write_execution_record", record)
    with pytest.raises(OSError) as caught:
        cli.invoke("grow", ["--json", "--write-config", str(tmp_path / "late.toml")])
    assert caught.value is failure
    payload = json.loads(capsys.readouterr().out)
    assert payload["scientific_status"] == "complete"
    assert payload["invocation_status"] == "failed"
    assert recorded == [({"status": "complete", "energy": -1.0}, failure)]


def test_nested_invocations_have_independent_config_and_provenance(
    task, monkeypatch, tmp_path
):
    from pyar.run_config import load_toml_document

    records = []
    outer = tmp_path / "outer.toml"
    inner = tmp_path / "inner.toml"

    def entry(argv, **kwargs):
        count = int(argv[argv.index("--count") + 1])
        prepare(count)
        started_run()
        if count == 2:
            cli.invoke("grow", ["--count", "3", "--write-config", str(inner)])
        return {"status": "complete", "count": count}

    def record(command, request, resolved, observation, **kwargs):
        records.append((resolved.values["count"], observation.result["count"]))

    task(entry)
    monkeypatch.setattr(cli, "_write_execution_record", record)
    assert (
        cli.invoke("grow", ["--count", "2", "--write-config", str(outer)])["count"] == 2
    )
    assert records == [(3, 3), (2, 2)]
    assert load_toml_document(inner, "pyar-run-config")["effective"]["count"] == 3
    assert load_toml_document(outer, "pyar-run-config")["effective"]["count"] == 2


def test_parserless_task_retains_its_arguments_and_result(task, monkeypatch):
    calls = []

    def entry(argv, **kwargs):
        calls.append((argv, kwargs["prog"]))
        return "analysis-result"

    task(entry, parser=False)
    assert cli.invoke("identify", ["--", "--profile", "literal"]) == "analysis-result"
    assert calls == [(["--", "--profile", "literal"], "pyar identify")]


def test_unknown_task_precedes_control_errors_and_task_import(monkeypatch, capsys):
    monkeypatch.setattr(
        cli, "import_module", lambda name: pytest.fail("imported unknown task")
    )
    with pytest.raises(SystemExit) as caught:
        cli.invoke("not-a-task", ["--quiet", "-v"])
    assert caught.value.code == 2
    assert "Unknown task: not-a-task" in capsys.readouterr().err
