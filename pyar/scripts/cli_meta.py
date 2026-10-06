"""Modern CLI environment, backend, and run-inspection commands."""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path
import shutil
import sys
from datetime import datetime, timezone

from pyar import __version__
from pyar.backend_capabilities import BACKEND_CAPABILITIES
from pyar.run_config import (
    PROFILE_KIND, RUN_CONFIG_KIND, RUN_RECORD_KIND, arguments_from_profile,
    load_toml_document, profile_directory, profile_path, profile_overrides_from_record,
    write_toml_document, arguments_from_values, without_execution_controls,
    explicit_destinations, discard_overridden,
)


def _emit(value, as_json):
    if as_json:
        print(json.dumps(value, indent=2, sort_keys=True))
    else:
        for key, item in value.items():
            print(f"{key}: {item}")


def info(argv=None, *, prog="pyar info"):
    parser = argparse.ArgumentParser(prog=prog, description="Show PyAR installation information")
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args(argv)
    from pyar.run_config import profile_directory
    _emit({"name": "pyar-chem", "version": __version__, "python": sys.version.split()[0],
           "executable": sys.executable, "package": str(Path(__file__).resolve().parents[1]),
           "profile_directory": str(profile_directory())}, args.json)


def backends(argv=None, *, prog="pyar backends"):
    parser = argparse.ArgumentParser(prog=prog, description="List registered PyAR backends and capabilities")
    parser.add_argument("--json", action="store_true")
    parser.add_argument('backend', nargs='?', help='Inspect one backend, including local requirement availability')
    args = parser.parse_args(argv)
    names = sorted(BACKEND_CAPABILITIES)
    if args.backend:
        from pyar.backend_capabilities import normalize_backend_name
        name = normalize_backend_name(args.backend.lower())
        if name not in BACKEND_CAPABILITIES:
            parser.error(f'Unknown backend {args.backend!r}; use pyar backends to list available names')
        names = [name]
    records = {}
    for name in names:
        cap = BACKEND_CAPABILITIES[name]
        records[name] = {
            "family": cap.family,
            "aliases": sorted(cap.aliases),
            "energy_gradient": cap.energy_gradient,
            "native_optimization": cap.native_optimization,
            "biased_optimization": cap.supports_biased_optimization,
            "charge": cap.supports_charge,
            "multiplicity": cap.supports_multiplicity,
            "executables": sorted(cap.required_executables),
            "python_modules": sorted(cap.required_python_modules),
            "optional_extra": cap.optional_extra_hint or None,
            "notes": cap.notes,
        }
        if args.backend:
            records[name]['availability'] = {
                'executables': {exe: shutil.which(exe) for exe in sorted(cap.required_executables)},
                'python_modules': {module: _module_available(module) for module in sorted(cap.required_python_modules)},
            }
    if args.json:
        print(json.dumps({"backends": records}, indent=2, sort_keys=True))
    else:
        for name, cap in records.items():
            print(f"{name} ({cap['family']})")
            print(f"  energy/gradient: {'yes' if cap['energy_gradient'] else 'no'}; "
                  f"native optimization: {'yes' if cap['native_optimization'] else 'no'}")
            if cap["executables"]:
                print("  executables: " + ", ".join(cap["executables"]))
            if cap["python_modules"]:
                print("  Python: " + ", ".join(cap["python_modules"]))
            if cap["notes"]:
                print("  " + cap["notes"])
            if 'availability' in cap:
                for exe, path in cap['availability']['executables'].items():
                    print(f"  {exe}: {path or 'not found on PATH'}")
                for module, found in cap['availability']['python_modules'].items():
                    print(f"  {module}: {'available' if found else 'not installed'}")
                if cap['optional_extra']:
                    print(f"  Python support: python -m pip install 'pyar-chem[{cap['optional_extra']}]'")


def _module_available(module):
    try:
        return importlib.util.find_spec(module) is not None
    except (ImportError, ValueError, AttributeError):
        return False


def doctor(argv=None, *, prog="pyar doctor"):
    parser = argparse.ArgumentParser(prog=prog, description="Check the PyAR installation")
    parser.add_argument("--all-backends", action="store_true", help="also check optional backend requirements")
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args(argv)
    checks = {}
    for module in ("numpy", "networkx"):
        checks[f"python:{module}"] = {"ok": _module_available(module)}
    if args.all_backends:
        for name, cap in sorted(BACKEND_CAPABILITIES.items()):
            for executable in sorted(cap.required_executables):
                path = shutil.which(executable)
                checks[f"backend:{name}:executable:{executable}"] = {"ok": bool(path), "path": path}
            for module in sorted(cap.required_python_modules):
                found = _module_available(module)
                checks[f"backend:{name}:python:{module}"] = {"ok": found}
    essential = all(check["ok"] for key, check in checks.items() if key.startswith("python:"))
    ok = essential and (not args.all_backends or all(check["ok"] for check in checks.values()))
    result = {"version": __version__, "ok": ok, "checks": checks}
    if args.json:
        print(json.dumps(result, indent=2, sort_keys=True))
    else:
        print(f"PyAR doctor: {'healthy' if ok else 'problems found'} (version {__version__})")
        for label, check in checks.items():
            suffix = f" ({check['path']})" if check.get("path") else ""
            print(f"  {'OK' if check['ok'] else 'MISSING'} {label}{suffix}")
    if not ok:
        raise SystemExit(1)


def inspect_run(argv=None, *, prog="pyar inspect"):
    parser = argparse.ArgumentParser(prog=prog, description="Inspect a PyAR run record or workflow state")
    parser.add_argument("run", help="run directory or pyar-run.toml")
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args(argv)
    path = Path(args.run).expanduser()
    if path.is_dir():
        candidates = _run_records(path) + [path / "pyar-run.json", path / "state.json"]
        record = next((candidate for candidate in candidates if candidate.is_file()), None)
        if record is None:
            parser.error(f"No PyAR run record or state.json found in {path}")
        path = record
    try:
        if path.suffix == ".toml":
            data = load_toml_document(path, RUN_RECORD_KIND)
        else:
            data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        parser.error(str(exc))
    if args.json:
        print(json.dumps(data, indent=2, sort_keys=True))
    else:
        for key in ("workflow", "status", "run_directory", "created_utc", "pyar_version"):
            if key in data:
                print(f"{key.replace('_', ' ').title()}: {data[key]}")
        metadata = data.get("metadata", {})
        for key in ("status", "run_directory", "record_directory", "created_utc", "pyar_version"):
            if key in metadata and key not in data:
                print(f"{key.replace('_', ' ').title()}: {metadata[key]}")
        effective = data.get("effective", {})
        if effective:
            print("Resolved request:")
            for key, value in effective.items():
                print(f"  {key}: {value}")
        state = data.get("state", {})
        if state:
            print("Workflow state:")
            for key, value in state.items():
                print(f"  {key}: {value}")


def run_config(argv=None, *, prog="pyar run", verbosity=0, quiet=False):
    parser = argparse.ArgumentParser(prog=prog, description="Run a saved calculation; task options override saved values", allow_abbrev=False)
    parser.add_argument("file")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--check", action="store_true")
    mode.add_argument("--dry-run", action="store_true")
    parser.add_argument('--json', action='store_true')
    args, overrides = parser.parse_known_args(argv)
    from pyar.modern_cli import _COMMANDS, _META_COMMANDS, _command_parser, invoke
    from importlib import import_module
    try:
        config = load_toml_document(args.file, RUN_CONFIG_KIND)
        workflow = config['workflow']
        if workflow not in _COMMANDS or workflow in _META_COMMANDS:
            raise ValueError(f'Not an executable task: {workflow}')
        task_parser = _command_parser(workflow, import_module(_COMMANDS[workflow][0]))
        if task_parser is None:
            raise ValueError(f'Task {workflow!r} cannot replay configurations')
        effective = config['effective']
        if not effective:
            # Read original schema-1 argument-only files without inheriting --check.
            effective = vars(task_parser.parse_args(without_execution_controls(config['arguments'])))
        override_dests = explicit_destinations(task_parser, overrides)
        positional = {a.dest for a in task_parser._actions if not a.option_strings}
        if override_dests & positional:
            raise ValueError('Override task options by name; edit the configuration to change positional inputs')
        values = discard_overridden(effective, task_parser, override_dests)
        vector = arguments_from_values(values, task_parser)
        boundary = vector.index('--') if '--' in vector else len(vector)
        vector[boundary:boundary] = overrides
        controls = ['--check'] if args.check else ['--dry-run'] if args.dry_run else []
        if args.json:
            controls.append('--json')
        controls.extend(['--verbose'] * verbosity)
        if quiet:
            controls.append('--quiet')
        vector[0:0] = controls
    except (ValueError, OSError) as exc:
        parser.error(str(exc))
    return invoke(workflow, vector, config_file=config['path'],
                  source_overrides={dest: 'cli' for dest in override_dests})


def config(argv=None, *, prog="pyar config"):
    parser = argparse.ArgumentParser(prog=prog, description="Create or validate versioned PyAR run configurations")
    sub = parser.add_subparsers(dest="action", required=True)
    create = sub.add_parser("create", help="Create a reproducible config from a run record")
    create.add_argument("--from", dest="source", default=".")
    create.add_argument("--output", "-o", default="pyar-run.toml")
    create.add_argument("--force", action="store_true", help="Replace an existing configuration")
    check = sub.add_parser("validate", help="Validate a run configuration")
    check.add_argument("file")
    args = parser.parse_args(argv)
    try:
        if args.action == "validate":
            data = load_toml_document(args.file, RUN_CONFIG_KIND)
            from pyar.modern_cli import _CHECK_COMMANDS
            if data['workflow'] not in _CHECK_COMMANDS:
                raise ValueError("Read-only configuration validation is available for calculation workflows; use task help for inspection utilities")
            run_config([args.file, '--check'])
            print(f"Valid PyAR run config: {data['workflow']} ({data['path']})")
            return
        source = _load_run_record(args.source)
        if Path(args.output).exists() and not args.force:
            raise ValueError(f"Configuration already exists: {args.output}; use --force to replace it")
        destination = write_toml_document(
            args.output, kind=RUN_CONFIG_KIND, workflow=source["workflow"],
            arguments=without_execution_controls(source.get("arguments", [])),
            effective={k: v for k, v in source.get("effective", {}).items() if k != "check"},
            metadata={"created_from": str(Path(args.source).resolve()), "created_utc": _now(),
                      "profile_values": source.get("metadata", {}).get("profile_values", {})},
        )
    except (ValueError, OSError, KeyError) as exc:
        parser.error(str(exc))
    print(f"Wrote run config: {destination}")


def profile(argv=None, *, prog="pyar profile"):
    parser = argparse.ArgumentParser(prog=prog, description="Manage reusable partial PyAR methodology profiles")
    sub = parser.add_subparsers(dest="action", required=True)
    sub.add_parser("list", help="List profiles")
    show = sub.add_parser("show", help="Show a profile")
    show.add_argument("name")
    create = sub.add_parser("create", help="Create a profile from a run record")
    create.add_argument("name")
    create.add_argument("--from", dest="source", default=".")
    create.add_argument("--force", action="store_true", help="Replace an existing profile")
    delete = sub.add_parser("delete", help="Delete a named profile")
    delete.add_argument("name")
    args = parser.parse_args(argv)
    directory = profile_directory()
    if args.action == "list":
        for path in sorted(directory.glob("*.toml")):
            try:
                document = load_toml_document(path, PROFILE_KIND)
                print(f"{path.stem}\t{document['workflow']}")
            except ValueError:
                print(f"{path.stem}\tinvalid profile")
        return
    try:
        if args.action == "delete":
            profile_path(args.name).unlink()
            print(f"Deleted profile {args.name}")
            return
        if args.action == "show":
            path = profile_path(args.name)
            value = load_toml_document(path, PROFILE_KIND)
            _emit({"name": args.name, "workflow": value["workflow"],
                   "arguments": value["arguments"], "effective": value.get("effective", {})}, False)
            return
        source = _load_run_record(args.source)
        workflow = source["workflow"]
        from pyar.modern_cli import _COMMANDS
        if workflow not in _COMMANDS:
            raise ValueError(f"Cannot create a profile for unsupported workflow {workflow!r}")
        module = __import__(_COMMANDS[workflow][0], fromlist=["build_parser"])
        parser_builder = getattr(module, "build_parser", None)
        if parser_builder is None:
            raise ValueError(f"Workflow {workflow!r} does not expose a parser for profile extraction")
        from pyar.modern_cli import _command_parser
        command_parser = _command_parser(workflow, module)
        profile_args = arguments_from_profile(source.get("effective", {}), command_parser, workflow=workflow)
        profile_values = profile_overrides_from_record(source)
        if profile_path(args.name).exists() and not args.force:
            raise ValueError(f"Profile {args.name!r} already exists; use --force to replace it")
        destination = write_toml_document(
            profile_path(args.name), kind=PROFILE_KIND, workflow=workflow, arguments=profile_args,
            effective=profile_values,
            metadata={"created_from": str(Path(args.source).resolve()), "created_utc": _now()},
        )
    except (ValueError, OSError, ImportError, TypeError) as exc:
        parser.error(str(exc))
    print(f"Created profile {args.name}: {destination}")


def _now():
    return datetime.now(timezone.utc).isoformat()


def _load_run_record(source):
    path = Path(source).expanduser()
    if path.is_dir():
        records = _run_records(path)
        if not records:
            raise ValueError(f"No run record found in {path}")
        path = records[0]
    value = load_toml_document(path, RUN_RECORD_KIND)
    if not value.get("workflow"):
        raise ValueError("Run record has no workflow name")
    return value


def _run_records(directory):
    """Prefer the latest invocation without overwriting previous run provenance."""
    return sorted((path for path in directory.glob('pyar-run*.toml') if path.is_file()),
                  key=lambda path: (path.stat().st_mtime_ns, path.name), reverse=True)
