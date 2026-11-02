"""Task-oriented CLI pilot, alongside the existing PyAR entry points.

The command tree owns task discovery; each command owns its argument parser
and execution. Command modules are imported only when their task is selected.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from typing import TYPE_CHECKING, Callable
import io
import json
from contextlib import ExitStack, contextmanager, redirect_stdout
from importlib import import_module
import sys
import logging
from datetime import datetime, timezone
from pathlib import Path
import uuid

from pyar import __version__

if TYPE_CHECKING:
    from pyar.cli_runtime import RunObservation
    from pyar.run_config import ResolvedRunSpec, RunSpec


# Add future tasks here without importing their scientific dependencies at startup.
_COMMANDS = {
    "explore": ("pyar.scripts.explore", "Generate formula-directed composite geometries"),
    "trial": ("pyar.scripts.trial_generation", "Generate trial directions and placement geometries"),
    "basin-memory": ("pyar.scripts.basin_memory", "Inspect and maintain basin registries"),
    "similarity": ("pyar.scripts.similarity", "Compare pools with Grigoryan–Springborg similarity"),
    "descriptor": ("pyar.scripts.descriptor", "Calculate cluster descriptors"),
    "mlopt": ("pyar.backends.mlopt", "Optimize with the MLatom AIQM1 driver"),
    "aimnet2-ase-opt": ("pyar.backends.aimnet2_ase_opt", "Optimize with a supplied AIMNet2 model"),
    "benchmark": ("pyar.scripts.benchmark", "Run developer benchmark campaigns"),
    "run": ("pyar.scripts.cli_meta", "Run a versioned PyAR task configuration"),
    "config": ("pyar.scripts.cli_meta", "Create or validate run configurations"),
    "profile": ("pyar.scripts.cli_meta", "Manage reusable methodology profiles"),
    "info": ("pyar.scripts.cli_meta", "Show PyAR installation information"),
    "backends": ("pyar.scripts.cli_meta", "List backend capabilities"),
    "doctor": ("pyar.scripts.cli_meta", "Check PyAR installation health"),
    "install": ("pyar.scripts.cli_meta", "Show installation status for Python packages and external programs"),
    "inspect": ("pyar.scripts.cli_meta", "Inspect a PyAR run or workflow state"),
    "conformer": ("pyar.scripts.conformer", "Generate and optionally refine molecular conformers"),
    "aggregate": ("pyar.scripts.aggregate", "Search structures for a requested final composition"),
    "grow": ("pyar.scripts.grow", "Sequentially add one species to a fixed seed"),
    "microsolvate": ("pyar.scripts.microsolvate", "Build explicit solvent shells around a fixed solute"),
    "neb": ("pyar.scripts.path_stage", "Search a reaction path between endpoint structures"),
    "ts": ("pyar.scripts.path_stage", "Optimize and validate a transition-state geometry"),
    "irc": ("pyar.scripts.path_stage", "Continue a validated TS along both IRC directions"),
    "react": ("pyar.scripts.react", "Search for reaction products with adaptive bias"),
    "orient": ("pyar.scripts.orient", "Generate trial encounter geometries for two fragments"),
    "deduplicate": ("pyar.scripts.deduplicate", "Remove structurally duplicate XYZ geometries"),
    "select": ("pyar.scripts.select", "Select low-energy structures"),
    "identify": ("pyar.scripts.identify", "Inspect composition and perceived chemical identity"),
    "split": ("pyar.scripts.split", "Split disconnected coordinate components"),
    "trace": ("pyar.scripts.reaction_trace", "Analyze a PyAR reaction trace"),
    "energies": ("pyar.scripts.energy_table", "Print a relative-energy table from XYZ files"),
    "compare": ("pyar.scripts.compare", "Compare two structures, energies, geometry, and connectivity"),
    "clustering": ("pyar.scripts.clustering", "Cluster or filter XYZ structures"),
    "energy": ("pyar.scripts.singlepoint", "Calculate single-point energies for fixed geometries"),
    "engrad": ("pyar.scripts.singlepoint", "Calculate energies and Cartesian gradients"),
    "gradient": ("pyar.scripts.singlepoint", "Alias for engrad"),
    "frequency": ("pyar.scripts.singlepoint", "Calculate harmonic frequencies and qualify minima"),
    "optimize": ("pyar.scripts.optimize", "Optimize each XYZ structure independently"),
    "scan-bond": ("pyar.scripts.scan_bond", "Relaxed bond scan with optional reaction-path continuation"),
}


_CHECK_COMMANDS = {
    "react", "scan-bond", "optimize", "conformer", "aggregate", "grow",
    "energy", "engrad", "gradient", "frequency",
    "microsolvate", "neb", "ts", "irc",
}
_META_COMMANDS = {"run", "config", "profile", "info", "backends", "doctor", "install", "inspect", "benchmark"}
_COMMAND_GROUPS = {
    "Structure exploration": ("conformer", "aggregate", "grow", "microsolvate", "optimize"),
    "Fixed-geometry calculations": ("energy", "engrad", "gradient", "frequency"),
    "Reaction exploration": ("react", "scan-bond", "neb", "ts", "irc"),
    "Analysis and geometry": ("clustering", "energies", "compare", "identify", "orient", "deduplicate", "select", "split", "trace"),
    "Specialist tools": ("explore", "trial", "basin-memory", "similarity", "descriptor", "mlopt", "aimnet2-ase-opt", "benchmark"),
    "Configuration and diagnostics": ("run", "config", "profile", "inspect", "info", "backends", "doctor", "install"),
}


class _TaskParser(argparse.ArgumentParser):
    def format_help(self):
        lines = ["usage: pyar [OPTIONS] COMMAND [ARGS]", "", self.description, ""]
        for group, commands in _COMMAND_GROUPS.items():
            lines.append(group + ":")
            for name in commands:
                lines.append(f"  {name:<14} {_COMMANDS[name][1]}")
            lines.append("")
        lines += ["Options: --version, -h/--help, -q/--quiet, -v/-vv/-vvv, --no-color",
                  "Use 'pyar COMMAND --help' for command options.", ""]
        return "\n".join(lines)


def build_parser():
    parser = _TaskParser(prog="pyar", description="PyAR molecular structure and reaction exploration", allow_abbrev=False)
    if hasattr(parser, 'color'):
        parser.color = False
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("-v", "--verbose", action="count", default=0)
    parser.add_argument("-q", "--quiet", action="store_true")
    parser.add_argument("--no-color", action="store_true")
    commands = parser.add_subparsers(dest="command", required=True, metavar="COMMAND")
    for name in _COMMANDS:
        commands.add_parser(name, add_help=False, allow_abbrev=False)
    return parser


def main(argv=None):
    """Console entry point: never pass a domain object to sys.exit()."""
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] == 'help':
        argv = [*argv[1:], '--help'] if len(argv) > 1 else ['--help']
    args, remainder = build_parser().parse_known_args(argv)
    controls = ['-v'] * args.verbose
    if args.quiet:
        controls.append('--quiet')
    if args.no_color:
        controls.append('--no-color')
    try:
        invoke(args.command, [*controls, *remainder])
    except KeyboardInterrupt:
        print('Interrupted.', file=sys.stderr)
        raise SystemExit(130) from None
    except BrokenPipeError:
        # Downstream consumers (e.g. head) may close stdout early.
        return


def _remove_option(arguments, name, *, takes_value=False):
    result, value = [], None
    index = 0
    while index < len(arguments):
        token = arguments[index]
        if token == '--':
            result.extend(arguments[index:])
            break
        if token == name:
            if takes_value:
                if index + 1 >= len(arguments) or arguments[index + 1].startswith('--'):
                    raise ValueError(f"{name} requires a value")
                value = arguments[index + 1]
                index += 2
            else:
                value = True
                index += 1
            continue
        if takes_value and token.startswith(name + "="):
            value = token.split("=", 1)[1]
            index += 1
            continue
        result.append(token)
        index += 1
    return result, value


def _command_parser(command, module):
    builder = getattr(module, 'build_parser', None) or getattr(module, '_build_parser', None)
    if builder is None:
        return None
    if command in {'scan-bond', 'energies'}:
        parser = builder(prog=f'pyar {command}')
    elif command in {'neb', 'ts', 'irc', 'energy', 'engrad', 'gradient', 'frequency'}:
        parser = builder(command, prog=f'pyar {command}')
    else:
        parser = builder(prog=f'pyar {command}')
    parser.allow_abbrev = False
    if hasattr(parser, 'color'):
        parser.color = False
    return parser


def _add_shared_help(parser, *, dry_run=True):
    group = parser.add_argument_group('Common execution controls')
    group.add_argument('--profile', metavar='NAME', help='Apply methodology settings; explicit CLI values win')
    group.add_argument('--write-config', metavar='FILE', help='Save a validated configuration (refuses existing files)')
    if dry_run:
        group.add_argument('--dry-run', action='store_true', help='Preflight and display the resolved execution plan')
        group.add_argument('--json', action='store_true', help='Emit a structured request/result; progress goes to stderr')
    group.add_argument('-v', '--verbose', action='count', default=0, help='Increase detail (-vv/-vvv for diagnostics)')
    group.add_argument('-q', '--quiet', action='store_true', help='Suppress workflow progress; preserve primary analysis data')
    group.add_argument('--no-color', action='store_true', help='Plain output (currently the default)')
    return parser


def _verbosity(arguments):
    remaining, level = [], 0
    for index, token in enumerate(arguments):
        if token == '--':
            remaining.extend(arguments[index:])
            break
        if token == '--verbose':
            level += 1
        elif token.startswith('-v') and set(token[1:]) == {'v'}:
            level += len(token) - 1
        else:
            remaining.append(token)
    return remaining, level


@contextmanager
def _logging(level, quiet):
    logger = logging.getLogger('pyar')
    previous = logger.level
    handler = None
    if level or quiet:
        logger.setLevel(logging.ERROR if quiet else logging.DEBUG if level > 1 else logging.INFO)
        handler = logging.StreamHandler(sys.stderr)
        handler.setFormatter(logging.Formatter('%(levelname)s: %(message)s'))
        logger.addHandler(handler)
    try:
        yield
    finally:
        if handler is not None:
            logger.removeHandler(handler)
            handler.close()
        logger.setLevel(previous)


def _print_execution_plan(plan, scientific):
    print(f'Execution plan: {plan.workflow}')
    if plan.inputs:
        print('Inputs: ' + ', '.join(plan.inputs))
    print(f"Backend: {plan.backend or 'none (geometry-only or RDKit generation)'}")
    print('Resolved command settings:')
    for key, value in sorted(plan.settings.items()):
        print(f"  {key.replace('_', '-')}: {value} [{plan.value_sources.get(key, 'resolved')}]")
    if scientific.get('backend'):
        print('Resolved backend settings:')
        for key, value in sorted(scientific['backend'].items()):
            print(f'  {key}: {value}')
    print('Planned stages:')
    for index, stage in enumerate(plan.stages, 1):
        print(f'  {index}. {stage}')
    if plan.outputs:
        print('Output locations: ' + ', '.join(plan.outputs))
    print('No calculation or workflow output was created.')


def _write_execution_record(command, request, resolved, observation, *, failure=None, config_file=None):
    from pyar.run_config import RUN_RECORD_KIND, write_toml_document
    result = observation.result or {}
    directory = result.get('run_directory') if isinstance(result, dict) else None
    if directory is None:
        outputs = resolved.scientific.get('outputs', [])
        directory = outputs[0] if command != 'optimize' and len(outputs) == 1 else None
    if directory is None or not Path(directory).is_dir():
        run_id = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ-') + uuid.uuid4().hex[:12]
        directory = Path.cwd() / '.pyar' / 'runs' / run_id
    directory = Path(directory).resolve()
    path = directory / 'pyar-run.toml'
    if path.exists():
        path = directory / f'pyar-run-{uuid.uuid4().hex}.toml'
    status = (
        result.get('status') if isinstance(result, dict) and result.get('status') else
        'failed' if failure is not None else 'completed'
    )
    metadata = {
        'pyar_version': __version__, 'created_utc': datetime.now(timezone.utc).isoformat(),
        'status': status, 'scientific_status': status,
        'invocation_status': 'failed' if failure is not None else 'completed',
        'record_id': str(uuid.uuid4()),
        'record_directory': str(directory), 'profile': request.profile,
        'value_sources': resolved.value_sources, 'resolved_request': resolved.scientific,
        'config_file': config_file, 'result': result,
    }
    if isinstance(result, dict) and result.get('run_directory'):
        metadata['run_directory'] = result['run_directory']
    if failure is not None:
        metadata['error'] = str(failure)
    return write_toml_document(path, kind=RUN_RECORD_KIND, workflow=command,
                               arguments=request.arguments, effective=resolved.values, metadata=metadata)


def _execution_json_payload(command, resolved, observation, *, checking, dry_run,
                            plan, failure, record_failure, provenance_only_failure):
    """Build the stable JSON response from resolved request and run evidence."""
    scientific_status = (
        observation.result.get('status', 'completed')
        if isinstance(observation.result, dict) else
        'ready' if checking and failure is None else
        'failed' if failure else 'completed'
    )
    invocation_status = (
        'provenance_write_failed' if provenance_only_failure else
        'failed' if failure else 'ready' if checking else 'completed'
    )
    payload = {'schema_version': 1, 'workflow': command,
               'status': scientific_status, 'scientific_status': scientific_status,
               'invocation_status': invocation_status,
               'mode': 'dry-run' if dry_run else 'check' if checking else 'run',
               'resolved': resolved.to_dict(), 'result': observation.result}
    if plan:
        payload['plan'] = plan.to_dict()
    if failure:
        payload['error'] = str(failure)
    if record_failure is not None:
        payload['provenance_error'] = str(record_failure)
    return payload


@dataclass(frozen=True, slots=True)
class _ConsoleControls:
    arguments: list[str]
    verbosity: int
    quiet: bool


@dataclass(frozen=True, slots=True)
class _InvocationOptions:
    command: str
    console: _ConsoleControls
    arguments: list[str]
    profile_name: str | None
    config_path: str | None
    dry_run: bool
    as_json: bool
    config_file: str | Path | None
    error_parser: argparse.ArgumentParser


@dataclass(slots=True)
class _ConfiguredInvocation:
    """Own the evolving prepared request and its one-time config snapshot."""

    options: _InvocationOptions
    request: RunSpec
    resolved: ResolvedRunSpec
    execution: list[str]
    checking: bool
    saved_config: bool = False

    def save_config(self):
        from pyar.run_config import RUN_CONFIG_KIND, write_toml_document

        path = self.options.config_path
        if not path or self.saved_config:
            return
        write_toml_document(
            path,
            kind=RUN_CONFIG_KIND,
            workflow=self.options.command,
            arguments=self.request.arguments,
            effective={**self.resolved.values, "check": False}
            if "check" in self.resolved.values
            else self.resolved.values,
            metadata={
                "profile": self.options.profile_name,
                "value_sources": self.resolved.value_sources,
                "resolved_request": self.resolved.scientific,
            },
        )
        self.saved_config = True
        print(
            f"Wrote run configuration: {Path(path).expanduser().resolve()}",
            file=sys.stderr,
        )

    def prepared(self, observation: RunObservation):
        from pyar.run_config import ResolvedRunSpec

        values = dict(self.resolved.values)
        sources = dict(self.resolved.value_sources)
        for key, value in observation.options.items():
            if key not in values:
                continue
            if values[key] != value and sources.get(key) == "default":
                sources[key] = "workflow-default"
            values[key] = value
        self.resolved = ResolvedRunSpec(
            self.options.command,
            values,
            sources,
            self.request.arguments,
            observation.prepared,
        )
        self.save_config()


@dataclass(slots=True)
class _InvocationOutcome:
    observation: RunObservation
    result: object = None
    failure: BaseException | None = None
    record_failure: OSError | ValueError | None = None
    provenance_only_failure: bool = False


def _console_controls(arguments, error_parser):
    try:
        args, verbosity = _verbosity(arguments)
        args, quiet_long = _remove_option(args, "--quiet")
        args, quiet_short = _remove_option(args, "-q")
        args, _ = _remove_option(args, "--no-color")
        quiet = bool(quiet_long or quiet_short)
        if quiet and verbosity:
            raise ValueError("--quiet and --verbose cannot be combined")
    except ValueError as exc:
        error_parser.error(str(exc))
    return _ConsoleControls(args, verbosity, quiet)


def _invoke_meta_command(command, console):
    module = import_module(_COMMANDS[command][0])
    if command == "benchmark":
        return module.main(console.arguments, prog="pyar benchmark")
    function = {
        "run": module.run_config,
        "config": module.config,
        "profile": module.profile,
        "info": module.info,
        "backends": module.backends,
        "doctor": module.doctor,
        "inspect": module.inspect_run,
        "install": module.install,
    }[command]
    if command == "run":
        return function(
            console.arguments,
            prog=f"pyar {command}",
            verbosity=console.verbosity,
            quiet=console.quiet,
        )
    with _logging(console.verbosity, console.quiet):
        return function(console.arguments, prog=f"pyar {command}")


def _workflow_options(command, console, error_parser, config_file):
    try:
        args, profile_name = _remove_option(
            console.arguments, "--profile", takes_value=True
        )
        args, config_path = _remove_option(args, "--write-config", takes_value=True)
        args, dry_run = _remove_option(args, "--dry-run")
        as_json = False
        if command in _CHECK_COMMANDS:
            args, as_json = _remove_option(args, "--json")
        if dry_run and command not in _CHECK_COMMANDS:
            raise ValueError(
                "--dry-run is supported only for calculation workflows with --check"
            )
    except ValueError as exc:
        error_parser.error(str(exc))
    return _InvocationOptions(
        command,
        console,
        args,
        profile_name,
        config_path,
        bool(dry_run),
        bool(as_json),
        config_file,
        error_parser,
    )


def _load_invocation_profile(options, profile_values):
    from pyar.run_config import PROFILE_KIND, load_toml_document, profile_path

    values = dict(profile_values or {})
    try:
        if options.profile_name:
            profile = load_toml_document(
                profile_path(options.profile_name), PROFILE_KIND
            )
            if profile["workflow"] != options.command:
                raise ValueError(
                    f"Profile {options.profile_name!r} is for {profile['workflow']!r}"
                )
            values.update(profile["effective"])
        if options.config_path and Path(options.config_path).expanduser().exists():
            raise ValueError(
                f"Configuration already exists: {options.config_path}. Choose a new path."
            )
    except (ValueError, OSError) as exc:
        options.error_parser.error(str(exc))
    return values


def _execution_arguments(options, parser, resolved, profile_values):
    from pyar.run_config import arguments_from_profile

    execution = list(options.arguments)
    if parser is not None:
        applied = {
            key: value
            for key, value in profile_values.items()
            if resolved.value_sources.get(key) == "profile"
        }
        extra = arguments_from_profile(
            applied, parser, workflow=options.command, skip_defaults=False
        )
        boundary = execution.index("--") if "--" in execution else len(execution)
        execution[boundary:boundary] = extra
    checking = bool(resolved.values.get("check") or options.dry_run)
    if options.dry_run and not resolved.values.get("check"):
        boundary = execution.index("--") if "--" in execution else len(execution)
        execution.insert(boundary, "--check")
    return execution, checking


def _configure_invocation(options, module, profile_values, source_overrides):
    from pyar.run_config import (
        ResolvedRunSpec,
        RunSpec,
        resolve_run_spec,
        without_execution_controls,
    )

    command, args = options.command, options.arguments
    parser = _command_parser(command, module)
    option_part = args[: args.index("--")] if "--" in args else args
    if any(token in {"-h", "--help"} for token in option_part) and parser is not None:
        _add_shared_help(parser, dry_run=command in _CHECK_COMMANDS).print_help()
        parser.exit()
    values = _load_invocation_profile(options, profile_values)
    if parser is not None:
        _, resolved = resolve_run_spec(
            command,
            parser,
            args,
            profile_values=values,
            profile=options.profile_name,
            argument_source="config" if options.config_file else "cli",
        )
    else:
        resolved = ResolvedRunSpec(command, {}, {}, tuple(args))
    if source_overrides:
        resolved.value_sources.update(source_overrides)
    execution, checking = _execution_arguments(options, parser, resolved, values)
    # Configs remain calculation specifications, never permanently validation-only.
    request = RunSpec(
        command, tuple(without_execution_controls(args)), options.profile_name
    )
    return _ConfiguredInvocation(options, request, resolved, execution, checking)


def _record_invocation(run, outcome):
    """Preserve the workflow failure when provenance publication also fails."""
    try:
        _write_execution_record(
            run.options.command,
            run.request,
            run.resolved,
            outcome.observation,
            failure=outcome.failure,
            config_file=run.options.config_file,
        )
    except (OSError, ValueError) as exc:
        outcome.record_failure = exc
        message = (
            f"pyar {run.options.command}: calculation outputs were retained, but run provenance "
            f"could not be written: {exc}. Check directory permissions and free space."
        )
        print(message, file=sys.stderr)
        if outcome.failure is None:
            outcome.provenance_only_failure = True
            outcome.failure = OSError(message)


def _execute_invocation(run, entry: Callable[..., object]):
    from pyar.cli_runtime import observe_run

    options, console = run.options, run.options.console
    with (
        observe_run(run.prepared) as observation,
        _logging(console.verbosity, console.quiet),
        ExitStack() as stack,
    ):
        outcome = _InvocationOutcome(observation)
        if options.as_json or (console.quiet and options.command in _CHECK_COMMANDS):
            stack.enter_context(
                redirect_stdout(
                    io.StringIO() if run.checking or console.quiet else sys.stderr
                )
            )
        try:
            outcome.result = entry(run.execution, prog=f"pyar {options.command}")
            if observation.result is None and outcome.result is not None:
                observation.result = (
                    outcome.result.to_dict()
                    if hasattr(outcome.result, "to_dict")
                    else outcome.result
                )
            run.save_config()
        except BaseException as exc:
            outcome.failure = exc
        finally:
            if observation.started and not run.checking:
                _record_invocation(run, outcome)
    return outcome


def _finish_invocation(run, outcome):
    from pyar.run_config import build_execution_plan

    options = run.options
    plan = (
        build_execution_plan(options.command, run.resolved)
        if options.dry_run and outcome.failure is None
        else None
    )
    if options.as_json:
        payload = _execution_json_payload(
            options.command,
            run.resolved,
            outcome.observation,
            checking=run.checking,
            dry_run=options.dry_run,
            plan=plan,
            failure=outcome.failure,
            record_failure=outcome.record_failure,
            provenance_only_failure=outcome.provenance_only_failure,
        )
        print(
            json.dumps(payload, indent=2, sort_keys=True, allow_nan=False, default=str)
        )
    elif plan and not options.console.quiet:
        _print_execution_plan(plan, run.resolved.scientific)
    if outcome.failure is not None:
        raise outcome.failure
    return outcome.result


def invoke(
    command, arguments, *, config_file=None, profile_values=None, source_overrides=None
):
    """Resolve controls, observe task execution, then publish its existing response."""
    # Preserve dependency import failures before command/control validation.
    from pyar.cli_runtime import observe_run  # noqa: F401
    from pyar import run_config  # noqa: F401

    error_parser = argparse.ArgumentParser(prog=f"pyar {command}", add_help=False)
    if command not in _COMMANDS:
        error_parser.error(f"Unknown task: {command}")
    console = _console_controls(list(arguments), error_parser)
    if command in _META_COMMANDS:
        return _invoke_meta_command(command, console)
    options = _workflow_options(command, console, error_parser, config_file)
    module = import_module(_COMMANDS[command][0])
    run = _configure_invocation(options, module, profile_values, source_overrides)
    entry = getattr(module, command) if command in {"neb", "ts", "irc"} else module.main
    outcome = _execute_invocation(run, entry)
    return _finish_invocation(run, outcome)


if __name__ == "__main__":
    main()
