"""Task-oriented CLI pilot, alongside the existing PyAR entry points.

The command tree owns task discovery; each command owns its argument parser
and execution. Command modules are imported only when their task is selected.
"""

import argparse
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


# Add future tasks here without importing their scientific dependencies at startup.
_COMMANDS = {
    "run": ("pyar.scripts.cli_meta", "Run a versioned PyAR task configuration"),
    "config": ("pyar.scripts.cli_meta", "Create or validate run configurations"),
    "profile": ("pyar.scripts.cli_meta", "Manage reusable methodology profiles"),
    "info": ("pyar.scripts.cli_meta", "Show PyAR installation information"),
    "backends": ("pyar.scripts.cli_meta", "List backend capabilities"),
    "doctor": ("pyar.scripts.cli_meta", "Check PyAR installation health"),
    "inspect": ("pyar.scripts.cli_meta", "Inspect a PyAR run or workflow state"),
    "conformer": ("pyar.scripts.modern_conformer", "Generate and optionally refine molecular conformers"),
    "aggregate": ("pyar.scripts.modern_aggregate", "Search structures for a requested final composition"),
    "grow": ("pyar.scripts.modern_grow", "Sequentially add one species to a fixed seed"),
    "microsolvate": ("pyar.scripts.modern_microsolvate", "Build explicit solvent shells around a fixed solute"),
    "solvate": ("pyar.scripts.modern_solvate", "Compatibility alias for solute-centred microsolvation"),
    "neb": ("pyar.scripts.modern_path_stage", "Search a reaction path between endpoint structures"),
    "ts": ("pyar.scripts.modern_path_stage", "Optimize and validate a transition-state geometry"),
    "irc": ("pyar.scripts.modern_path_stage", "Continue a validated TS along both IRC directions"),
    "react": ("pyar.scripts.modern_react", "Search for reaction products with adaptive bias"),
    "orient": ("pyar.scripts.orient", "Generate trial encounter geometries for two fragments"),
    "deduplicate": ("pyar.scripts.deduplicate", "Remove structurally duplicate XYZ geometries"),
    "select": ("pyar.scripts.select", "Select low-energy structures"),
    "identify": ("pyar.scripts.identify", "Inspect composition and perceived chemical identity"),
    "split": ("pyar.scripts.split", "Split disconnected coordinate components"),
    "trace": ("pyar.scripts.reaction_trace", "Analyze a PyAR reaction trace"),
    "energies": ("pyar.scripts.energy_table", "Print a relative-energy table from XYZ files"),
    "compare": ("pyar.scripts.compare", "Compare two structures, energies, geometry, and connectivity"),
    "clustering": ("pyar.scripts.clustering", "Cluster or filter XYZ structures"),
    "optimize": ("pyar.scripts.optimize", "Optimize each XYZ structure independently"),
    "scan-bond": ("pyar.scripts.scan_bond", "Relaxed bond scan with optional reaction-path continuation"),
}


_CHECK_COMMANDS = {
    "react", "scan-bond", "optimize", "conformer", "aggregate", "grow",
    "microsolvate", "solvate", "neb", "ts", "irc",
}
_META_COMMANDS = {"run", "config", "profile", "info", "backends", "doctor", "inspect"}
_COMMAND_GROUPS = {
    "Structure exploration": ("conformer", "aggregate", "grow", "microsolvate", "solvate", "optimize"),
    "Reaction exploration": ("react", "scan-bond", "neb", "ts", "irc"),
    "Analysis and geometry": ("clustering", "energies", "compare", "identify", "orient", "deduplicate", "select", "split", "trace"),
    "Configuration and diagnostics": ("run", "config", "profile", "inspect", "info", "backends", "doctor"),
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
        parser = builder(modern=True, prog=f'pyar {command}')
    elif command in {'neb', 'ts', 'irc'}:
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
    status = result.get('status', 'completed') if isinstance(result, dict) else 'completed'
    metadata = {
        'pyar_version': __version__, 'created_utc': datetime.now(timezone.utc).isoformat(),
        'status': 'failed' if failure is not None else status, 'record_id': str(uuid.uuid4()),
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


def invoke(command, arguments, *, config_file=None, profile_values=None, source_overrides=None):
    """Dispatch a task, observing its own validated request and actual result."""
    from pyar.cli_runtime import observe_run
    from pyar.run_config import (PROFILE_KIND, RUN_CONFIG_KIND, ResolvedRunSpec, RunSpec,
                                 arguments_from_profile, build_execution_plan, load_toml_document,
                                 profile_path, resolve_run_spec, without_execution_controls,
                                 write_toml_document)
    error_parser = argparse.ArgumentParser(prog=f'pyar {command}', add_help=False)
    if command not in _COMMANDS:
        error_parser.error(f'Unknown task: {command}')
    args = list(arguments)
    try:
        args, verbosity = _verbosity(args)
        args, quiet_long = _remove_option(args, '--quiet')
        args, quiet_short = _remove_option(args, '-q')
        args, _ = _remove_option(args, '--no-color')
        quiet = bool(quiet_long or quiet_short)
        if quiet and verbosity:
            raise ValueError('--quiet and --verbose cannot be combined')
    except ValueError as exc:
        error_parser.error(str(exc))
    if command in _META_COMMANDS:
        module = import_module(_COMMANDS[command][0])
        function = {'run': module.run_config, 'config': module.config, 'profile': module.profile,
                    'info': module.info, 'backends': module.backends, 'doctor': module.doctor,
                    'inspect': module.inspect_run}[command]
        if command == 'run':
            return function(args, prog=f'pyar {command}', verbosity=verbosity, quiet=quiet)
        with _logging(verbosity, quiet):
            return function(args, prog=f'pyar {command}')
    try:
        args, profile_name = _remove_option(args, '--profile', takes_value=True)
        args, config_path = _remove_option(args, '--write-config', takes_value=True)
        args, dry_run = _remove_option(args, '--dry-run')
        as_json = False
        if command in _CHECK_COMMANDS:
            args, as_json = _remove_option(args, '--json')
        if dry_run and command not in _CHECK_COMMANDS:
            raise ValueError('--dry-run is supported only for calculation workflows with --check')
    except ValueError as exc:
        error_parser.error(str(exc))
    module = import_module(_COMMANDS[command][0])
    parser = _command_parser(command, module)
    option_part = args[:args.index('--')] if '--' in args else args
    if any(token in {'-h', '--help'} for token in option_part) and parser is not None:
        _add_shared_help(parser, dry_run=command in _CHECK_COMMANDS).print_help()
        parser.exit()
    profile_values = dict(profile_values or {})
    try:
        if profile_name:
            profile = load_toml_document(profile_path(profile_name), PROFILE_KIND)
            if profile['workflow'] != command:
                raise ValueError(f"Profile {profile_name!r} is for {profile['workflow']!r}")
            profile_values.update(profile['effective'])
        if config_path and Path(config_path).expanduser().exists():
            raise ValueError(f'Configuration already exists: {config_path}. Choose a new path.')
    except (ValueError, OSError) as exc:
        error_parser.error(str(exc))
    if parser is not None:
        request, resolved = resolve_run_spec(command, parser, args, profile_values=profile_values,
                                             profile=profile_name, argument_source='config' if config_file else 'cli')
    else:
        request = RunSpec(command, tuple(args), profile_name)
        resolved = ResolvedRunSpec(command, {}, {}, tuple(args))
    if source_overrides:
        resolved.value_sources.update(source_overrides)
    execution = list(args)
    if parser is not None:
        applied = {key: value for key, value in profile_values.items() if resolved.value_sources.get(key) == 'profile'}
        extra = arguments_from_profile(applied, parser, workflow=command, skip_defaults=False)
        boundary = execution.index('--') if '--' in execution else len(execution)
        execution[boundary:boundary] = extra
    checking = bool(resolved.values.get('check') or dry_run)
    if dry_run and not resolved.values.get('check'):
        boundary = execution.index('--') if '--' in execution else len(execution)
        execution.insert(boundary, '--check')
    # Configs are calculation specifications, never permanently validation-only.
    request = RunSpec(command, tuple(without_execution_controls(args)), profile_name)
    saved_config = False

    def save_config():
        nonlocal saved_config
        if not config_path or saved_config:
            return
        write_toml_document(config_path, kind=RUN_CONFIG_KIND, workflow=command,
                            arguments=request.arguments, effective={**resolved.values, 'check': False} if 'check' in resolved.values else resolved.values,
                            metadata={'profile': profile_name, 'value_sources': resolved.value_sources,
                                      'resolved_request': resolved.scientific})
        saved_config = True
        print(f'Wrote run configuration: {Path(config_path).expanduser().resolve()}', file=sys.stderr)

    def prepared(observation):
        nonlocal resolved
        values = dict(resolved.values)
        sources = dict(resolved.value_sources)
        for key, value in observation.options.items():
            if key not in values:
                continue
            if values[key] != value and sources.get(key) == 'default':
                sources[key] = 'workflow-default'
            values[key] = value
        resolved = ResolvedRunSpec(command, values, sources, request.arguments, observation.prepared)
        save_config()

    if command in {'neb', 'ts', 'irc'}:
        entry = getattr(module, 'modern_' + command)
    else:
        entry = getattr(module, 'modern_main', None) or module.main
    result, failure = None, None
    with observe_run(prepared) as observation, _logging(verbosity, quiet), ExitStack() as stack:
        if as_json or (quiet and command in _CHECK_COMMANDS):
            stack.enter_context(redirect_stdout(io.StringIO() if checking or quiet else sys.stderr))
        try:
            result = entry(execution, prog=f'pyar {command}')
            if observation.result is None and result is not None:
                observation.result = result.to_dict() if hasattr(result, 'to_dict') else result
            save_config()
        except BaseException as exc:
            failure = exc
        finally:
            if observation.started and not checking:
                try:
                    _write_execution_record(command, request, resolved, observation, failure=failure, config_file=config_file)
                except (OSError, ValueError) as exc:
                    print(f'pyar {command}: could not write run record: {exc}', file=sys.stderr)
    plan = build_execution_plan(command, resolved) if dry_run and failure is None else None
    if as_json:
        payload = {'schema_version': 1, 'workflow': command,
                   'status': 'failed' if failure else 'ready' if checking else
                   observation.result.get('status', 'completed') if isinstance(observation.result, dict) else 'completed',
                   'mode': 'dry-run' if dry_run else 'check' if checking else 'run',
                   'resolved': resolved.to_dict(), 'result': observation.result}
        if plan:
            payload['plan'] = plan.to_dict()
        if failure:
            payload['error'] = str(failure)
        print(json.dumps(payload, indent=2, sort_keys=True, allow_nan=False, default=str))
    elif plan and not quiet:
        _print_execution_plan(plan, resolved.scientific)
    if failure is not None:
        raise failure
    return result


if __name__ == '__main__':
    main()
