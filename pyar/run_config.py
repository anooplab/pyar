"""Versioned resolved run configurations and reusable CLI profiles.

Configurations retain validated settings as editable TOML plus their argument
provenance. Replay uses each task's normal request resolver and preflight.
"""

from __future__ import annotations

import argparse
import copy
import json
import os
from pathlib import Path
import re
import tempfile
try:
    import tomllib
except ModuleNotFoundError:  # Python 3.10, supported by pyar-chem.
    import tomli as tomllib
from dataclasses import asdict, dataclass, field


SCHEMA_VERSION = 1
RUN_CONFIG_KIND = "pyar-run-config"
RUN_RECORD_KIND = "pyar-run-record"
PROFILE_KIND = "pyar-profile"
PROFILE_EXCLUDED_DESTS = {
    "input_files", "inputs", "input", "input_file", "solute", "solvent",
    "seed", "monomer", "start", "end", "ts_guess", "geometry", "run_directory",
    "run", "first", "second", "path", "formula",
    "output", "output_dir", "labels_output", "report_output", "structure_report",
    "plot_directory", "check", "dry_run", "write_config", "profile",
    "verbosity", "verbose", "quiet", "no_color",
}

EXECUTION_CONTROL_DESTS = {"check", "dry_run", "write_config", "profile", "verbosity", "verbose", "quiet", "no_color"}


def without_execution_controls(arguments):
    """A saved calculation must not permanently inherit a validation-only mode."""
    result = []
    for index, token in enumerate(arguments):
        if token == "--":
            result.extend(arguments[index:])
            break
        if token not in {"--check", "--dry-run"}:
            result.append(token)
    return result


@dataclass(frozen=True)
class RunSpec:
    """User's task choice and argument sources before task-specific resolution."""

    workflow: str
    arguments: tuple[str, ...]
    profile: str | None = None


@dataclass(frozen=True)
class ResolvedRunSpec:
    """Parser-resolved task request with provenance for each option.

    Scientific request objects remain owned by each workflow. This shared
    envelope records the exact CLI-level values that were handed to that
    resolver and where each value came from.
    """

    workflow: str
    values: dict
    value_sources: dict[str, str]
    arguments: tuple[str, ...]
    scientific: dict = field(default_factory=dict)

    def to_dict(self):
        return asdict(self)


@dataclass(frozen=True)
class ExecutionPlan:
    """Side-effect-free outline of an already validated workflow request.

    ``--dry-run`` uses this envelope only after the command's normal resolver
    and preflight have succeeded.  It is intentionally an execution outline,
    not a second scientific request model.
    """

    workflow: str
    inputs: tuple[str, ...]
    backend: str | None
    stages: tuple[str, ...]
    outputs: tuple[str, ...]
    settings: dict
    value_sources: dict[str, str]

    def to_dict(self):
        return asdict(self)


_PLAN_STAGES = {
    "optimize": ("Validate structures", "Optimize each structure", "Write optimized geometries"),
    "conformer": ("Generate conformers", "Optionally refine with the backend", "Select and write conformers"),
    "aggregate": ("Generate composition pathways", "Sample and select candidates", "Write selected final structures"),
    "grow": ("Initialize the seed pool", "Repeat add, filter, and bounded selection", "Retain every growth stage"),
    "microsolvate": ("Construct original-solute surface", "Place solvent candidates against that surface", "Select first-shell structures by coverage and quality"),
    "solvate": ("Construct original-solute surface", "Place solvent candidates against that surface", "Select first-shell structures by coverage and quality"),
    "react": ("Generate reactant orientations", "Run biased reaction search", "Relax and identify candidate products"),
    "scan-bond": ("Generate orientations", "Run relaxed bond scans", "Run requested continuation stages"),
    "neb": ("Relax and validate endpoints", "Interpolate endpoint path", "Optimize and validate the reaction path"),
    "ts": ("Optimize transition-state candidate", "Validate the first-order saddle"),
    "irc": ("Validate saved first-order saddle and artifacts", "Run both IRC directions"),
}

_INPUT_DESTINATIONS = (
    "inputs", "input_files", "input", "input_file", "seed", "monomer",
    "solute", "solvent", "start", "end", "ts_guess", "geometry", "run", "run_directory",
)


def build_execution_plan(workflow, resolved):
    """Build a generic plan envelope from shared parser resolution.

    Workflow-specific ``--check`` remains the authority for scientific
    resolution and capability checks. This function adds the common, stable
    summary without importing workflow engines or causing filesystem writes.
    """
    values = dict(resolved.values)
    inputs = []
    for name in _INPUT_DESTINATIONS:
        if workflow == 'conformer' and name == 'seed':
            continue
        value = values.get(name)
        if value is None:
            continue
        items = value if isinstance(value, (list, tuple)) else [value]
        inputs.extend(str(item) for item in items if item is not None)

    backend = resolved.scientific.get('backend', {}).get('software') or values.get("backend") or values.get("software")
    output = values.get("output")
    if output:
        outputs = (str(Path(output).expanduser()),)
    elif workflow == "optimize" and inputs:
        outputs = tuple(f"job_{Path(item).stem}/" for item in inputs)
    else:
        default_outputs = {
            "react": "reaction/", "scan-bond": "scan_bond/", "aggregate": "aggregates/",
            "grow": "grow/", "microsolvate": "microsolvation/", "solvate": "microsolvation/",
            "conformer": "conformers/",
        }
        outputs = (default_outputs[workflow],) if workflow in default_outputs else ()

    if resolved.scientific.get('outputs'):
        outputs = tuple(resolved.scientific['outputs'])
    excluded = set(_INPUT_DESTINATIONS) | {
        "check", "dry_run", "verbose", "verbosity", "profile", "write_config",
        "output", "json", "help",
    }
    if workflow == 'conformer':
        excluded.discard('seed')
    settings = {
        key: value for key, value in values.items()
        if key not in excluded and value is not None
    }
    return ExecutionPlan(
        workflow=workflow,
        inputs=tuple(inputs),
        backend=backend,
        stages=_PLAN_STAGES.get(workflow, ("Run validated workflow",)),
        outputs=outputs,
        settings=settings,
        value_sources=dict(resolved.value_sources),
    )


def explicit_destinations(parser, arguments):
    """Use argparse itself to recognize aliases, attached values, and ``--``.

    A probe with suppressed defaults distinguishes an explicit option from an
    omitted one, including explicit values equal to defaults. Required options
    may be supplied by another source, so disable that constraint on the probe.
    """
    probe = copy.deepcopy(parser)
    probe._defaults.clear()
    for action in probe._actions:
        action.default = argparse.SUPPRESS
        if action.option_strings:
            action.required = False
        elif action.nargs not in ('*', '?'):
            action.nargs = '*'
    for group in probe._mutually_exclusive_groups:
        group.required = False
    values = vars(probe.parse_args(arguments))
    return {key for key, value in values.items() if value != []}


def discard_overridden(values, parser, destinations):
    """Explicit options also override profile/config alternatives in a mutex group."""
    discarded = set(destinations)
    for group in parser._mutually_exclusive_groups:
        members = {action.dest for action in group._group_actions}
        if members & discarded:
            discarded.update(members)
    return {key: value for key, value in values.items() if key not in discarded}


def resolve_run_spec(workflow, parser, arguments, *, profile_values=None, profile=None,
                     argument_source="cli"):
    """Validate saved methodology through the same argparse actions as CLI values."""
    arguments = tuple(arguments)
    profile_values = dict(profile_values or {})
    known = {action.dest for action in parser._actions} | set(parser._defaults)
    unknown = set(profile_values) - known
    if unknown:
        parser.error("Unknown profile settings: " + ", ".join(sorted(unknown)))
    action_dests = {action.dest for action in parser._actions}
    for key in set(profile_values) & (set(parser._defaults) - action_dests):
        if profile_values[key] != parser._defaults[key]:
            parser.error(f"Profile cannot change fixed setting {key!r}")
    explicit = explicit_destinations(parser, arguments)
    applied = discard_overridden(profile_values, parser, explicit)
    try:
        extra = arguments_from_profile(applied, parser, workflow=workflow, skip_defaults=False)
    except ValueError as exc:
        parser.error(str(exc))
    vector = list(arguments)
    boundary = vector.index('--') if '--' in vector else len(vector)
    vector[boundary:boundary] = extra
    namespace = vars(parser.parse_args(vector))
    sources = {key: argument_source if key in explicit else
               'profile' if key in applied else 'default' for key in namespace}
    request = RunSpec(workflow, arguments, profile)
    return request, ResolvedRunSpec(workflow, namespace, sources, request.arguments)


def write_toml_document(path, *, kind, workflow, arguments, effective=None, metadata=None):
    """Atomically write the stable, deliberately small TOML envelope."""
    path = Path(path).expanduser().resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    document = {
        "schema_version": SCHEMA_VERSION,
        "kind": kind,
        "workflow": workflow,
        "arguments": list(arguments),
        "metadata_json": json.dumps(metadata or {}, sort_keys=True, default=str),
    }
    if kind == RUN_CONFIG_KIND:
        document['unset_settings'] = sorted(key for key, value in (effective or {}).items() if value is None)
    else:
        document['effective_json'] = json.dumps(effective or {}, sort_keys=True, default=str)
    lines = []
    for key, value in document.items():
        encoded = json.dumps(value, ensure_ascii=False)
        lines.append(f"{key} = {encoded}")
    if kind == RUN_CONFIG_KIND:
        lines.extend(['', '[settings]'])
        for key, value in sorted((effective or {}).items()):
            if value is not None:
                lines.append(f'{json.dumps(key)} = {_toml_value(value)}')
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent, text=True)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as stream:
            stream.write("\n".join(lines) + "\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)
    return path


def _toml_value(value):
    """Encode the flat argparse types supported by user-editable run settings."""
    if isinstance(value, (list, tuple)):
        return '[' + ', '.join(_toml_value(item) for item in value) + ']'
    if isinstance(value, (str, bool, int, float)):
        return json.dumps(value, ensure_ascii=False, allow_nan=False)
    raise ValueError(f'Unsupported configuration value: {value!r}')


def load_toml_document(path, expected_kind):
    path = Path(path).expanduser().resolve()
    try:
        with path.open("rb") as stream:
            value = tomllib.load(stream)
    except (OSError, tomllib.TOMLDecodeError) as exc:
        raise ValueError(f"Cannot read {path}: {exc}") from exc
    if value.get("schema_version") != SCHEMA_VERSION or value.get("kind") != expected_kind:
        raise ValueError(f"Unsupported or incompatible PyAR document: {path}")
    if not isinstance(value.get("workflow"), str) or not isinstance(value.get("arguments"), list):
        raise ValueError(f"Invalid PyAR document contents: {path}")
    if any(not isinstance(arg, str) for arg in value["arguments"]):
        raise ValueError(f"Invalid argument vector in {path}")
    for key in ("effective_json", "metadata_json"):
        raw = value.get(key, "{}")
        try:
            value[key.removesuffix("_json")] = json.loads(raw)
        except (TypeError, json.JSONDecodeError) as exc:
            raise ValueError(f"Invalid {key} in {path}") from exc
    if not isinstance(value.get("effective", {}), dict) or not isinstance(value.get("metadata", {}), dict):
        raise ValueError(f"Invalid settings or metadata object in {path}")
    if 'settings' in value:
        unset = value.get('unset_settings', [])
        if (not isinstance(value['settings'], dict) or not isinstance(unset, list)
                or any(not isinstance(key, str) for key in unset)):
            raise ValueError(f'Invalid settings table in {path}')
        value['effective'] = {**dict.fromkeys(unset), **value['settings']}
    value["path"] = str(path)
    return value


def profile_directory():
    root = Path(os.environ.get("XDG_CONFIG_HOME", Path.home() / ".config"))
    return root / "pyar" / "profiles"


def profile_path(name):
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", name):
        raise ValueError("Profile names may contain letters, numbers, dot, underscore, and dash")
    return profile_directory() / f"{name}.toml"


def _option_arguments(action, value):
    option = next((name for name in action.option_strings if name.startswith('--') and not name.startswith('--no-')),
                  action.option_strings[0])
    if isinstance(action, argparse.BooleanOptionalAction):
        if not isinstance(value, bool):
            raise ValueError(f"{option} must be a Boolean in saved settings")
        return [option if value else '--no-' + option[2:]]
    if isinstance(action, (argparse._StoreTrueAction, argparse._StoreFalseAction)):
        if not isinstance(value, bool):
            raise ValueError(f"{option} must be a Boolean in saved settings")
        if value == action.const:
            return [option]
        if value != action.default:
            raise ValueError(f"{option} cannot represent the saved value {value!r}")
        return []
    if isinstance(action, argparse._CountAction):
        if not isinstance(value, int) or value < 0:
            raise ValueError(f"{option} must be a nonnegative integer")
        return [option] * value
    if action.nargs in ('+', '*') or isinstance(action.nargs, int):
        if not isinstance(value, (list, tuple)):
            raise ValueError(f"{option} requires a list in saved settings")
        return [option, *(str(item) for item in value)] if value else []
    if isinstance(value, (list, tuple, dict, bool)):
        raise ValueError(f"{option} requires a scalar value in saved settings")
    # An equals token preserves values starting with '-' (e.g. custom keywords).
    return [f'{option}={value}'] if option.startswith('--') else [option, str(value)]


def arguments_from_profile(effective, parser, *, workflow, skip_defaults=True):
    """Serialize methodology options, including negative Boolean switches."""
    result = []
    excluded = PROFILE_EXCLUDED_DESTS - ({'seed'} if workflow == 'conformer' else set())
    for action in parser._actions:
        if (action.dest in excluded or action.dest not in effective or
                not action.option_strings or action.dest == 'help'):
            continue
        value = effective[action.dest]
        if value is None or (skip_defaults and value == action.default):
            continue
        result.extend(_option_arguments(action, value))
    return result


def arguments_from_values(effective, parser):
    """Replay a complete saved parser snapshot, with positional paths after --."""
    action_dests = {action.dest for action in parser._actions}
    known = action_dests | set(parser._defaults)
    unknown = set(effective) - known
    if unknown:
        raise ValueError('Unknown configuration settings: ' + ', '.join(sorted(unknown)))
    for key in set(parser._defaults) - action_dests:
        if key in effective and effective[key] != parser._defaults[key]:
            raise ValueError(f'Internal setting {key!r} is fixed for this task')
    options, positional = [], []
    for action in parser._actions:
        if action.dest in EXECUTION_CONTROL_DESTS | {'help'}:
            continue
        value = effective.get(action.dest)
        if value is None:
            continue
        if action.option_strings:
            options.extend(_option_arguments(action, value))
        else:
            positional.extend(str(item) for item in (value if isinstance(value, (list, tuple)) else [value]))
    return options + (['--', *positional] if positional else [])


def profile_overrides_from_record(record):
    """Keep methodology, including defaults that actually ran, for reproducibility."""
    workflow = record.get('workflow')
    excluded = PROFILE_EXCLUDED_DESTS - ({'seed'} if workflow == 'conformer' else set())
    return {key: value for key, value in record.get('effective', {}).items() if key not in excluded}
