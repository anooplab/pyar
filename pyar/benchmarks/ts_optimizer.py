"""Paired, manifest-driven transition-state optimizer benchmark support."""

from __future__ import annotations

import csv
import hashlib
import importlib.metadata
import json
import math
import os
import platform
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

from pyar.backends.xtb_utils import canonical_xtb_model
from pyar.data import defualt_parameters
from pyar.neb import _aligned_rmsd, read_xyz, run_neb


OPTIMIZERS = ("geometric", "sella")
SUCCESS_OUTCOMES = {"validated_first_order_saddle", "reaction_connected_success"}
SUMMARY_COLUMNS = (
    "case_id", "difficulty", "optimizer", "status", "outcome",
    "ts_optimizer_steps", "ts_backend_evaluations", "ts_wall_seconds",
    "validation_backend_evaluations", "validation_wall_seconds",
    "validation_hessian_evaluations", "validation_hessian_wall_seconds",
    "total_backend_evaluations", "reference_rmsd_angstrom",
    "reference_energy_difference_hartree", "input_sha256", "error",
)


class TSOptimizerBenchmarkError(ValueError):
    """Raised when a TS-optimizer benchmark manifest or run is invalid."""


@dataclass(frozen=True)
class TSOptimizerCase:
    """One paired TS-optimization case with a fixed starting geometry."""

    id: str
    ts_guess: str
    reference_ts: str
    reactant: str | None = None
    product: str | None = None
    charge: int = 0
    multiplicity: int = 1
    difficulty: str = "unspecified"
    reference_energy_hartree: float | None = None
    seed: int | None = None
    perturbation: dict | None = None


@dataclass(frozen=True)
class TSOptimizerBenchmark:
    """Validated benchmark specification, with paths resolved at load time."""

    name: str
    source: dict
    qc_model: dict
    settings: dict
    cases: tuple[TSOptimizerCase, ...]
    manifest_path: str
    manifest_sha256: str
    raw: dict


def _sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _require_number(value, name, *, positive=False):
    if isinstance(value, bool):
        raise TSOptimizerBenchmarkError(f"{name} must be a finite number")
    try:
        number = float(value)
    except (TypeError, ValueError, OverflowError):
        raise TSOptimizerBenchmarkError(f"{name} must be a finite number") from None
    if not math.isfinite(number) or (positive and number <= 0):
        qualifier = "positive and finite" if positive else "finite"
        raise TSOptimizerBenchmarkError(f"{name} must be {qualifier}")
    return number


def _required_text(mapping, key, context):
    value = mapping.get(key)
    if not isinstance(value, str) or not value.strip():
        raise TSOptimizerBenchmarkError(f"{context} must define non-empty {key!r}")
    return value.strip()


def _resolve_input(raw_path, root, case_id, label, *, required):
    if raw_path is None and not required:
        return None
    if not isinstance(raw_path, str) or not raw_path.strip():
        raise TSOptimizerBenchmarkError(f"case {case_id!r} must define {label!r}")
    path = Path(raw_path)
    if not path.is_absolute():
        path = root / path
    path = path.resolve()
    if not path.is_file():
        raise TSOptimizerBenchmarkError(
            f"case {case_id!r} {label} file does not exist: {path}"
        )
    return str(path)


def load_ts_optimizer_benchmark(path):
    """Load and validate a JSON manifest for paired geometric/Sella runs."""
    manifest_path = Path(path).expanduser().resolve()
    try:
        raw_bytes = manifest_path.read_bytes()
        raw = json.loads(raw_bytes)
    except OSError as exc:
        raise TSOptimizerBenchmarkError(f"Cannot read benchmark manifest {path}: {exc}") from exc
    except json.JSONDecodeError as exc:
        raise TSOptimizerBenchmarkError(f"Invalid benchmark JSON: {exc}") from exc
    if not isinstance(raw, dict):
        raise TSOptimizerBenchmarkError("benchmark manifest must contain a JSON object")

    name = _required_text(raw, "name", "benchmark")
    source = raw.get("source", {})
    qc_model = raw.get("qc_model")
    settings = raw.get("settings")
    raw_cases = raw.get("cases")
    if not isinstance(source, dict):
        raise TSOptimizerBenchmarkError("source must be an object")
    if not isinstance(qc_model, dict):
        raise TSOptimizerBenchmarkError("qc_model must be an object")
    if not isinstance(settings, dict):
        raise TSOptimizerBenchmarkError("settings must be an object")
    software = _required_text(qc_model, "software", "qc_model").lower()
    qc_model = dict(qc_model, software=software)
    for key in ("name", "version", "license"):
        _required_text(source, key, "source")
    if not any(isinstance(source.get(key), str) and source[key].strip() for key in ("doi", "url")):
        raise TSOptimizerBenchmarkError("source must define a non-empty 'doi' or 'url'")
    if software.lower() == "xtb":
        if "xtb_model" not in qc_model:
            raise TSOptimizerBenchmarkError("qc_model must explicitly define xtb_model for xTB")
        try:
            qc_model = dict(qc_model, xtb_model=canonical_xtb_model(qc_model["xtb_model"]))
        except ValueError as exc:
            raise TSOptimizerBenchmarkError(str(exc)) from exc
    elif "xtb_model" in qc_model:
        raise TSOptimizerBenchmarkError("qc_model.xtb_model is valid only when software is 'xtb'")
    if software == "xtb":
        unused = {key for key in ("method", "basis") if key in qc_model}
        if unused:
            raise TSOptimizerBenchmarkError(
                "qc_model.method and qc_model.basis are not used by software='xtb'; "
                "select its Hamiltonian with xtb_model"
            )
    else:
        qc_model.setdefault("method", defualt_parameters.values["method"])
        qc_model.setdefault("basis", defualt_parameters.values["basis"])
    qc_model.setdefault("nprocs", 1)
    for key in ("method", "basis"):
        if key in qc_model and (not isinstance(qc_model[key], str) or not qc_model[key].strip()):
            raise TSOptimizerBenchmarkError(f"qc_model.{key} must be a non-empty string when provided")
    unsupported_qc = set(qc_model) - {"software", "method", "basis", "xtb_model", "nprocs"}
    if unsupported_qc:
        raise TSOptimizerBenchmarkError(
            "unsupported qc_model key(s): " + ", ".join(sorted(unsupported_qc))
        )
    if "nprocs" in qc_model and (
        isinstance(qc_model["nprocs"], bool) or not isinstance(qc_model["nprocs"], int)
        or qc_model["nprocs"] < 1
    ):
        raise TSOptimizerBenchmarkError("qc_model.nprocs must be a positive integer")

    settings = dict(settings)
    allowed_settings = {
        "ts_fmax", "ts_max_cycles", "imaginary_frequency_threshold", "irc_max_cycles",
        "endpoint_max_cycles", "product_relaxation_fmax", "product_relaxation_max_steps",
        "irc_endpoint_rmsd_tolerance", "sella_internal_coordinates",
    }
    unsupported_settings = set(settings) - allowed_settings
    if unsupported_settings:
        raise TSOptimizerBenchmarkError(
            "unsupported settings key(s): " + ", ".join(sorted(unsupported_settings))
        )
    if "ts_fmax" not in settings or "ts_max_cycles" not in settings:
        raise TSOptimizerBenchmarkError(
            "settings must define common ts_fmax and ts_max_cycles for paired runs"
        )
    settings["ts_fmax"] = _require_number(settings["ts_fmax"], "settings.ts_fmax", positive=True)
    cycles = settings["ts_max_cycles"]
    if isinstance(cycles, bool) or not isinstance(cycles, int) or cycles < 1:
        raise TSOptimizerBenchmarkError("settings.ts_max_cycles must be a positive integer")
    if "sella_internal_coordinates" in settings and not isinstance(
        settings["sella_internal_coordinates"], bool
    ):
        raise TSOptimizerBenchmarkError("settings.sella_internal_coordinates must be a boolean")
    settings.setdefault("imaginary_frequency_threshold", 20.0)
    settings["imaginary_frequency_threshold"] = _require_number(
        settings["imaginary_frequency_threshold"],
        "settings.imaginary_frequency_threshold", positive=True,
    )
    settings.setdefault("irc_max_cycles", 200)
    settings.setdefault("endpoint_max_cycles", 300)
    settings.setdefault("product_relaxation_fmax", 0.05)
    settings.setdefault("product_relaxation_max_steps", 200)
    settings.setdefault("irc_endpoint_rmsd_tolerance", 0.5)
    for key in ("irc_max_cycles", "endpoint_max_cycles", "product_relaxation_max_steps"):
        if isinstance(settings[key], bool) or not isinstance(settings[key], int) or settings[key] < 1:
            raise TSOptimizerBenchmarkError(f"settings.{key} must be a positive integer")
    for key in ("product_relaxation_fmax", "irc_endpoint_rmsd_tolerance"):
        settings[key] = _require_number(settings[key], f"settings.{key}", positive=True)

    if not isinstance(raw_cases, list) or not raw_cases:
        raise TSOptimizerBenchmarkError("benchmark manifest must define a non-empty cases list")
    root = manifest_path.parent
    cases = []
    seen_ids = set()
    for index, raw_case in enumerate(raw_cases):
        if not isinstance(raw_case, dict):
            raise TSOptimizerBenchmarkError(f"case {index} must be an object")
        case_id = _required_text(raw_case, "id", f"case {index}")
        if case_id in seen_ids:
            raise TSOptimizerBenchmarkError(f"duplicate case id: {case_id}")
        if case_id in {".", ".."} or any(char not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_.-" for char in case_id):
            raise TSOptimizerBenchmarkError(
                f"case id {case_id!r} may contain only letters, numbers, period, underscore, and hyphen"
            )
        seen_ids.add(case_id)
        reactant_raw, product_raw = raw_case.get("reactant"), raw_case.get("product")
        if (reactant_raw is None) != (product_raw is None):
            raise TSOptimizerBenchmarkError(
                f"case {case_id!r} must provide both reactant and product endpoints, or neither"
            )
        charge = raw_case.get("charge", 0)
        multiplicity = raw_case.get("multiplicity", 1)
        if isinstance(charge, bool) or not isinstance(charge, int):
            raise TSOptimizerBenchmarkError(f"case {case_id!r} charge must be an integer")
        if isinstance(multiplicity, bool) or not isinstance(multiplicity, int) or multiplicity < 1:
            raise TSOptimizerBenchmarkError(f"case {case_id!r} multiplicity must be a positive integer")
        reference_energy = raw_case.get("reference_energy_hartree")
        if reference_energy is not None:
            reference_energy = _require_number(
                reference_energy, f"case {case_id!r} reference_energy_hartree"
            )
        seed = raw_case.get("seed")
        if seed is not None and (isinstance(seed, bool) or not isinstance(seed, int)):
            raise TSOptimizerBenchmarkError(f"case {case_id!r} seed must be an integer")
        perturbation = raw_case.get("perturbation")
        if perturbation is not None and not isinstance(perturbation, dict):
            raise TSOptimizerBenchmarkError(f"case {case_id!r} perturbation must be an object")
        cases.append(TSOptimizerCase(
            id=case_id,
            ts_guess=_resolve_input(raw_case.get("ts_guess"), root, case_id, "ts_guess", required=True),
            reference_ts=_resolve_input(raw_case.get("reference_ts"), root, case_id, "reference_ts", required=True),
            reactant=_resolve_input(reactant_raw, root, case_id, "reactant", required=False),
            product=_resolve_input(product_raw, root, case_id, "product", required=False),
            charge=charge,
            multiplicity=multiplicity,
            difficulty=str(raw_case.get("difficulty", "unspecified") or "unspecified"),
            reference_energy_hartree=reference_energy,
            seed=seed,
            perturbation=perturbation,
        ))
        case = cases[-1]
        ordered_symbols, _, _ = read_xyz(case.ts_guess)
        for label in ("reference_ts", "reactant", "product"):
            geometry_path = getattr(case, label)
            if geometry_path is not None:
                symbols, _, _ = read_xyz(geometry_path)
                if symbols != ordered_symbols:
                    raise TSOptimizerBenchmarkError(
                        f"case {case.id!r} {label} must preserve the ordered atoms of ts_guess"
                    )

    return TSOptimizerBenchmark(
        name=name, source=dict(source), qc_model=dict(qc_model), settings=settings,
        cases=tuple(cases), manifest_path=str(manifest_path),
        manifest_sha256=hashlib.sha256(raw_bytes).hexdigest(), raw=raw,
    )


def _metadata_version(distribution):
    try:
        return importlib.metadata.version(distribution)
    except importlib.metadata.PackageNotFoundError:
        return None


def _backend_version(qc_model):
    """Return a best-effort version string without making it a run prerequisite."""
    software = str(qc_model["software"]).lower()
    executable = "xtb" if software == "xtb" else software
    executable_path = shutil.which(executable)
    if executable_path is None:
        return {"executable": executable, "path": None, "version": None}
    version_text = None
    if software == "xtb":
        try:
            process = subprocess.run(
                [executable_path, "--version"], check=False, capture_output=True,
                text=True, timeout=15,
            )
            version_text = (process.stdout or process.stderr).strip() or None
        except (OSError, subprocess.TimeoutExpired):
            pass
    return {"executable": executable, "path": executable_path, "version": version_text}


def _git_revision():
    repository = Path(__file__).resolve().parents[2]
    try:
        completed = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=repository, check=False,
            capture_output=True, text=True, timeout=5,
        )
    except (OSError, subprocess.TimeoutExpired):
        return None
    return completed.stdout.strip() if completed.returncode == 0 else None


def _runtime_payload(qc_model):
    """Capture the environment where a job actually executes."""
    return {
        "pyar_version": _metadata_version("pyar-chem"),
        "python": sys.version,
        "platform": platform.platform(),
        "git_commit": _git_revision(),
        "backend": _backend_version(qc_model),
        "geometric_version": _metadata_version("geometric"),
        "sella_version": _metadata_version("sella"),
    }


def _json_safe(value):
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def _write_json(path, payload):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w", encoding="utf-8", dir=path.parent,
        prefix=f".{path.name}.", suffix=".tmp", delete=False,
    ) as stream:
        temporary = Path(stream.name)
        stream.write(json.dumps(_json_safe(payload), indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)


def _run_manifest_payload(spec):
    expected = []
    for case in spec.cases:
        input_hashes = {"ts_guess": _sha256(case.ts_guess), "reference_ts": _sha256(case.reference_ts)}
        guess_symbols, _, _ = read_xyz(case.ts_guess)
        for label in ("reference_ts", "reactant", "product"):
            candidate_path = getattr(case, label)
            if candidate_path is not None:
                symbols, _, _ = read_xyz(candidate_path)
                if symbols != guess_symbols:
                    raise TSOptimizerBenchmarkError(
                        f"case {case.id!r} {label} must preserve the ordered atoms of ts_guess"
                    )
        if case.reactant:
            input_hashes["reactant"] = _sha256(case.reactant)
            input_hashes["product"] = _sha256(case.product)
        expected.append({
            "case_id": case.id,
            "difficulty": case.difficulty,
            "ts_guess_sha256": input_hashes["ts_guess"],
            "input_hashes": input_hashes,
            "charge": case.charge,
            "multiplicity": case.multiplicity,
        })
    return {
        "schema_version": 1,
        "benchmark_name": spec.name,
        "manifest_path": spec.manifest_path,
        "manifest_sha256": spec.manifest_sha256,
        "source": spec.source,
        "qc_model": spec.qc_model,
        "settings": spec.settings,
        "expected_optimizers": list(OPTIMIZERS),
        "expected_runs": expected,
        "runtime": _runtime_payload(spec.qc_model),
    }


def _prepare_run_root(spec, output):
    output = Path(output).expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    run_manifest_path = output / "run_manifest.json"
    if run_manifest_path.exists():
        existing = json.loads(run_manifest_path.read_text(encoding="utf-8"))
        if existing.get("manifest_sha256") != spec.manifest_sha256:
            raise TSOptimizerBenchmarkError(
                "output directory belongs to a different benchmark manifest"
            )
        expected_case_ids = [case.id for case in spec.cases]
        stored_case_ids = [item.get("case_id") for item in existing.get("expected_runs", [])]
        if (
            existing.get("benchmark_name") != spec.name
            or existing.get("source") != spec.source
            or existing.get("qc_model") != spec.qc_model
            or existing.get("settings") != spec.settings
            or existing.get("expected_optimizers") != list(OPTIMIZERS)
            or stored_case_ids != expected_case_ids
        ):
            raise TSOptimizerBenchmarkError("run_manifest.json does not match the benchmark manifest")
        source_manifest_copy = output / "benchmark_manifest.json"
        if not source_manifest_copy.is_file() or _sha256(source_manifest_copy) != spec.manifest_sha256:
            raise TSOptimizerBenchmarkError(
                "benchmark_manifest.json is missing or differs from the original manifest"
            )
    else:
        shutil.copy2(spec.manifest_path, output / "benchmark_manifest.json")
        _write_json(run_manifest_path, _run_manifest_payload(spec))
    return output


VALIDATION_SETTINGS = (
    "imaginary_frequency_threshold", "product_relaxation_fmax",
    "product_relaxation_max_steps", "irc_max_cycles", "endpoint_max_cycles",
    "irc_endpoint_rmsd_tolerance",
)


def _stage_arguments(spec, case, output):
    args = {
        "software": spec.qc_model["software"],
        "charge": case.charge,
        "multiplicity": case.multiplicity,
        "nprocs": spec.qc_model.get("nprocs", 1),
        "output": output,
        # Later stages recursively verify their dependencies against these
        # settings. Carry the same protocol through every stage, including IRC.
        **{key: spec.settings[key] for key in VALIDATION_SETTINGS},
    }
    for key in ("method", "basis"):
        if key in spec.qc_model:
            args[key] = spec.qc_model[key]
    if spec.settings.get("sella_internal_coordinates") is True:
        args["sella_internal_coordinates"] = True
    if str(spec.qc_model["software"]).lower() == "xtb":
        args["xtb_model"] = spec.qc_model["xtb_model"]
    return args


def _stage_counts(run_directory, names):
    evaluation_count = 0
    wall_seconds = 0.0
    found = False
    for name in names:
        path = Path(run_directory) / f"{name}_summary.json"
        if not path.is_file():
            continue
        try:
            result = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        if isinstance(result.get("backend_energy_gradient_evaluations"), int):
            evaluation_count += result["backend_energy_gradient_evaluations"]
            found = True
        if isinstance(result.get("wall_seconds"), (float, int)):
            wall_seconds += float(result["wall_seconds"])
    return evaluation_count if found else None, wall_seconds if found else None


def _hessian_metrics(run_directory, stage_names):
    """Collect Hessian work without inferring it from total provider calls."""
    evaluations, wall_seconds, sources = [], [], set()
    unknown_evaluations = unknown_wall_time = False
    for stage in stage_names:
        summary = _read_summary(run_directory, stage)
        if not summary:
            continue
        candidates = [summary]
        if stage == "endpoints":
            candidates.extend(
                value for key, value in summary.items()
                if key.endswith("_frequency") and isinstance(value, dict)
            )
        for candidate in candidates:
            source = candidate.get("hessian_source")
            if source:
                sources.add(str(source))
            count = candidate.get("hessian_evaluations")
            elapsed = candidate.get("hessian_wall_seconds")
            if source and count is None:
                unknown_evaluations = True
            elif isinstance(count, int):
                evaluations.append(count)
            if source and elapsed is None:
                unknown_wall_time = True
            elif isinstance(elapsed, (int, float)):
                wall_seconds.append(float(elapsed))
    return {
        "evaluations": sum(evaluations) if evaluations and not unknown_evaluations else None,
        "wall_seconds": sum(wall_seconds) if wall_seconds and not unknown_wall_time else None,
        "sources": sorted(sources),
    }


def _run_ts_optimizer_case_locked(benchmark, *, case_id, optimizer, output):
    """Run one optimizer for one manifest case and retain all workflow files."""
    spec = benchmark if isinstance(benchmark, TSOptimizerBenchmark) else load_ts_optimizer_benchmark(benchmark)
    if optimizer not in OPTIMIZERS:
        raise TSOptimizerBenchmarkError(
            f"optimizer must be one of {', '.join(OPTIMIZERS)}; no optimizer fallback is performed"
        )
    case = next((entry for entry in spec.cases if entry.id == case_id), None)
    if case is None:
        raise TSOptimizerBenchmarkError(f"unknown case id: {case_id}")
    output_root = _prepare_run_root(spec, output)
    run_directory = output_root / "cases" / case.id / optimizer
    result_path = run_directory / "benchmark_result.json"
    if result_path.exists():
        raise TSOptimizerBenchmarkError(
            f"result already exists at {result_path}; choose a fresh output directory to preserve raw runs"
        )
    existing_run_files = [
        path for path in run_directory.iterdir()
        if path.name != ".benchmark-job.lock"
    ]
    if existing_run_files:
        raise TSOptimizerBenchmarkError(
            f"run directory already contains partial or uncollected artifacts: {run_directory}; "
            "choose a fresh output directory to preserve raw runs"
        )
    inputs = {"ts_guess": _sha256(case.ts_guess), "reference_ts": _sha256(case.reference_ts)}
    if case.reactant:
        inputs["reactant"] = _sha256(case.reactant)
        inputs["product"] = _sha256(case.product)
    run_manifest = json.loads((output_root / "run_manifest.json").read_text(encoding="utf-8"))
    expected_run = next(item for item in run_manifest["expected_runs"] if item["case_id"] == case.id)
    if inputs != expected_run["input_hashes"]:
        raise TSOptimizerBenchmarkError(
            f"case {case.id!r} input files changed after the benchmark run manifest was created; "
            "choose a fresh output directory"
        )
    run_directory.mkdir(parents=True, exist_ok=True)
    input_paths = {}
    for label, input_path in (("ts_guess", case.ts_guess), ("reference_ts", case.reference_ts),
                              ("reactant", case.reactant), ("product", case.product)):
        if input_path is not None:
            copied_path = run_directory / f"input_{label}.xyz"
            shutil.copy2(input_path, copied_path)
            input_paths[label] = copied_path
            if _sha256(copied_path) != inputs[label]:
                raise TSOptimizerBenchmarkError(
                    f"case {case.id!r} input {label} changed while being copied; choose a fresh output directory"
                )

    guess_symbols, _, _ = read_xyz(input_paths["ts_guess"])
    reference_symbols, reference_coordinates, _ = read_xyz(input_paths["reference_ts"])
    if guess_symbols != reference_symbols:
        raise TSOptimizerBenchmarkError(
            f"case {case.id!r} reference_ts atom order differs from ts_guess"
        )

    base = _stage_arguments(spec, case, run_directory)
    runtime = _runtime_payload(spec.qc_model)
    settings = spec.settings
    stages = {}
    outcome = "optimizer_exception"
    failure_stage = "ts"
    error = None
    ts_result = None
    frequency_result = None
    endpoint_result = None
    try:
        ts_result = run_neb(
            ts_geometry=input_paths["ts_guess"], stage="ts", ts_optimizer=optimizer,
            ts_fmax=settings["ts_fmax"], ts_max_cycles=settings["ts_max_cycles"],
            **base,
        )
        stages["ts"] = ts_result
        if ts_result.get("ts_optimization_converged") is not True:
            outcome = "optimizer_not_converged"
        else:
            failure_stage = "frequency"
            frequency_result = run_neb(
                ts_geometry=run_directory / "ts_optimized.xyz", stage="frequency",
                **base,
            )
            stages["frequency"] = frequency_result
            if frequency_result.get("stationary") is not True:
                outcome = "converged_not_stationary"
            elif frequency_result.get("first_order_saddle_confirmed") is not True:
                outcome = "stationary_not_first_order_saddle"
            elif case.reactant is None:
                outcome = "validated_first_order_saddle"
            else:
                failure_stage = "relax"
                run_neb(
                    start=input_paths["reactant"], end=input_paths["product"], stage="relax",
                    **base,
                )
                stages["relax"] = json.loads(
                    (run_directory / "relax_summary.json").read_text(encoding="utf-8")
                )
                failure_stage = "irc"
                run_neb(stage="irc", **base)
                stages["irc"] = json.loads(
                    (run_directory / "irc_summary.json").read_text(encoding="utf-8")
                )
                failure_stage = "endpoints"
                endpoint_result = run_neb(
                    start=input_paths["reactant"], end=input_paths["product"], stage="endpoints",
                    **base,
                )
                stages["endpoints"] = endpoint_result
                outcome = (
                    "reaction_connected_success"
                    if endpoint_result.get("reactant_product_connection_confirmed") is True
                    else "first_order_saddle_wrong_connection"
                )
    except Exception as exc:  # Keep the failed job as an analyzable result.
        error = f"{type(exc).__name__}: {exc}"
        outcome = {
            "ts": "optimizer_exception",
            "frequency": "frequency_exception",
            "relax": "endpoint_relaxation_exception",
            "irc": "irc_exception",
            "endpoints": "endpoint_validation_exception",
        }.get(failure_stage, "optimizer_exception")

    ts_summary = _read_summary(run_directory, "ts")
    reference_rmsd = None
    energy_difference = None
    final_geometry = run_directory / "ts_optimized.xyz"
    if ts_summary and final_geometry.is_file():
        try:
            final_symbols, final_coordinates, _ = read_xyz(final_geometry)
            if final_symbols != reference_symbols:
                raise ValueError("optimized geometry atom identities/order differ from reference")
            reference_rmsd = _aligned_rmsd(reference_coordinates, final_coordinates)
        except (OSError, ValueError) as exc:
            if error is None:
                error = f"Invalid optimized geometry: {exc}"
                outcome, failure_stage = "optimizer_exception", "ts"
        final_energy = ts_summary.get("ts_energy_hartree")
        if case.reference_energy_hartree is not None and isinstance(final_energy, (int, float)):
            energy_difference = float(final_energy) - case.reference_energy_hartree

    ts_evaluations, _ = _stage_counts(run_directory, ("ts",))
    validation_evaluations, validation_seconds = _stage_counts(
        run_directory, ("frequency", "relax", "irc", "endpoints")
    )
    hessian_metrics = _hessian_metrics(run_directory, ("frequency", "endpoints"))
    result = {
        "schema_version": 1,
        "benchmark_name": spec.name,
        "case_id": case.id,
        "difficulty": case.difficulty,
        "seed": case.seed,
        "perturbation": case.perturbation,
        "optimizer": optimizer,
        "optimizer_class": "transition_state",
        "status": "failed" if error is not None else "complete",
        "outcome": outcome,
        "failure_stage": None if error is None else failure_stage,
        "error": error,
        "source": spec.source,
        "qc_model": spec.qc_model,
        "case_qc_settings": {"charge": case.charge, "multiplicity": case.multiplicity},
        "settings": settings,
        "runtime": runtime,
        "input_hashes": inputs,
        "input_sha256": inputs["ts_guess"],
        "reference_rmsd_angstrom": reference_rmsd,
        "reference_energy_difference_hartree": energy_difference,
        "ts_optimizer_steps": None if ts_summary is None else ts_summary.get("optimizer_steps"),
        "ts_backend_energy_gradient_evaluations": ts_evaluations,
        "ts_wall_seconds": None if ts_summary is None else ts_summary.get("ts_optimization_wall_seconds"),
        "validation_backend_energy_gradient_evaluations": validation_evaluations,
        "validation_wall_seconds": validation_seconds,
        "validation_hessian_evaluations": hessian_metrics["evaluations"],
        "validation_hessian_wall_seconds": hessian_metrics["wall_seconds"],
        "hessian_sources": hessian_metrics["sources"],
        "total_backend_energy_gradient_evaluations": (
            None if ts_evaluations is None and validation_evaluations is None
            else (ts_evaluations or 0) + (validation_evaluations or 0)
        ),
        "ts_result": ts_result,
        "frequency_result": frequency_result,
        "endpoint_result": endpoint_result,
        "stages": list(stages),
        "run_directory": str(run_directory),
        "finished_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    _write_json(result_path, result)
    return result


def run_ts_optimizer_case(benchmark, *, case_id, optimizer, output):
    """Reserve one case/optimizer pair before running, safe for array jobs."""
    spec = benchmark if isinstance(benchmark, TSOptimizerBenchmark) else load_ts_optimizer_benchmark(benchmark)
    if optimizer not in OPTIMIZERS:
        raise TSOptimizerBenchmarkError(
            f"optimizer must be one of {', '.join(OPTIMIZERS)}; no optimizer fallback is performed"
        )
    if not any(case.id == case_id for case in spec.cases):
        raise TSOptimizerBenchmarkError(f"unknown case id: {case_id}")
    output_root = _prepare_run_root(spec, output)
    run_directory = output_root / "cases" / case_id / optimizer
    run_directory.mkdir(parents=True, exist_ok=True)
    lock = run_directory / ".benchmark-job.lock"
    try:
        descriptor = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    except FileExistsError:
        raise TSOptimizerBenchmarkError(
            f"job is already running for case {case_id!r} with optimizer {optimizer}"
        ) from None
    os.close(descriptor)
    try:
        return _run_ts_optimizer_case_locked(
            spec, case_id=case_id, optimizer=optimizer, output=output_root,
        )
    finally:
        lock.unlink(missing_ok=True)


def _read_summary(run_directory, stage):
    path = Path(run_directory) / f"{stage}_summary.json"
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None


def _result_row(result):
    return {
        "case_id": result.get("case_id"),
        "difficulty": result.get("difficulty"),
        "optimizer": result.get("optimizer"),
        "status": result.get("status"),
        "outcome": result.get("outcome"),
        "ts_optimizer_steps": result.get("ts_optimizer_steps"),
        "ts_backend_evaluations": result.get("ts_backend_energy_gradient_evaluations"),
        "ts_wall_seconds": result.get("ts_wall_seconds"),
        "validation_backend_evaluations": result.get("validation_backend_energy_gradient_evaluations"),
        "validation_wall_seconds": result.get("validation_wall_seconds"),
        "validation_hessian_evaluations": result.get("validation_hessian_evaluations"),
        "validation_hessian_wall_seconds": result.get("validation_hessian_wall_seconds"),
        "total_backend_evaluations": result.get("total_backend_energy_gradient_evaluations"),
        "reference_rmsd_angstrom": result.get("reference_rmsd_angstrom"),
        "reference_energy_difference_hartree": result.get("reference_energy_difference_hartree"),
        "input_sha256": result.get("input_sha256"),
        "error": result.get("error"),
    }


def _write_csv(path, columns, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def collect_ts_optimizer_benchmark(output):
    """Collect complete and missing paired jobs into machine-readable tables."""
    output = Path(output).expanduser().resolve()
    manifest_path = output / "run_manifest.json"
    if not manifest_path.is_file():
        raise TSOptimizerBenchmarkError(f"missing {manifest_path.name} in {output}")
    run_manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    expected_runs = run_manifest.get("expected_runs")
    optimizers = run_manifest.get("expected_optimizers")
    if not isinstance(expected_runs, list) or not isinstance(optimizers, list):
        raise TSOptimizerBenchmarkError("run_manifest.json lacks the expected run matrix")
    result_map = {}
    rows = []
    for expected in expected_runs:
        case_id = expected["case_id"]
        for optimizer in optimizers:
            result_path = output / "cases" / case_id / optimizer / "benchmark_result.json"
            key = (case_id, optimizer)
            if result_path.is_file():
                try:
                    result = json.loads(result_path.read_text(encoding="utf-8"))
                except json.JSONDecodeError as exc:
                    result = {
                        "case_id": case_id, "optimizer": optimizer, "difficulty": expected["difficulty"],
                        "status": "invalid_result", "outcome": "invalid_result",
                        "error": f"invalid result JSON: {exc}", "input_sha256": expected["ts_guess_sha256"],
                    }
            else:
                result = {
                    "case_id": case_id, "optimizer": optimizer, "difficulty": expected["difficulty"],
                    "status": "incomplete", "outcome": "incomplete_run",
                    "error": None, "input_sha256": expected["ts_guess_sha256"],
                }
            result_map[key] = result
            rows.append(_result_row(result))
    paired_rows = []
    by_optimizer = {optimizer: [] for optimizer in optimizers}
    for row in rows:
        by_optimizer.setdefault(row["optimizer"], []).append(row)
    for expected in expected_runs:
        case_id = expected["case_id"]
        geometric = result_map.get((case_id, "geometric"))
        sella = result_map.get((case_id, "sella"))
        expected_qc = {"charge": expected["charge"], "multiplicity": expected["multiplicity"]}

        def identity_matches_manifest(result, optimizer):
            return bool(
                result
                and result.get("schema_version") == 1
                and result.get("case_id") == case_id
                and result.get("optimizer") == optimizer
                and result.get("source") == run_manifest.get("source")
                and result.get("input_hashes") == expected["input_hashes"]
                and result.get("qc_model") == run_manifest.get("qc_model")
                and result.get("case_qc_settings") == expected_qc
                and result.get("settings") == run_manifest.get("settings")
                and result.get("runtime") == run_manifest.get("runtime")
            )

        identity_matches = bool(
            identity_matches_manifest(geometric, "geometric")
            and identity_matches_manifest(sella, "sella")
        )
        g_success = bool(geometric and geometric.get("outcome") in SUCCESS_OUTCOMES)
        s_success = bool(sella and sella.get("outcome") in SUCCESS_OUTCOMES)
        g_evals = None if not geometric else geometric.get("total_backend_energy_gradient_evaluations")
        s_evals = None if not sella else sella.get("total_backend_energy_gradient_evaluations")
        g_wall = None if not geometric else geometric.get("ts_wall_seconds")
        s_wall = None if not sella else sella.get("ts_wall_seconds")
        paired_rows.append({
            "case_id": case_id,
            "difficulty": expected["difficulty"],
            "input_sha256": expected["ts_guess_sha256"],
            "pair_identity_matches": identity_matches,
            "geometric_outcome": None if not geometric else geometric.get("outcome"),
            "sella_outcome": None if not sella else sella.get("outcome"),
            "geometric_success": g_success,
            "sella_success": s_success,
            "discordant_success": g_success != s_success if identity_matches else None,
            "total_evaluation_difference_geometric_minus_sella": (
                g_evals - s_evals if identity_matches and isinstance(g_evals, int)
                and isinstance(s_evals, int) else None
            ),
            "ts_wall_difference_seconds_geometric_minus_sella": (
                float(g_wall) - float(s_wall) if identity_matches
                and isinstance(g_wall, (int, float)) and isinstance(s_wall, (int, float)) else None
            ),
            "ts_evaluation_difference_geometric_minus_sella": (
                geometric.get("ts_backend_energy_gradient_evaluations")
                - sella.get("ts_backend_energy_gradient_evaluations")
                if identity_matches and isinstance(geometric.get("ts_backend_energy_gradient_evaluations"), int)
                and isinstance(sella.get("ts_backend_energy_gradient_evaluations"), int) else None
            ),
            "validation_evaluation_difference_geometric_minus_sella": (
                geometric.get("validation_backend_energy_gradient_evaluations")
                - sella.get("validation_backend_energy_gradient_evaluations")
                if identity_matches and isinstance(geometric.get("validation_backend_energy_gradient_evaluations"), int)
                and isinstance(sella.get("validation_backend_energy_gradient_evaluations"), int) else None
            ),
            "validation_wall_difference_seconds_geometric_minus_sella": (
                float(geometric.get("validation_wall_seconds"))
                - float(sella.get("validation_wall_seconds"))
                if identity_matches
                and isinstance(geometric.get("validation_wall_seconds"), (int, float))
                and isinstance(sella.get("validation_wall_seconds"), (int, float)) else None
            ),
            "validation_hessian_evaluation_difference_geometric_minus_sella": (
                geometric.get("validation_hessian_evaluations")
                - sella.get("validation_hessian_evaluations")
                if identity_matches and isinstance(geometric.get("validation_hessian_evaluations"), int)
                and isinstance(sella.get("validation_hessian_evaluations"), int) else None
            ),
            "validation_hessian_wall_difference_seconds_geometric_minus_sella": (
                float(geometric.get("validation_hessian_wall_seconds"))
                - float(sella.get("validation_hessian_wall_seconds"))
                if identity_matches
                and isinstance(geometric.get("validation_hessian_wall_seconds"), (int, float))
                and isinstance(sella.get("validation_hessian_wall_seconds"), (int, float)) else None
            ),
            "pair_status": "paired" if identity_matches else "incomplete_or_mismatched",
        })
    _write_csv(output / "runs.csv", SUMMARY_COLUMNS, rows)
    _write_csv(output / "paired.csv", tuple(paired_rows[0]) if paired_rows else (
        "case_id", "difficulty", "input_sha256", "pair_identity_matches", "geometric_outcome",
        "sella_outcome", "geometric_success", "sella_success", "discordant_success",
        "total_evaluation_difference_geometric_minus_sella",
        "ts_wall_difference_seconds_geometric_minus_sella",
        "ts_evaluation_difference_geometric_minus_sella",
        "validation_evaluation_difference_geometric_minus_sella",
        "validation_wall_difference_seconds_geometric_minus_sella",
        "validation_hessian_evaluation_difference_geometric_minus_sella",
        "validation_hessian_wall_difference_seconds_geometric_minus_sella", "pair_status",
    ), paired_rows)

    summaries = {}
    for optimizer, optimizer_rows in by_optimizer.items():
        outcomes = {}
        for row in optimizer_rows:
            outcomes[row["outcome"]] = outcomes.get(row["outcome"], 0) + 1
        successes = sum(
            row["outcome"] in SUCCESS_OUTCOMES for row in optimizer_rows
        )
        completed = sum(row["status"] in {"complete", "failed"} for row in optimizer_rows)
        evaluations = [row["total_backend_evaluations"] for row in optimizer_rows
                       if isinstance(row["total_backend_evaluations"], int)]
        timings = [row["ts_wall_seconds"] for row in optimizer_rows
                   if isinstance(row["ts_wall_seconds"], (int, float))]
        summaries[optimizer] = {
            "expected_runs": len(expected_runs),
            "completed_runs": completed,
            "successful_runs": successes,
            "success_fraction_of_completed": successes / completed if completed else None,
            "outcomes": outcomes,
            "mean_total_backend_evaluations": (
                sum(evaluations) / len(evaluations) if evaluations else None
            ),
            "mean_ts_wall_seconds": sum(timings) / len(timings) if timings else None,
        }
    summary = {
        "schema_version": 1,
        "benchmark_name": run_manifest.get("benchmark_name"),
        "manifest_sha256": run_manifest.get("manifest_sha256"),
        "expected_runs": len(expected_runs) * len(optimizers),
        "completed_runs": sum(value["completed_runs"] for value in summaries.values()),
        "incomplete_runs": sum(
            row["outcome"] == "incomplete_run" for row in rows
        ),
        "optimizers": summaries,
        "paired_cases": sum(row["pair_status"] == "paired" for row in paired_rows),
        "paired": paired_rows,
    }
    _write_json(output / "runs.json", rows)
    _write_json(output / "summary.json", summary)
    return summary
