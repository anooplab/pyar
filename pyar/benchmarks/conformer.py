"""Conformer-search benchmark and failure diagnosis utilities."""

from __future__ import annotations

import csv
import json
import math
import shutil
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from pyar.core.molecule import Molecule
from pyar.data import defualt_parameters
from pyar.workflow_results import ConformerResult
from pyar.workflows import conformer as conformer_workflow


FAILURE_CLASSES = {
    "found": "reference/global minimum found in selected conformers",
    "generated_lost_selection": "reference-like conformer generated but not selected",
    "generated_lost_backend_or_dedup": "reference-like conformer generated but lost after backend refinement or deduplication",
    "never_generated": "reference-like conformer was not generated",
    "wrong_ranking": "reference-like conformer generated but ranked outside the requested energy window",
    "input_chemistry_issue": "input or reference chemistry prevents a meaningful comparison",
    "uncertain": "available data are insufficient for a conservative diagnosis",
}

SUMMARY_COLUMNS = [
    "case_id",
    "input",
    "input_format",
    "status",
    "failure_class",
    "reference_found_generated",
    "reference_found_selected",
    "best_generated_rmsd",
    "best_selected_rmsd",
    "best_generated_rank",
    "best_selected_rank",
    "lowest_rdkit_energy",
    "lowest_backend_energy",
    "selected_count",
    "generated_count",
    "backend_requested",
    "notes",
]


class ConformerBenchmarkError(ValueError):
    """Raised when a conformer benchmark specification is invalid."""


@dataclass(frozen=True)
class ConformerBenchmarkCase:
    """One molecule in a conformer benchmark specification."""

    id: str
    input: str
    input_format: str = "auto"
    charge: int | None = None
    multiplicity: int = 1
    reference_xyz: str | None = None
    reference_energy: float | None = None
    notes: str = ""


@dataclass(frozen=True)
class ConformerBenchmarkSpec:
    """A loaded conformer benchmark specification."""

    name: str
    method: str
    cases: tuple[ConformerBenchmarkCase, ...]
    path: str | None = None


def _optional_int(value):
    if value in {None, ""}:
        return None
    return int(value)


def _optional_float(value):
    if value in {None, ""}:
        return None
    return float(value)


def _require_string(mapping, key, context):
    value = mapping.get(key)
    if not isinstance(value, str) or not value.strip():
        raise ConformerBenchmarkError(f"{context} must define non-empty {key!r}")
    return value.strip()


def _load_case(raw_case, index):
    if not isinstance(raw_case, dict):
        raise ConformerBenchmarkError(f"case {index} must be an object")
    return ConformerBenchmarkCase(
        id=_require_string(raw_case, "id", f"case {index}"),
        input=_require_string(raw_case, "input", f"case {index}"),
        input_format=str(raw_case.get("input_format", "auto") or "auto"),
        charge=_optional_int(raw_case.get("charge")),
        multiplicity=int(raw_case.get("multiplicity", 1) or 1),
        reference_xyz=raw_case.get("reference_xyz"),
        reference_energy=_optional_float(raw_case.get("reference_energy")),
        notes=str(raw_case.get("notes", "") or ""),
    )


def load_conformer_benchmark(path):
    """Load and validate a JSON conformer benchmark specification."""
    benchmark_path = Path(path)
    with benchmark_path.open(encoding="utf-8") as fp:
        raw = json.load(fp)
    if not isinstance(raw, dict):
        raise ConformerBenchmarkError("benchmark file must contain a JSON object")
    raw_cases = raw.get("cases")
    if not isinstance(raw_cases, list):
        raise ConformerBenchmarkError("benchmark file must define a cases list")
    cases = tuple(_load_case(raw_case, index) for index, raw_case in enumerate(raw_cases))
    return ConformerBenchmarkSpec(
        name=_require_string(raw, "name", "benchmark"),
        method=str(raw.get("method", "reference") or "reference"),
        cases=cases,
        path=str(benchmark_path),
    )


def _resolve_case_path(case, benchmark_root):
    if case.reference_xyz is None:
        return None
    path = Path(case.reference_xyz)
    if not path.is_absolute():
        path = benchmark_root / path
    return path


def _json_safe(value):
    if isinstance(value, np.generic):
        return value.item()
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
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fp:
        json.dump(_json_safe(payload), fp, indent=2, sort_keys=True)
        fp.write("\n")


def _write_csv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=SUMMARY_COLUMNS)
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column) for column in SUMMARY_COLUMNS})


def _heavy_coordinates(molecule):
    indices = [
        index
        for index, atom in enumerate(molecule.atoms_list)
        if str(atom).upper() != "H"
    ]
    if not indices:
        indices = list(range(len(molecule.atoms_list)))
    if not indices or molecule.coordinates is None:
        return None
    return np.asarray(molecule.coordinates, dtype=float)[np.asarray(indices, dtype=int)]


def _heavy_symbols(molecule):
    symbols = [str(atom).capitalize() for atom in molecule.atoms_list if str(atom).upper() != "H"]
    return symbols or [str(atom).capitalize() for atom in molecule.atoms_list]


def _kabsch_rmsd(reference, mobile):
    reference = np.asarray(reference, dtype=float)
    mobile = np.asarray(mobile, dtype=float)
    if reference.shape != mobile.shape:
        return float("inf")
    reference_centered = reference - np.mean(reference, axis=0)
    mobile_centered = mobile - np.mean(mobile, axis=0)
    covariance = mobile_centered.T @ reference_centered
    left_vectors, _, right_vectors = np.linalg.svd(covariance)
    correction = np.eye(3)
    if np.linalg.det(left_vectors @ right_vectors) < 0:
        correction[-1, -1] = -1.0
    rotation = left_vectors @ correction @ right_vectors
    aligned = mobile_centered @ rotation
    difference = aligned - reference_centered
    return float(np.sqrt(np.mean(np.sum(difference * difference, axis=1))))


def _read_xyz(path):
    return Molecule.from_xyz(str(path))


def _geometry_rmsd(reference, candidate):
    if _heavy_symbols(reference) != _heavy_symbols(candidate):
        return float("inf")
    reference_coordinates = _heavy_coordinates(reference)
    candidate_coordinates = _heavy_coordinates(candidate)
    if reference_coordinates is None or candidate_coordinates is None:
        return float("inf")
    return _kabsch_rmsd(reference_coordinates, candidate_coordinates)


def _read_summary_rows(summary_path):
    if not summary_path.is_file():
        return []
    with summary_path.open(encoding="utf-8", newline="") as fp:
        return list(csv.DictReader(fp))


def _float_or_none(value):
    if value in {None, ""}:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _int_or_none(value):
    if value in {None, ""}:
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _candidate_entries(summary_rows, run_directory, field):
    entries = []
    for row_index, row in enumerate(summary_rows):
        raw_path = row.get(field)
        if not raw_path:
            continue
        path = Path(raw_path)
        if not path.is_absolute():
            path = run_directory / path
        if not path.is_file():
            continue
        rank = _int_or_none(row.get("rank"))
        if rank is None:
            rank = row_index
        entries.append(
            {
                "path": path,
                "rank": rank,
                "rdkit_energy": _float_or_none(row.get("rdkit_energy")),
                "backend_energy": _float_or_none(row.get("backend_energy")),
                "summary_row": row,
            }
        )
    return entries


def _selected_entries(summary_rows, selected_directory):
    selected_paths = sorted(selected_directory.glob("*.xyz")) if selected_directory.is_dir() else []
    by_path = {}
    for row in summary_rows:
        selected_path = row.get("selected_path")
        if selected_path:
            by_path[Path(selected_path).name] = row
    entries = []
    for rank, path in enumerate(selected_paths):
        row = by_path.get(path.name, {})
        entries.append(
            {
                "path": path,
                "rank": rank,
                "rdkit_energy": _float_or_none(row.get("rdkit_energy")),
                "backend_energy": _float_or_none(row.get("backend_energy")),
                "summary_row": row,
            }
        )
    return entries


def _nearest(reference, entries):
    best = None
    for entry in entries:
        try:
            candidate = _read_xyz(entry["path"])
            rmsd = _geometry_rmsd(reference, candidate)
        except Exception as exc:
            rmsd = float("inf")
            entry = dict(entry)
            entry["read_error"] = str(exc)
        key = (rmsd, int(entry.get("rank", 0)))
        if best is None or key < best["key"]:
            best = {
                "key": key,
                "path": str(entry["path"]),
                "rank": entry.get("rank"),
                "rmsd": rmsd,
                "rdkit_energy": entry.get("rdkit_energy"),
                "backend_energy": entry.get("backend_energy"),
            }
    if best is None:
        return {
            "path": None,
            "rank": None,
            "rmsd": float("inf"),
            "rdkit_energy": None,
            "backend_energy": None,
        }
    best.pop("key", None)
    return best


def _copy_if_available(source, target):
    if source:
        source_path = Path(source)
        if source_path.is_file():
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source_path, target)


def _load_state(run_directory):
    state_path = run_directory / "state.json"
    if not state_path.is_file():
        return {}
    with state_path.open(encoding="utf-8") as fp:
        return json.load(fp)


def _lowest(values):
    numeric_values = [value for value in values if value is not None]
    return min(numeric_values) if numeric_values else None


def _classify(
    *,
    reference_path,
    reference_error,
    generated_hit,
    selected_hit,
    backend_requested,
    generated_entries,
    nearest_generated,
    rms_hit_threshold,
    energy_window,
):
    if reference_path is None:
        return "input_chemistry_issue", "case has no reference_xyz"
    if reference_error:
        return "input_chemistry_issue", reference_error
    if selected_hit:
        return "found", "reference-like conformer is present in final selected outputs"
    if generated_hit:
        if energy_window is not None and nearest_generated.get("rdkit_energy") is not None:
            lowest = _lowest(entry.get("rdkit_energy") for entry in generated_entries)
            if lowest is not None and (nearest_generated["rdkit_energy"] - lowest) * 627.509474 > float(energy_window):
                return "wrong_ranking", "reference-like conformer is outside the RDKit energy window"
        if backend_requested:
            return (
                "generated_lost_backend_or_dedup",
                "reference-like conformer was generated but absent from final backend/deduplicated selection",
            )
        return "generated_lost_selection", "reference-like conformer was generated but not selected"
    if math.isfinite(nearest_generated.get("rmsd", float("inf"))):
        return "never_generated", "no generated conformer is within the RMSD hit threshold"
    return "uncertain", "no generated conformer geometries were available for comparison"


def diagnose_conformer_case(
    case,
    conformer_result,
    *,
    rms_hit_threshold=0.5,
    energy_window=None,
):
    """Diagnose why one conformer benchmark case did or did not find a reference."""
    run_directory = Path(conformer_result.run_directory)
    state = _load_state(run_directory)
    summary_rows = _read_summary_rows(run_directory / "summary.csv")
    reference_path = Path(case.reference_xyz) if case.reference_xyz else None
    backend_requested = bool(state.get("backend_requested"))
    if not state and getattr(conformer_result, "metadata", None):
        backend_requested = bool(conformer_result.metadata.get("backend_requested"))

    reference = None
    reference_error = None
    if reference_path is None:
        reference_error = "case has no reference_xyz"
    elif not reference_path.is_file():
        reference_error = f"reference geometry is missing: {reference_path}"
    else:
        try:
            reference = _read_xyz(reference_path)
        except Exception as exc:
            reference_error = f"could not read reference geometry: {exc}"

    generated_entries = _candidate_entries(summary_rows, run_directory, "rdkit_path")
    selected_entries = _selected_entries(summary_rows, run_directory / "selected")
    if reference is not None:
        nearest_generated = _nearest(reference, generated_entries)
        nearest_selected = _nearest(reference, selected_entries)
    else:
        nearest_generated = _nearest_empty()
        nearest_selected = _nearest_empty()

    generated_hit = nearest_generated["rmsd"] <= float(rms_hit_threshold)
    selected_hit = nearest_selected["rmsd"] <= float(rms_hit_threshold)
    failure_class, reason = _classify(
        reference_path=reference_path,
        reference_error=reference_error,
        generated_hit=generated_hit,
        selected_hit=selected_hit,
        backend_requested=backend_requested,
        generated_entries=generated_entries,
        nearest_generated=nearest_generated,
        rms_hit_threshold=rms_hit_threshold,
        energy_window=energy_window,
    )

    lowest_backend = _lowest(_float_or_none(row.get("backend_energy")) for row in summary_rows)
    lowest_rdkit = _lowest(_float_or_none(row.get("rdkit_energy")) for row in summary_rows)
    return {
        "case_id": case.id,
        "input": case.input,
        "input_format": case.input_format,
        "status": state.get("status", getattr(conformer_result, "status", "unknown")),
        "failure_class": failure_class,
        "failure_class_description": FAILURE_CLASSES.get(failure_class),
        "reason": reason,
        "reference_found_generated": bool(generated_hit),
        "reference_found_selected": bool(selected_hit),
        "best_generated_rmsd": nearest_generated["rmsd"],
        "best_selected_rmsd": nearest_selected["rmsd"],
        "best_generated_rank": nearest_generated["rank"],
        "best_selected_rank": nearest_selected["rank"],
        "lowest_rdkit_energy": lowest_rdkit,
        "lowest_backend_energy": lowest_backend,
        "selected_count": len(selected_entries),
        "generated_count": len(generated_entries),
        "backend_requested": backend_requested,
        "notes": case.notes,
        "reference_xyz": str(reference_path) if reference_path else None,
        "nearest_generated": nearest_generated,
        "nearest_selected": nearest_selected,
    }


def _nearest_empty():
    return {
        "path": None,
        "rank": None,
        "rmsd": float("inf"),
        "rdkit_energy": None,
        "backend_energy": None,
    }


def _summary_row(diagnosis):
    return {column: diagnosis.get(column) for column in SUMMARY_COLUMNS}


def _case_output_directory(output_directory, case):
    safe_id = "".join(char if char.isalnum() or char in {"-", "_"} else "_" for char in case.id)
    return Path(output_directory) / "cases" / safe_id


def _default_backend_qc_params(software):
    return {
        "basis": defualt_parameters.values["basis"],
        "method": defualt_parameters.values["method"],
        "software": software,
        "geometry_optimizer": defualt_parameters.values["geometry_optimizer"],
        "opt_target": defualt_parameters.values["opt_target"],
        "opt_cycles": defualt_parameters.values["opt_cycles"],
        "opt_threshold": defualt_parameters.values["opt_threshold"],
        "scf_cycles": defualt_parameters.values["scf_cycles"],
        "scf_threshold": defualt_parameters.values["scf_threshold"],
        "nprocs": defualt_parameters.values["nprocs"],
        "gamma": None,
        "custom_keywords": None,
        "custom_keyword": None,
        "model": defualt_parameters.values["model"],
    }


def run_conformer_benchmark(
    benchmark,
    *,
    output="conformer_benchmark_results",
    num_conformers=100,
    num_seeds=3,
    top_n=10,
    backend_top_n=None,
    software=None,
    rms_hit_threshold=0.5,
    energy_window=None,
    conformer_options=None,
    qc_params=None,
):
    """Run a JSON conformer benchmark and write summary/diagnosis artifacts."""
    spec = load_conformer_benchmark(benchmark) if not isinstance(benchmark, ConformerBenchmarkSpec) else benchmark
    benchmark_root = Path(spec.path).resolve().parent if spec.path else Path.cwd()
    output_directory = Path(output)
    output_directory.mkdir(parents=True, exist_ok=True)
    rows = []
    diagnoses = []
    conformer_options = dict(conformer_options or {})
    backend_params = dict(qc_params or {}) if qc_params is not None else None
    if backend_params is None and software:
        backend_params = _default_backend_qc_params(software)

    for case in spec.cases:
        case_directory = _case_output_directory(output_directory, case)
        run_root = case_directory / "run"
        case_directory.mkdir(parents=True, exist_ok=True)
        reference_path = _resolve_case_path(case, benchmark_root)
        case_for_diagnosis = case
        if reference_path is not None:
            case_for_diagnosis = ConformerBenchmarkCase(
                **{**asdict(case), "reference_xyz": str(reference_path)}
            )
            _copy_if_available(reference_path, case_directory / "reference.xyz")

        try:
            result = conformer_workflow.conformer_search(
                case.input,
                input_format=case.input_format,
                num_conformers=num_conformers,
                top_n=top_n,
                backend_top_n=backend_top_n,
                num_seeds=num_seeds,
                charge=case.charge,
                multiplicity=case.multiplicity,
                qc_params=backend_params,
                root_directory=str(run_root),
                **conformer_options,
            )
            diagnosis = diagnose_conformer_case(
                case_for_diagnosis,
                result,
                rms_hit_threshold=rms_hit_threshold,
                energy_window=energy_window,
            )
        except Exception as exc:
            result = ConformerResult(
                workflow="conformer",
                status="failed",
                run_directory=str(run_root / "conformers"),
                state_path=str(run_root / "conformers" / "state.json"),
                selected_paths=(),
                metadata={"backend_requested": bool(backend_params)},
            )
            diagnosis = diagnose_conformer_case(
                case_for_diagnosis,
                result,
                rms_hit_threshold=rms_hit_threshold,
                energy_window=energy_window,
            )
            diagnosis["status"] = "failed"
            diagnosis["failure_class"] = (
                "input_chemistry_issue"
                if diagnosis["failure_class"] == "input_chemistry_issue"
                else "uncertain"
            )
            diagnosis["reason"] = f"conformer workflow failed: {exc}"

        _copy_if_available(diagnosis["nearest_generated"].get("path"), case_directory / "nearest_generated.xyz")
        _copy_if_available(diagnosis["nearest_selected"].get("path"), case_directory / "nearest_selected.xyz")
        _write_json(case_directory / "diagnosis.json", diagnosis)
        diagnoses.append(diagnosis)
        rows.append(_summary_row(diagnosis))

    _write_csv(output_directory / "benchmark_summary.csv", rows)
    payload = {
        "name": spec.name,
        "method": spec.method,
        "failure_classes": FAILURE_CLASSES,
        "cases": diagnoses,
    }
    _write_json(output_directory / "benchmark_summary.json", payload)
    return payload
