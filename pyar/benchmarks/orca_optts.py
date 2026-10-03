"""Run and validate ORCA OptTS reference calculations for TS benchmarks."""

from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
import subprocess
import time
from pathlib import Path

from pyar.benchmarks.ts_optimizer import (
    TSOptimizerBenchmarkError,
    VALIDATION_SETTINGS,
    _hessian_metrics,
    _runtime_payload,
    _stage_arguments,
    _stage_counts,
    _write_json,
    load_ts_optimizer_benchmark,
)
from pyar.neb import _aligned_rmsd, read_xyz, run_neb


def _sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def build_orca_optts_input(
    symbols, coordinates, *, charge, multiplicity, nprocs=1, max_cycles=200,
    ts_fmax_ev_per_angstrom=None,
):
    """Format ORCA OptTS with the GFN2-xTB external-method interface."""
    lines = [
        "! GFN2-xTB OptTS",
        f"%pal nprocs {int(nprocs)} end",
        "%geom",
        "  CoordSys redundant",
        f"  MaxIter {int(max_cycles)}",
        "  InHess XTB2",
        "  Update Bofill",
    ]
    if ts_fmax_ev_per_angstrom is not None:
        from ase.units import Bohr, Hartree

        if ts_fmax_ev_per_angstrom <= 0:
            raise ValueError("ts_fmax_ev_per_angstrom must be positive")
        max_gradient_hartree_bohr = ts_fmax_ev_per_angstrom * Bohr / Hartree
        rms_gradient_hartree_bohr = max_gradient_hartree_bohr / 1.5
        lines.extend((
            f"  TolMaxG {max_gradient_hartree_bohr:.10g}",
            f"  TolRMSG {rms_gradient_hartree_bohr:.10g}",
        ))
    lines.extend(("end", f"*xyz {int(charge)} {int(multiplicity)}"))
    lines.extend(
        f"{symbol:>3} {x: .10f} {y: .10f} {z: .10f}"
        for symbol, (x, y, z) in zip(symbols, coordinates)
    )
    lines.append("*")
    return "\n".join(lines) + "\n"


def _classify_orca_failure(output_text, returncode, converged):
    lowered = output_text.lower()
    if "scf not converged" in lowered or "scf convergence failure" in lowered:
        return "scf_failure"
    if any(token in lowered for token in (
        "maximum number of geometry optimization cycles reached",
        "maximum number of optimization cycles reached",
        "the optimization did not converge",
        "maximum number of cycles reached",
    )):
        return "max_iterations"
    if returncode != 0:
        if "hessian" in lowered and any(
            token in lowered for token in ("failed", "error", "could not", "not available")
        ):
            return "hessian_failure"
        return "orca_failure"
    if not converged:
        return "optimizer_not_converged"
    return None


def parse_orca_optts_output(output_text, returncode=0):
    """Extract native optimizer and explicitly sourceable Hessian diagnostics."""
    normally_terminated = "****ORCA TERMINATED NORMALLY****" in output_text
    converged = "***        THE OPTIMIZATION HAS CONVERGED     ***" in output_text
    cycles = len(re.findall(r"GEOMETRY OPTIMIZATION CYCLE\s+\d+", output_text))
    # ORCA 6.1's GFN2-xTB interface does not print the SCF-gradient banner.
    # It does print one timing record for each optimizer energy/gradient cycle.
    gradient_events = len(re.findall(r"Time for energy\+gradient\s*:", output_text, re.I))
    hessian_times = [
        float(value)
        for value in re.findall(
            r"Hessian update/contruction\s*:\s*([0-9]*\.?[0-9]+)\s*s",
            output_text,
            flags=re.IGNORECASE,
        )
    ]
    initial_hessian_events = len(re.findall(
        r"Evaluating the initial hessian", output_text, flags=re.IGNORECASE,
    ))
    failure = _classify_orca_failure(output_text, returncode, converged)
    if not normally_terminated and failure is None:
        failure = "orca_failure"
    return {
        "normally_terminated": normally_terminated,
        "optimizer_converged": converged and normally_terminated and returncode == 0,
        "failure_class": failure,
        "optimizer_steps": cycles,
        "orca_gradient_events": gradient_events,
        "initial_hessian_source": "ORCA GFN2-xTB model Hessian (InHess XTB2)",
        "initial_hessian_events": initial_hessian_events,
        "hessian_update_policy": "Bofill update each geometry cycle; no scheduled full Hessian refresh",
        "hessian_update_construction_wall_seconds": sum(hessian_times),
        "hessian_update_timing_count": len(hessian_times),
        "convergence_criteria": "ORCA default OptTS criteria; exact thresholds retained in .out",
    }


def _matching_line(text, pattern):
    return next((line.strip() for line in text.splitlines() if re.search(pattern, line, re.I)), None)


def _xtb_version(executable):
    process = subprocess.run(
        [str(Path(executable).resolve()), "--version"], capture_output=True,
        text=True, check=False, timeout=15,
    )
    text = process.stdout or process.stderr or ""
    line = _matching_line(text, r"xtb version\s+[0-9.]+")
    match = re.search(r"xtb version\s+(\d+)\.(\d+)\.(\d+)", line or "", re.I)
    if process.returncode != 0 or not match:
        raise TSOptimizerBenchmarkError("could not determine the selected xTB executable version")
    version = tuple(int(part) for part in match.groups())
    if version < (6, 7, 1):
        raise TSOptimizerBenchmarkError(
            f"ORCA's GFN2-xTB interface requires xTB 6.7.1 or later; found {match.group(1)}."
            f"{match.group(2)}.{match.group(3)}"
        )
    return line


def run_orca_optts_case(
    benchmark, *, case_id, output, orca_executable=None, xtb_executable=None,
    match_ts_fmax=False,
):
    """Run native ORCA OptTS and validate the resulting saddle with PyAR stages."""
    spec = benchmark if hasattr(benchmark, "cases") else load_ts_optimizer_benchmark(benchmark)
    if spec.qc_model.get("xtb_model") != "gfn2":
        raise TSOptimizerBenchmarkError(
            "ORCA OptTS reference currently requires qc_model.xtb_model='gfn2'"
        )
    case = next((item for item in spec.cases if item.id == case_id), None)
    if case is None:
        raise TSOptimizerBenchmarkError(f"unknown case id: {case_id}")
    orca_executable = orca_executable or shutil.which("orca")
    xtb_executable = xtb_executable or shutil.which("xtb")
    if not orca_executable or not Path(orca_executable).is_file():
        raise TSOptimizerBenchmarkError("ORCA executable not found; pass --orca-executable")
    if not xtb_executable or not Path(xtb_executable).is_file():
        raise TSOptimizerBenchmarkError("xTB executable not found; pass --xtb-executable")
    provider_xtb = shutil.which("xtb")
    if not provider_xtb or Path(provider_xtb).resolve() != Path(xtb_executable).resolve():
        raise TSOptimizerBenchmarkError(
            "--xtb-executable must match xtb on PATH, which PyAR uses for validation; "
            "put the selected executable on PATH before running this benchmark"
        )
    xtb_version = _xtb_version(xtb_executable)

    run_directory = Path(output).resolve() / "cases" / case.id / "orca_optts"
    result_path = run_directory / "benchmark_result.json"
    if result_path.exists():
        raise TSOptimizerBenchmarkError(f"ORCA result already exists: {result_path}")
    if run_directory.exists() and any(run_directory.iterdir()):
        raise TSOptimizerBenchmarkError(
            f"ORCA case directory contains artifacts from an incomplete run: {run_directory}; "
            "use a fresh output directory rather than reusing potentially stale optimizer files"
        )
    run_directory.mkdir(parents=True, exist_ok=True)
    copied = {}
    for label in ("ts_guess", "reference_ts", "reactant", "product"):
        source = getattr(case, label)
        if source:
            destination = run_directory / f"input_{label}.xyz"
            shutil.copy2(source, destination)
            copied[label] = destination
    symbols, coordinates, _ = read_xyz(copied["ts_guess"])
    reference_symbols, reference_coordinates, _ = read_xyz(copied["reference_ts"])
    if symbols != reference_symbols:
        raise TSOptimizerBenchmarkError(f"case {case.id} reference atom order differs from input")

    input_path = run_directory / "orca_optts.inp"
    output_path = run_directory / "orca_optts.out"
    input_path.write_text(build_orca_optts_input(
        symbols, coordinates, charge=case.charge, multiplicity=case.multiplicity,
        nprocs=spec.qc_model.get("nprocs", 1), max_cycles=spec.settings["ts_max_cycles"],
        ts_fmax_ev_per_angstrom=(spec.settings["ts_fmax"] if match_ts_fmax else None),
    ), encoding="utf-8")
    environment = os.environ.copy()
    environment["XTBEXE"] = str(Path(xtb_executable).resolve())
    started = time.perf_counter()
    try:
        process = subprocess.run(
            [str(Path(orca_executable).resolve()), str(input_path)],
            cwd=run_directory, env=environment, capture_output=True, text=True,
            check=False,
        )
        output_text = (process.stdout or "") + "\n" + (process.stderr or "")
        returncode = process.returncode
    except OSError as exc:
        output_text = f"{type(exc).__name__}: {exc}"
        returncode = 127
    elapsed = time.perf_counter() - started
    output_path.write_text(output_text, encoding="utf-8")
    parsed = parse_orca_optts_output(output_text, returncode)
    optimized_path = run_directory / "orca_optts.xyz"
    geometry_available = optimized_path.is_file()
    geometry_error = None
    if geometry_available:
        try:
            final_symbols, final_coordinates, final_comment = read_xyz(optimized_path)
            if final_symbols != symbols:
                raise ValueError("optimized geometry atom identities/order differ from input")
        except (OSError, ValueError) as exc:
            geometry_available = False
            geometry_error = str(exc)
    if parsed["optimizer_converged"] and not geometry_available:
        parsed["failure_class"] = "optimized_geometry_invalid" if geometry_error else "optimized_geometry_missing"
    optimizer_converged = parsed["optimizer_converged"] and geometry_available

    ts_summary = {
        **parsed,
        "geometry_error": geometry_error,
        "convergence_criteria": (
            "matched gradient thresholds from benchmark ts_fmax plus ORCA default "
            "energy/displacement criteria"
            if match_ts_fmax else parsed["convergence_criteria"]
        ),
        "optimizer_converged": optimizer_converged,
        "optimizer_wall_seconds": elapsed,
        "input_sha256": _sha256(copied["ts_guess"]),
        "optimized_geometry": optimized_path.name if geometry_available else None,
        "orca_version": _matching_line(output_text, r"Program Version\s+[0-9.]+"),
        "orca_executable": str(Path(orca_executable).resolve()),
        "xtb_executable": str(Path(xtb_executable).resolve()),
        "xtb_version": xtb_version,
    }

    runtime = _runtime_payload(spec.qc_model)
    outcome = parsed["failure_class"] or "optimizer_not_converged"
    failure_stage = "ts"
    frequency_result = endpoint_result = None
    base = _stage_arguments(spec, case, run_directory)
    if optimizer_converged:
        outcome = "converged_not_stationary"
        failure_stage = "frequency"
        try:
            frequency_result = run_neb(
                ts_geometry=optimized_path, stage="frequency",
                **base,
            )
            if frequency_result.get("stationary") is True:
                if frequency_result.get("first_order_saddle_confirmed") is True:
                    outcome = "validated_first_order_saddle"
                    if case.reactant and case.product:
                        failure_stage = "relax"
                        run_neb(
                            start=copied["reactant"], end=copied["product"], stage="relax",
                            **base,
                        )
                        failure_stage = "irc"
                        run_neb(stage="irc", **base)
                        failure_stage = "endpoints"
                        endpoint_result = run_neb(
                            start=copied["reactant"], end=copied["product"], stage="endpoints",
                            **base,
                        )
                        outcome = (
                            "reaction_connected_success"
                            if endpoint_result.get("reactant_product_connection_confirmed") is True
                            else "first_order_saddle_wrong_connection"
                        )
                else:
                    outcome = "stationary_not_first_order_saddle"
        except Exception as exc:
            outcome = {
                "frequency": "frequency_exception",
                "relax": "endpoint_relaxation_exception",
                "irc": "irc_exception",
                "endpoints": "endpoint_validation_exception",
            }.get(failure_stage, "orca_validation_exception")
            ts_summary["validation_error"] = f"{type(exc).__name__}: {exc}"

    reference_rmsd = energy_difference = None
    if geometry_available:
        if final_symbols == reference_symbols:
            reference_rmsd = _aligned_rmsd(reference_coordinates, final_coordinates)
        energy_match = re.search(r"\bE\s+([-+0-9.Ee]+)", final_comment)
        if energy_match and case.reference_energy_hartree is not None:
            energy_difference = float(energy_match.group(1)) - case.reference_energy_hartree
    validation_evaluations, validation_wall = _stage_counts(
        run_directory, ("frequency", "relax", "irc", "endpoints"),
    )
    validation_hessian = _hessian_metrics(run_directory, ("frequency", "endpoints"))
    successful_or_validated = outcome in {
        "reaction_connected_success", "first_order_saddle_wrong_connection",
        "validated_first_order_saddle", "stationary_not_first_order_saddle",
        "converged_not_stationary",
    }
    result = {
        "schema_version": 1,
        "benchmark_name": spec.name,
        "case_id": case.id,
        "difficulty": case.difficulty,
        "optimizer": "orca_optts",
        "optimizer_class": "transition_state",
        "status": "complete" if successful_or_validated else "failed",
        "outcome": outcome,
        "failure_stage": None if outcome in {
            "reaction_connected_success", "first_order_saddle_wrong_connection",
            "validated_first_order_saddle",
        } else failure_stage,
        "error": ts_summary.get("validation_error", parsed["failure_class"]),
        "source": spec.source,
        "qc_model": {
            "software": "orca",
            "method": "GFN2-xTB via ORCA xtb interface",
            "xtb_model": "gfn2",
            "xtb_executable": ts_summary["xtb_executable"],
            "xtb_version": ts_summary["xtb_version"],
        },
        "case_qc_settings": {"charge": case.charge, "multiplicity": case.multiplicity},
        "settings": {
            **{key: spec.settings[key] for key in VALIDATION_SETTINGS},
            "orca_coordinate_system": "redundant_internal",
            "orca_optts_max_cycles": spec.settings["ts_max_cycles"],
            "orca_convergence_criteria": (
                "gradient thresholds matched to benchmark ts_fmax; ORCA energy and "
                "displacement thresholds remain at ORCA defaults"
                if match_ts_fmax else
                "ORCA default OptTS criteria; thresholds recorded in .out"
            ),
            "reference_pyAR_ts_fmax_ev_per_angstrom": spec.settings["ts_fmax"],
            "reference_pyAR_ts_fmax_applied_to_orca": bool(match_ts_fmax),
            "imaginary_frequency_threshold": spec.settings["imaginary_frequency_threshold"],
            "irc_max_cycles": spec.settings["irc_max_cycles"],
            "endpoint_max_cycles": spec.settings["endpoint_max_cycles"],
        },
        "runtime": {
            **runtime,
            "orca_version": ts_summary["orca_version"],
            "orca_executable": ts_summary["orca_executable"],
        },
        "input_sha256": _sha256(copied["ts_guess"]),
        "input_hashes": {key: _sha256(path) for key, path in copied.items()},
        "reference_rmsd_angstrom": reference_rmsd,
        "reference_energy_difference_hartree": energy_difference,
        "ts_optimizer_steps": ts_summary["optimizer_steps"],
        "ts_orca_gradient_events": ts_summary["orca_gradient_events"],
        "ts_backend_energy_gradient_evaluations": None,
        "ts_energy_gradient_counter": {
            "value": ts_summary["orca_gradient_events"],
            "source": "count of ORCA optimizer energy+gradient timing records",
            "equivalent_to_provider_calls": False,
        },
        "ts_wall_seconds": elapsed,
        "ts_hessian": {key: ts_summary[key] for key in (
            "initial_hessian_source", "initial_hessian_events", "hessian_update_policy",
            "hessian_update_construction_wall_seconds", "hessian_update_timing_count",
        )},
        "validation_backend_energy_gradient_evaluations": validation_evaluations,
        "validation_wall_seconds": validation_wall,
        "validation_hessian_evaluations": validation_hessian["evaluations"],
        "validation_hessian_wall_seconds": validation_hessian["wall_seconds"],
        "hessian_sources": validation_hessian["sources"],
        "ts_result": ts_summary,
        "frequency_result": frequency_result,
        "endpoint_result": endpoint_result,
        "run_directory": str(run_directory),
    }
    _write_json(result_path, result)
    return result
