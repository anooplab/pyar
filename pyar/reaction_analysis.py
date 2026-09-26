"""AFIR reaction-path analysis for recorded geomeTRIC traces.

The analysis layer reads the JSONL/XYZ trace emitted by the reaction
workflow, writes a CSV summary, exports candidate geometries for the most
interesting path points, and can generate standalone plots from the same
recorded trace.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from pyar.data.units import atomic_unit2kilo_calories
from pyar.reaction_trace import load_trace_records


def _sorted_records(records):
    """Return records ordered by their recorded step index."""
    return sorted(records, key=lambda record: int(record.get("step_index", 0)))


def _validate_sorted_trace(records):
    """Fail clearly when the trace contains duplicate or malformed step indices."""
    previous_step = None
    seen_steps = set()
    for record in records:
        step_index = int(record["step_index"])
        if step_index in seen_steps:
            raise ValueError(f"Duplicate trace step_index {step_index} detected")
        if previous_step is not None and step_index < previous_step:
            raise ValueError("Trace records are not in ascending step_index order")
        seen_steps.add(step_index)
        previous_step = step_index


def _relative_kcalmol(value, baseline):
    """Return a Hartree difference converted to kcal/mol."""
    return atomic_unit2kilo_calories(float(value) - float(baseline))


def _record_index_with_max(records, key):
    """Return the index of the record with the largest value for ``key``."""
    values = [float(record[key]) for record in records]
    return int(np.argmax(values))


def _bond_sets(records):
    """Return bond sets normalized as tuples for every record."""
    bond_sets = []
    for record in records:
        current_bonds = {
            tuple(pair)
            for pair in record.get("current_bonds", [])
        }
        bond_sets.append(current_bonds)
    return bond_sets


def _persistent_transition_index(records):
    """Return the first persistent topology change after the reactant frame."""
    if len(records) < 2:
        return 0

    bond_sets = _bond_sets(records)
    for index in range(1, len(records)):
        if bond_sets[index] == bond_sets[index - 1]:
            continue
        if index == len(records) - 1:
            return index
        if bond_sets[index] == bond_sets[index + 1]:
            return index
        if index + 2 < len(records) and bond_sets[index] == bond_sets[index + 1] == bond_sets[index + 2]:
            return index

    for index in range(1, len(records)):
        if bond_sets[index] != bond_sets[index - 1]:
            return index
    return 0


def _pre_product_index(records, eligible_indices=None):
    """Return the highest electronic-energy eligible frame before the event."""
    transition_index = _persistent_transition_index(records)
    eligible = range(len(records)) if eligible_indices is None else eligible_indices
    return max(
        (index for index in eligible if index < transition_index),
        key=lambda index: float(records[index]["backend_energy_hartree"]), default=None,
    )


def _energy_outlier_scores(records):
    """Robust scores from unbiased backend energies, including zero-MAD data."""
    energies = np.asarray(
        [float(record["backend_energy_hartree"]) for record in records], dtype=float
    )
    median = float(np.median(energies))
    deviations = np.abs(energies - median)
    mad = float(np.median(deviations))
    if mad > 0.0:
        return 0.67448975 * deviations / mad
    # A repeated/quantized majority has MAD=0. Treat values equal within
    # floating-point noise as inliers and every distinct value as an outlier.
    tolerance = max(1.0e-12, 8.0 * np.finfo(float).eps * max(1.0, abs(median)))
    return np.where(deviations <= tolerance, 0.0, np.inf)


def _physical_max_force(record):
    """Recover physical force from legacy vectors; never substitute biased force."""
    value = record.get("backend_max_force")
    if value is not None:
        return float(value)
    forces = record.get("backend_forces_hartree_per_bohr")
    if forces is None:
        return None
    forces = np.asarray(forces, dtype=float)
    if forces.ndim != 2 or forces.shape[1] != 3 or not forces.size or not np.all(np.isfinite(forces)):
        return None
    return float(np.max(np.linalg.norm(forces, axis=1)))


def _candidate_record_indices(records, max_force=None, energy_outlier_z=None):
    """Return trace-record indices eligible for candidate geometry selection."""
    eligible = list(range(len(records)))
    if max_force is not None:
        max_force = float(max_force)
        if not np.isfinite(max_force) or max_force <= 0.0:
            raise ValueError("max_force must be a finite positive value in Hartree/bohr")
        eligible = [
            index for index in eligible
            if _physical_max_force(records[index]) is not None
            and float(_physical_max_force(records[index])) <= max_force
        ]
    if energy_outlier_z is not None:
        energy_outlier_z = float(energy_outlier_z)
        if not np.isfinite(energy_outlier_z) or energy_outlier_z <= 0.0:
            raise ValueError("energy_outlier_z must be a finite positive value")
        robust_z = _energy_outlier_scores(records)
        eligible = [index for index in eligible if robust_z[index] <= energy_outlier_z]
    if not eligible:
        raise ValueError("Candidate filters excluded every reaction-trace geometry")
    return eligible


def _candidate_exclusion_reasons(records, max_force=None, energy_outlier_z=None):
    """Describe why optional force/outlier filters reject individual frames."""
    reasons = {}
    if max_force is not None:
        for index, record in enumerate(records):
            force = _physical_max_force(record)
            if force is None:
                reasons.setdefault(index, []).append("physical_force_unavailable")
            elif float(force) > float(max_force):
                reasons.setdefault(index, []).append("max_force")
    if energy_outlier_z is not None:
        robust_z = _energy_outlier_scores(records)
        for index, score in enumerate(robust_z):
            if score > float(energy_outlier_z):
                reasons.setdefault(index, []).append("backend_energy_outlier")
    return reasons


def _first_persistent_transition_index(records, eligible_indices=None):
    """Select the first persistent topology event, optionally requiring eligibility."""
    eligible = None if eligible_indices is None else set(eligible_indices)
    transition_index = _persistent_transition_index(records)
    if transition_index == 0:
        return None
    if eligible is None or transition_index in eligible:
        return transition_index
    bond_sets = _bond_sets(records)
    transition_bonds = bond_sets[transition_index]
    for index in range(transition_index + 1, len(records)):
        if bond_sets[index] != transition_bonds:
            break
        if index in eligible:
            return index
    return None


def _step_xyz_path(record):
    """Return the relative path for a trace snapshot XYZ file."""
    step_index = int(record["step_index"])
    return str(Path("reaction_trace") / "steps" / f"step_{step_index:06d}.xyz")


def _resolve_trace_directory(path):
    """Return the job directory and trace directory for a trace path."""
    path = Path(path)
    if (path / "trace.jsonl").exists() or (path / "steps").exists():
        return path.parent, path

    trace_directory = path / "reaction_trace"
    if (trace_directory / "trace.jsonl").exists() or (trace_directory / "steps").exists():
        return path, trace_directory

    return path, trace_directory


def _xyz_comment(record, label, reference_record=None, energy_key=None):
    """Build a compact but informative XYZ comment line for candidate files."""
    reference_record = reference_record or record
    backend_rel = _relative_kcalmol(
        record["backend_energy_hartree"],
        reference_record["backend_energy_hartree"],
    )
    bias_rel = _relative_kcalmol(
        record["bias_energy_hartree"],
        reference_record["bias_energy_hartree"],
    )
    total_rel = _relative_kcalmol(
        record["total_energy_hartree"],
        reference_record["total_energy_hartree"],
    )
    step_index = int(record["step_index"])
    current_bonds = len(record.get("current_bonds", []))
    bond_changes = int(record.get("bond_change_count", 0))
    min_distance = record.get("min_interfragment_distance_angstrom")
    parts = [
        f"candidate={label}",
        f"step_index={step_index}",
        f"source_energy_key={energy_key or 'total_energy_hartree'}",
        f"backend={float(record['backend_energy_hartree']):.12f} Ha",
        f"bias={float(record['bias_energy_hartree']):.12f} Ha",
        f"total={float(record['total_energy_hartree']):.12f} Ha",
        f"backend_rel={backend_rel:.6f} kcal/mol",
        f"bias_rel={bias_rel:.6f} kcal/mol",
        f"total_rel={total_rel:.6f} kcal/mol",
        f"bond_changes={bond_changes}",
        f"bonds={current_bonds}",
    ]
    if min_distance is not None:
        parts.append(f"min_interfragment_distance={float(min_distance):.6f} angstrom")
    coordinate = record.get("collective_coordinate_bohr")
    if coordinate is not None:
        parts.append(f"collective_coordinate={float(coordinate):.6f} bohr")
    diagnostics = record.get("contact_diagnostics") or {}
    if diagnostics:
        parts.append(f"contacts={int(diagnostics.get('contact_count', 0))}")
    return " ".join(parts)


def _write_xyz_record(path, record, label, reference_record=None, energy_key="total_energy_hartree"):
    """Write one trace record as an XYZ file with a rich comment line."""
    coordinates = np.asarray(record["coordinates_angstrom"], dtype=float)
    symbols = list(record["symbols"])
    comment = _xyz_comment(record, label, reference_record=reference_record, energy_key=energy_key)
    with Path(path).open("w", encoding="utf-8") as fp:
        fp.write(f"{len(symbols)}\n")
        fp.write(f"{comment}\n")
        for symbol, (x, y, z) in zip(symbols, coordinates):
            fp.write(f"{symbol:<2}{x:12.5f}{y:12.5f}{z:12.5f}\n")


def _candidate_metadata(
    job_directory,
    trace_records,
    highest_backend_index,
    pre_product_index,
    topology_change_index,
    highest_total_index,
    transition_index,
    baseline_record,
):
    """Return metadata for the candidate TS geometries."""
    candidate_directory = Path(job_directory) / "candidate_ts"
    candidate_files = {}
    for label, index in (
        ("highest_backend_energy", highest_backend_index),
        ("pre_product_geometry", pre_product_index),
        ("first_topology_change", topology_change_index),
        ("highest_total_energy", highest_total_index),
    ):
        candidate_files[label] = {
            "available": index is not None,
            "xyz_file": f"{label}.xyz" if index is not None else None,
            "source_step_index": int(trace_records[index]["step_index"]) if index is not None else None,
            "source_trace_xyz": _step_xyz_path(trace_records[index]) if index is not None else None,
        }
        if index is None:
            candidate_files[label]["reason"] = (
                "no_topology_change" if transition_index == 0 else "no_eligible_event_geometry"
            )
    return {
        "path_summary_csv": str(Path(job_directory) / "path_summary.csv"),
        "trace_file": str(Path(job_directory) / "reaction_trace" / "trace.jsonl"),
        "candidate_directory": str(candidate_directory),
        "relative_energy_reference": {
            "step_index": int(baseline_record["step_index"]),
            "backend_energy_hartree": float(baseline_record["backend_energy_hartree"]),
            "afir_energy_hartree": float(baseline_record["afir_energy_hartree"]),
            "total_energy_hartree": float(baseline_record["total_energy_hartree"]),
        },
        "candidate_files": candidate_files,
        "selected_indices": {
            "highest_backend_energy_index": highest_backend_index,
            "pre_product_index": pre_product_index,
            "first_topology_change_index": topology_change_index,
            "highest_total_energy_index": highest_total_index,
            "persistent_transition_index": transition_index,
        },
    }


def analyse_reaction_trace(job_directory, max_force=None, energy_outlier_z=None):
    """Write a compact summary and candidate geometries for one reaction path.

    ``job_directory`` may point either to the workflow root or directly to a
    ``reaction_trace/`` directory. The function writes ``path_summary.csv`` in
    the job root, exports the candidate XYZ files under ``candidate_ts/``, and
    returns a metadata dictionary describing the selected points.
    """
    job_directory, trace_directory = _resolve_trace_directory(job_directory)
    trace_records = _sorted_records(load_trace_records(trace_directory))
    if not trace_records:
        return None

    _validate_sorted_trace(trace_records)

    path_summary = job_directory / "path_summary.csv"
    candidate_directory = job_directory / "candidate_ts"
    candidate_directory.mkdir(parents=True, exist_ok=True)

    baseline_record = trace_records[0]

    with path_summary.open("w", newline="", encoding="utf-8") as fp:
        fieldnames = [
            "step_index",
            "backend_energy_hartree",
            "backend_relative_kcalmol",
            "bias_energy_hartree",
            "bias_relative_kcalmol",
            "afir_energy_hartree",
            "afir_relative_kcalmol",
            "total_energy_hartree",
            "total_relative_kcalmol",
            "backend_force_norm",
            "bias_force_norm",
            "afir_force_norm",
            "total_force_norm",
            "max_force",
            "backend_max_force",
            "min_interfragment_distance_angstrom",
            "collective_coordinate_bohr",
            "bias_controller_policy",
            "bias_potential",
            "bias_gamma_kj_mol",
            "bias_alpha_max_hartree_per_bohr",
            "bias_distance_power",
            "softmin_beta_per_bohr",
            "bias_alpha",
            "bias_alpha_critical",
            "bias_alpha_target",
            "bias_segment_alpha_critical",
            "bias_segment_alpha_target",
            "bias_energy_offset_hartree",
            "bias_segment_index",
            "contact_count",
            "minimum_contact_gap_bohr",
            "closest_contact_pair",
            "bond_change_count",
            "candidate_eligible",
            "candidate_exclusion_reason",
            "formed_bonds",
            "broken_bonds",
            "xyz_file",
        ]
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        exclusion_reasons = _candidate_exclusion_reasons(
            trace_records, max_force=max_force, energy_outlier_z=energy_outlier_z
        )
        for record_index, record in enumerate(trace_records):
            controller = record.get("bias_controller") or {}
            decision = controller.get("decision") or {}
            bias_parameters = record.get("bias_parameters") or {}
            row = {
                "step_index": int(record["step_index"]),
                "backend_energy_hartree": float(record["backend_energy_hartree"]),
                "backend_relative_kcalmol": _relative_kcalmol(
                    record["backend_energy_hartree"],
                    baseline_record["backend_energy_hartree"],
                ),
                "bias_energy_hartree": float(record["bias_energy_hartree"]),
                "bias_relative_kcalmol": _relative_kcalmol(
                    record["bias_energy_hartree"],
                    baseline_record["bias_energy_hartree"],
                ),
                "afir_energy_hartree": float(record["afir_energy_hartree"]),
                "afir_relative_kcalmol": _relative_kcalmol(
                    record["afir_energy_hartree"],
                    baseline_record["afir_energy_hartree"],
                ),
                "total_energy_hartree": float(record["total_energy_hartree"]),
                "total_relative_kcalmol": _relative_kcalmol(
                    record["total_energy_hartree"],
                    baseline_record["total_energy_hartree"],
                ),
                "backend_force_norm": record.get("backend_force_norm"),
                "bias_force_norm": record.get("bias_force_norm"),
                "afir_force_norm": record.get("afir_force_norm"),
                "total_force_norm": record.get("total_force_norm"),
                "max_force": record.get("max_force"),
                "backend_max_force": _physical_max_force(record),
                "min_interfragment_distance_angstrom": record.get(
                    "min_interfragment_distance_angstrom"
                ),
                "collective_coordinate_bohr": record.get("collective_coordinate_bohr"),
                "bias_controller_policy": decision.get("policy"),
                "bias_potential": bias_parameters.get("potential"),
                "bias_gamma_kj_mol": bias_parameters.get("gamma_kj_mol"),
                "bias_alpha_max_hartree_per_bohr": bias_parameters.get("alpha_max_hartree_per_bohr"),
                "bias_distance_power": bias_parameters.get("distance_power"),
                "softmin_beta_per_bohr": bias_parameters.get("beta_per_bohr"),
                "bias_alpha": decision.get("alpha"),
                "bias_alpha_critical": record.get("alpha_critical", decision.get("alpha_critical")),
                "bias_alpha_target": record.get("alpha_target", decision.get("alpha_target")),
                "bias_segment_alpha_critical": decision.get("alpha_critical"),
                "bias_segment_alpha_target": decision.get("alpha_target"),
                "bias_energy_offset_hartree": controller.get("energy_offset_hartree"),
                "bias_segment_index": controller.get("segment_index"),
                "contact_count": (record.get("contact_diagnostics") or {}).get("contact_count"),
                "minimum_contact_gap_bohr": (record.get("contact_diagnostics") or {}).get("minimum_gap_bohr"),
                "closest_contact_pair": json.dumps(
                    (record.get("contact_diagnostics") or {}).get("closest_pair"),
                    sort_keys=True,
                ),
                "bond_change_count": int(record["bond_change_count"]),
                "candidate_eligible": record_index not in exclusion_reasons,
                "candidate_exclusion_reason": ";".join(exclusion_reasons.get(record_index, [])),
                "formed_bonds": json.dumps(record.get("formed_bonds", []), sort_keys=True),
                "broken_bonds": json.dumps(record.get("broken_bonds", []), sort_keys=True),
                "xyz_file": _step_xyz_path(record),
            }
            writer.writerow(row)

    eligible_indices = _candidate_record_indices(
        trace_records, max_force=max_force, energy_outlier_z=energy_outlier_z
    )
    highest_backend_index = max(
        eligible_indices, key=lambda index: float(trace_records[index]["backend_energy_hartree"])
    )
    highest_total_index = max(
        eligible_indices, key=lambda index: float(trace_records[index]["total_energy_hartree"])
    )
    topology_change_index = _first_persistent_transition_index(trace_records, eligible_indices)
    transition_index = _persistent_transition_index(trace_records)
    pre_product_index = _pre_product_index(trace_records, eligible_indices)

    _write_xyz_record(
        candidate_directory / "highest_backend_energy.xyz",
        trace_records[highest_backend_index],
        "highest_backend_energy",
        reference_record=baseline_record,
        energy_key="backend_energy_hartree",
    )
    for label, index in (("pre_product_geometry", pre_product_index),
                         ("first_topology_change", topology_change_index)):
        path = candidate_directory / f"{label}.xyz"
        if index is None:
            # Reanalysis with stricter filters must not leave an old candidate.
            path.unlink(missing_ok=True)
        else:
            _write_xyz_record(path, trace_records[index], label, reference_record=baseline_record)
    _write_xyz_record(
        candidate_directory / "highest_total_energy.xyz",
        trace_records[highest_total_index],
        "highest_total_energy",
        reference_record=baseline_record,
        energy_key="total_energy_hartree",
    )

    metadata = _candidate_metadata(
        job_directory,
        trace_records,
        highest_backend_index,
        pre_product_index,
        topology_change_index,
        highest_total_index,
        transition_index,
        baseline_record,
    )
    metadata["candidate_filters"] = {
        "max_force_hartree_per_bohr": max_force,
        "max_force_metric": "backend_max_force or maximum norm of backend force vectors; missing physical force excluded",
        "backend_energy_outlier_robust_z": energy_outlier_z,
        "eligible_geometry_count": len(eligible_indices),
        "excluded_geometry_count": len(trace_records) - len(eligible_indices),
        "excluded_steps": {
            str(trace_records[index]["step_index"]): reasons
            for index, reasons in exclusion_reasons.items()
        },
    }
    metadata_path = candidate_directory / "metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    return {
        "path_summary_csv": str(path_summary),
        "candidate_ts_directory": str(candidate_directory),
        "candidate_metadata_json": str(metadata_path),
        "highest_backend_energy_index": highest_backend_index,
        "pre_product_index": pre_product_index,
        "first_topology_change_index": topology_change_index,
        "highest_total_energy_index": highest_total_index,
        "persistent_transition_index": transition_index,
        "eligible_geometry_count": len(eligible_indices),
        "excluded_geometry_count": len(trace_records) - len(eligible_indices),
    }


def _reaction_trace_plot_data(trace_records):
    """Return arrays used by the reaction trace plots."""
    baseline_record = trace_records[0]
    step_indices = np.asarray([int(record["step_index"]) for record in trace_records], dtype=int)
    backend_relative = np.asarray(
        [
            _relative_kcalmol(record["backend_energy_hartree"], baseline_record["backend_energy_hartree"])
            for record in trace_records
        ],
        dtype=float,
    )
    bias_relative = np.asarray(
        [
            _relative_kcalmol(record["bias_energy_hartree"], baseline_record["bias_energy_hartree"])
            for record in trace_records
        ],
        dtype=float,
    )
    total_relative = np.asarray(
        [
            _relative_kcalmol(record["total_energy_hartree"], baseline_record["total_energy_hartree"])
            for record in trace_records
        ],
        dtype=float,
    )
    bond_change_count = np.asarray(
        [int(record.get("bond_change_count", 0)) for record in trace_records],
        dtype=int,
    )
    min_interfragment_distance = np.asarray(
        [
            np.nan if record.get("min_interfragment_distance_angstrom") is None else float(record["min_interfragment_distance_angstrom"])
            for record in trace_records
        ],
        dtype=float,
    )
    backend_force_norm = np.asarray(
        [
            np.nan if record.get("backend_force_norm") is None else float(record["backend_force_norm"])
            for record in trace_records
        ],
        dtype=float,
    )
    bias_force_norm = np.asarray(
        [
            np.nan if record.get("bias_force_norm") is None else float(record["bias_force_norm"])
            for record in trace_records
        ],
        dtype=float,
    )
    total_force_norm = np.asarray(
        [
            np.nan if record.get("total_force_norm") is None else float(record["total_force_norm"])
            for record in trace_records
        ],
        dtype=float,
    )
    return {
        "step_indices": step_indices,
        "backend_relative": backend_relative,
        "bias_relative": bias_relative,
        "total_relative": total_relative,
        "bond_change_count": bond_change_count,
        "min_interfragment_distance": min_interfragment_distance,
        "backend_force_norm": backend_force_norm,
        "bias_force_norm": bias_force_norm,
        "total_force_norm": total_force_norm,
    }


def plot_reaction_trace(job_directory, output_directory=None):
    """Create PNG plots from a recorded reaction trace.

    The plots summarize the relative energies, bond-change count, inter-
    fragment distance, and force norms across the recorded steps. When no
    output directory is supplied, the figures are written to ``trace_plots/``
    in the job root.
    """
    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    job_directory, trace_directory = _resolve_trace_directory(job_directory)
    trace_records = _sorted_records(load_trace_records(trace_directory))
    if not trace_records:
        return None

    _validate_sorted_trace(trace_records)
    data = _reaction_trace_plot_data(trace_records)
    if output_directory is None:
        output_directory = job_directory / "trace_plots"
    output_directory = Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)

    plot_paths = {}

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(data["step_indices"], data["backend_relative"], marker="o", label="backend")
    ax.plot(data["step_indices"], data["bias_relative"], marker="o", label="bias")
    ax.plot(data["step_indices"], data["total_relative"], marker="o", label="total")
    ax.set_xlabel("Step index")
    ax.set_ylabel("Relative energy (kcal/mol)")
    ax.set_title("Reaction trace energies")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best")
    energy_path = output_directory / "reaction_trace_energies.png"
    fig.tight_layout()
    fig.savefig(energy_path, dpi=200)
    plt.close(fig)
    plot_paths["reaction_trace_energies"] = str(energy_path)

    fig, axes = plt.subplots(3, 1, figsize=(8, 9), sharex=True)
    axes[0].plot(data["step_indices"], data["bond_change_count"], marker="o", color="tab:purple")
    axes[0].set_ylabel("Bond changes")
    axes[0].set_title("Reaction trace diagnostics")
    axes[0].grid(True, alpha=0.25)

    axes[1].plot(
        data["step_indices"],
        data["min_interfragment_distance"],
        marker="o",
        color="tab:green",
    )
    axes[1].set_ylabel("Min interfragment distance (Å)")
    axes[1].grid(True, alpha=0.25)

    axes[2].plot(data["step_indices"], data["backend_force_norm"], marker="o", label="backend")
    axes[2].plot(data["step_indices"], data["bias_force_norm"], marker="o", label="bias")
    axes[2].plot(data["step_indices"], data["total_force_norm"], marker="o", label="total")
    axes[2].set_xlabel("Step index")
    axes[2].set_ylabel("Force norm")
    axes[2].grid(True, alpha=0.25)
    axes[2].legend(loc="best")

    metrics_path = output_directory / "reaction_trace_metrics.png"
    fig.tight_layout()
    fig.savefig(metrics_path, dpi=200)
    plt.close(fig)
    plot_paths["reaction_trace_metrics"] = str(metrics_path)

    return {
        "plot_directory": str(output_directory),
        "plots": plot_paths,
    }
