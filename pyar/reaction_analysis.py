"""Simple AFIR reaction-path analysis built on recorded geomeTRIC traces."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from pyar.backends import write_xyz
from pyar.reaction_trace import load_trace_records


def _sorted_records(records):
    """Return records ordered by their recorded step index."""
    return sorted(records, key=lambda record: int(record.get("step_index", 0)))


def _write_xyz_record(path, record, energy_key):
    """Write one trace record as an XYZ file."""
    write_xyz(
        record["symbols"],
        np.asarray(record["coordinates_angstrom"], dtype=float),
        str(path),
        job_name=Path(path).stem,
        energy=float(record[energy_key]),
    )


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
    """Return the first bond-change step that persists beyond one frame."""
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


def _pre_product_index(records):
    """Return the geometry immediately preceding the first persistent bond event."""
    transition_index = _persistent_transition_index(records)
    if transition_index <= 0:
        return _record_index_with_max(records, "backend_energy_hartree")
    return transition_index - 1


def analyse_reaction_trace(job_directory):
    """Write a compact summary and TS-candidate geometries for one path."""
    job_directory = Path(job_directory)
    trace_records = _sorted_records(load_trace_records(job_directory / "reaction_trace"))
    if not trace_records:
        return None

    path_summary = job_directory / "path_summary.csv"
    candidate_directory = job_directory / "candidate_ts"
    candidate_directory.mkdir(parents=True, exist_ok=True)

    with path_summary.open("w", newline="", encoding="utf-8") as fp:
        fieldnames = [
            "step_index",
            "backend_energy_hartree",
            "afir_energy_hartree",
            "total_energy_hartree",
            "backend_force_norm",
            "afir_force_norm",
            "total_force_norm",
            "max_force",
            "min_interfragment_distance_angstrom",
            "bond_change_count",
            "formed_bonds",
            "broken_bonds",
        ]
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        for record in trace_records:
            writer.writerow({key: record.get(key) for key in fieldnames})

    highest_backend_index = _record_index_with_max(trace_records, "backend_energy_hartree")
    highest_total_index = _record_index_with_max(trace_records, "total_energy_hartree")
    max_bond_change_index = _record_index_with_max(trace_records, "bond_change_count")
    pre_product_index = _pre_product_index(trace_records)
    transition_index = _persistent_transition_index(trace_records)

    _write_xyz_record(
        candidate_directory / "highest_backend_energy.xyz",
        trace_records[highest_backend_index],
        "backend_energy_hartree",
    )
    _write_xyz_record(
        candidate_directory / "pre_product_geometry.xyz",
        trace_records[pre_product_index],
        "total_energy_hartree",
    )
    _write_xyz_record(
        candidate_directory / "max_bond_change.xyz",
        trace_records[max_bond_change_index],
        "total_energy_hartree",
    )
    _write_xyz_record(
        candidate_directory / "highest_total_energy.xyz",
        trace_records[highest_total_index],
        "total_energy_hartree",
    )

    return {
        "path_summary_csv": str(path_summary),
        "candidate_ts_directory": str(candidate_directory),
        "highest_backend_energy_index": highest_backend_index,
        "pre_product_index": pre_product_index,
        "max_bond_change_index": max_bond_change_index,
        "highest_total_energy_index": highest_total_index,
        "persistent_transition_index": transition_index,
    }
