"""Reaction-path trace recording for geomeTRIC-backed AFIR optimizations."""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from pyar.backends import write_xyz
from pyar.data import new_atomic_data

_BOND_ON_SCALE = 1.15
_BOND_OFF_SCALE = 1.35
_VALID_TRACE_MODES = {"write", "append"}


def _symbol_radius(symbol):
    """Return the covalent radius in Angstrom for ``symbol``."""
    key = str(symbol).strip().capitalize()
    try:
        return float(new_atomic_data.covalent_radius[key])
    except KeyError as exc:
        raise KeyError(f"Unknown covalent radius for atomic symbol {symbol!r}") from exc


def infer_bonds(symbols, coordinates_angstrom, previous_bonds=None):
    """Infer a conservative bond set from coordinates and covalent radii."""
    coordinates = np.asarray(coordinates_angstrom, dtype=float)
    symbols = list(symbols)
    previous_bonds = set(previous_bonds or ())
    bonds = set()

    for i in range(len(symbols)):
        radius_i = _symbol_radius(symbols[i])
        for j in range(i + 1, len(symbols)):
            radius_j = _symbol_radius(symbols[j])
            radii_sum = radius_i + radius_j
            if radii_sum <= 0.0:
                continue
            distance = float(np.linalg.norm(coordinates[i] - coordinates[j]))
            ratio = distance / radii_sum
            pair = (i, j)
            if pair in previous_bonds:
                if ratio <= _BOND_OFF_SCALE:
                    bonds.add(pair)
            elif ratio <= _BOND_ON_SCALE:
                bonds.add(pair)

    return bonds


def bond_changes(previous_bonds, current_bonds):
    """Return formed and broken bonds between consecutive bond sets."""
    previous_bonds = set(previous_bonds or ())
    current_bonds = set(current_bonds or ())
    formed = sorted(current_bonds - previous_bonds)
    broken = sorted(previous_bonds - current_bonds)
    return formed, broken


def min_interfragment_distance(coordinates_angstrom, fragment_indices):
    """Return the minimum distance between atoms in distinct fragments."""
    if not fragment_indices or len(fragment_indices) < 2:
        return None

    coordinates = np.asarray(coordinates_angstrom, dtype=float)
    fragments = [list(fragment) for fragment in fragment_indices if fragment]
    if len(fragments) < 2:
        return None

    minimum = math.inf
    for left in range(len(fragments)):
        for right in range(left + 1, len(fragments)):
            for i in fragments[left]:
                for j in fragments[right]:
                    distance = float(np.linalg.norm(coordinates[i] - coordinates[j]))
                    if distance < minimum:
                        minimum = distance
    return None if not np.isfinite(minimum) else float(minimum)


def _read_xyz_frame(path):
    """Read one XYZ frame and return symbols plus coordinates."""
    lines = Path(path).read_text().splitlines()
    if len(lines) < 2:
        raise ValueError(f"Invalid XYZ trace file: {path}")

    natoms = int(lines[0].strip())
    if len(lines) < natoms + 2:
        raise ValueError(f"Invalid XYZ trace file: {path}")

    symbols = []
    coordinates = []
    for line in lines[2:2 + natoms]:
        parts = line.split()
        if len(parts) < 4:
            raise ValueError(f"Invalid XYZ trace file: {path}")
        symbols.append(parts[0])
        coordinates.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return symbols, np.asarray(coordinates, dtype=float)


def _discover_next_step_index(step_directory):
    """Return the next available step index from XYZ snapshots."""
    step_numbers = []
    for path in Path(step_directory).glob("step_*.xyz"):
        stem = path.stem
        try:
            step_numbers.append(int(stem.split("_")[-1]))
        except ValueError:
            continue
    return (max(step_numbers) + 1) if step_numbers else 0


class ReactionTraceRecorder:
    """Append geomeTRIC evaluation data to JSONL and XYZ trace files."""

    def __init__(self, job_directory, trace_name="reaction_trace", mode="write"):
        self.job_directory = Path(job_directory)
        self.trace_directory = self.job_directory / trace_name
        self.step_directory = self.trace_directory / "steps"
        self.trace_file = self.trace_directory / "trace.jsonl"
        self.mode = str(mode).lower()
        if self.mode not in _VALID_TRACE_MODES:
            raise ValueError(
                f"Unsupported trace mode {mode!r}; expected one of {sorted(_VALID_TRACE_MODES)}"
            )
        self.step_index = 0
        self.previous_bonds = None

        self.trace_directory.mkdir(parents=True, exist_ok=True)
        self.step_directory.mkdir(parents=True, exist_ok=True)
        if self.mode == "write":
            self.trace_file.unlink(missing_ok=True)
            for existing in self.step_directory.glob("step_*.xyz"):
                existing.unlink()
            return

        existing_records = load_trace_records(self.trace_file)
        if existing_records:
            last_record = max(
                existing_records,
                key=lambda record: int(record.get("step_index", -1)),
            )
            self.step_index = int(last_record.get("step_index", len(existing_records) - 1)) + 1
            previous_bonds = last_record.get("current_bonds") or []
            self.previous_bonds = {tuple(pair) for pair in previous_bonds}
            return

        self.step_index = _discover_next_step_index(self.step_directory)
        if self.step_index > 0:
            last_step_path = self.step_directory / f"step_{self.step_index - 1:06d}.xyz"
            try:
                symbols, coordinates = _read_xyz_frame(last_step_path)
            except (OSError, ValueError):
                self.previous_bonds = None
            else:
                self.previous_bonds = infer_bonds(symbols, coordinates)

    def record(
        self,
        *,
        symbols,
        coordinates_angstrom,
        backend_energy_hartree,
        afir_energy_hartree,
        total_energy_hartree,
        backend_forces_hartree_per_bohr,
        afir_forces_hartree_per_bohr,
        total_forces_hartree_per_bohr,
        backend_force_norm,
        afir_force_norm,
        total_force_norm,
        max_force,
        fragment_indices=None,
    ):
        """Write one trace record and its corresponding XYZ snapshot."""
        coordinates = np.asarray(coordinates_angstrom, dtype=float)
        backend_forces = np.asarray(backend_forces_hartree_per_bohr, dtype=float)
        afir_forces = np.asarray(afir_forces_hartree_per_bohr, dtype=float)
        total_forces = np.asarray(total_forces_hartree_per_bohr, dtype=float)
        current_bonds = infer_bonds(symbols, coordinates, self.previous_bonds)
        formed_bonds, broken_bonds = bond_changes(self.previous_bonds, current_bonds)
        min_distance = min_interfragment_distance(coordinates, fragment_indices)

        record = {
            "step_index": self.step_index,
            "symbols": list(symbols),
            "coordinates_angstrom": coordinates.tolist(),
            "backend_energy_hartree": float(backend_energy_hartree),
            "afir_energy_hartree": float(afir_energy_hartree),
            "total_energy_hartree": float(total_energy_hartree),
            "backend_force_norm": float(backend_force_norm),
            "afir_force_norm": float(afir_force_norm),
            "total_force_norm": float(total_force_norm),
            "max_force": float(max_force),
            "backend_forces_hartree_per_bohr": backend_forces.tolist(),
            "afir_forces_hartree_per_bohr": afir_forces.tolist(),
            "total_forces_hartree_per_bohr": total_forces.tolist(),
            "current_bonds": [list(pair) for pair in sorted(current_bonds)],
            "formed_bonds": [list(pair) for pair in formed_bonds],
            "broken_bonds": [list(pair) for pair in broken_bonds],
            "bond_change_count": len(formed_bonds) + len(broken_bonds),
            "min_interfragment_distance_angstrom": min_distance,
        }

        with self.trace_file.open("a", encoding="utf-8") as fp:
            json.dump(record, fp, sort_keys=True)
            fp.write("\n")

        write_xyz(
            list(symbols),
            coordinates,
            str(self.step_directory / f"step_{self.step_index:06d}.xyz"),
            job_name=f"step_{self.step_index:06d}",
            energy=float(total_energy_hartree),
        )

        self.previous_bonds = current_bonds
        self.step_index += 1
        return record


def load_trace_records(trace_file):
    """Load JSONL trace records from ``trace_file``."""
    trace_path = Path(trace_file)
    if trace_path.is_dir():
        trace_path = trace_path / "trace.jsonl"
    if not trace_path.exists():
        return []
    records = []
    with trace_path.open(encoding="utf-8") as fp:
        for line in fp:
            line = line.strip()
            if not line:
                continue
            records.append(json.loads(line))
    return records
