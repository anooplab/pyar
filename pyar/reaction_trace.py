"""Reaction-path trace recording for geomeTRIC-backed AFIR optimizations.

The recorder writes a JSONL trace plus per-step XYZ snapshots for each
geomeTRIC-backed reaction evaluation. The resulting trace is the input for
the analysis and plotting helpers in :mod:`pyar.reaction_analysis`.
"""

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
_REQUIRED_TRACE_KEYS = {
    "step_index",
    "symbols",
    "coordinates_angstrom",
    "backend_energy_hartree",
    "afir_energy_hartree",
    "total_energy_hartree",
    "current_bonds",
    "formed_bonds",
    "broken_bonds",
    "bond_change_count",
    "min_interfragment_distance_angstrom",
}


def _symbol_radius(symbol):
    """Return the covalent radius in Angstrom for ``symbol``."""
    key = str(symbol).strip().capitalize()
    try:
        return float(new_atomic_data.covalent_radius[key])
    except KeyError as exc:
        raise KeyError(f"Unknown covalent radius for atomic symbol {symbol!r}") from exc


def infer_bonds(symbols, coordinates_angstrom, previous_bonds=None):
    """Infer a conservative bond set from coordinates and covalent radii.

    The heuristic is intentionally cautious: new bonds are only introduced
    when the interatomic distance is comfortably below the covalent-radius
    threshold, while previously observed bonds are allowed a wider hysteresis
    window so small numerical changes do not flicker the bond list.
    """
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
    """Return the minimum distance between atoms in distinct fragments.

    The value is only computed when at least two non-empty fragments are
    available. Otherwise the trace record stores ``None`` for this field.
    """
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


def _ensure_finite_scalar(value, field_name, record_index):
    """Return ``value`` as float or raise a clear validation error."""
    try:
        scalar = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"Trace record {record_index} field {field_name!r} must be a finite scalar; "
            f"got {value!r}"
        ) from exc
    if not np.isfinite(scalar):
        raise ValueError(
            f"Trace record {record_index} field {field_name!r} must be a finite scalar; "
            f"got {value!r}"
        )
    return scalar


def _validate_bond_list(value, field_name, natoms, record_index):
    """Validate a bond list of zero-based atom index pairs."""
    if value is None:
        raise ValueError(f"Trace record {record_index} field {field_name!r} is missing")
    if not isinstance(value, list):
        raise ValueError(
            f"Trace record {record_index} field {field_name!r} must be a list of pairs; "
            f"got {type(value).__name__}"
        )
    bonds = []
    for pair in value:
        if not isinstance(pair, (list, tuple)) or len(pair) != 2:
            raise ValueError(
                f"Trace record {record_index} field {field_name!r} must contain 2-item pairs; "
                f"got {pair!r}"
            )
        try:
            i = int(pair[0])
            j = int(pair[1])
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"Trace record {record_index} field {field_name!r} must contain integer pairs; "
                f"got {pair!r}"
            ) from exc
        if i < 0 or j < 0 or i >= natoms or j >= natoms or i >= j:
            raise ValueError(
                f"Trace record {record_index} field {field_name!r} contains invalid pair {pair!r} "
                f"for a system with {natoms} atoms"
            )
        bonds.append((i, j))
    return bonds


def validate_trace_record(record, record_index):
    """Validate the structure of one trace record.

    The validator enforces the schema used by the analysis and plotting code
    and returns a normalized record with numeric fields converted to native
    Python types.
    """
    if not isinstance(record, dict):
        raise ValueError(
            f"Trace record {record_index} must be a JSON object; got {type(record).__name__}"
        )

    missing = sorted(_REQUIRED_TRACE_KEYS - record.keys())
    if missing:
        raise ValueError(
            f"Trace record {record_index} is missing required fields: {', '.join(missing)}"
        )

    step_index = int(record["step_index"])
    if step_index < 0:
        raise ValueError(f"Trace record {record_index} has negative step_index {step_index}")

    symbols = record["symbols"]
    if not isinstance(symbols, list) or not symbols:
        raise ValueError(
            f"Trace record {record_index} field 'symbols' must be a non-empty list"
        )
    normalized_symbols = [str(symbol) for symbol in symbols]

    coordinates = np.asarray(record["coordinates_angstrom"], dtype=float)
    if coordinates.ndim != 2 or coordinates.shape != (len(normalized_symbols), 3):
        raise ValueError(
            f"Trace record {record_index} field 'coordinates_angstrom' must have shape "
            f"({len(normalized_symbols)}, 3); got {coordinates.shape!r}"
        )
    if not np.all(np.isfinite(coordinates)):
        raise ValueError(
            f"Trace record {record_index} field 'coordinates_angstrom' must contain finite values"
        )

    _ensure_finite_scalar(record["backend_energy_hartree"], "backend_energy_hartree", record_index)
    _ensure_finite_scalar(record["afir_energy_hartree"], "afir_energy_hartree", record_index)
    _ensure_finite_scalar(record["total_energy_hartree"], "total_energy_hartree", record_index)
    _ensure_finite_scalar(record["bond_change_count"], "bond_change_count", record_index)

    min_distance = record["min_interfragment_distance_angstrom"]
    if min_distance is not None:
        _ensure_finite_scalar(
            min_distance,
            "min_interfragment_distance_angstrom",
            record_index,
        )

    if "backend_force_norm" in record:
        _ensure_finite_scalar(record["backend_force_norm"], "backend_force_norm", record_index)
    if "afir_force_norm" in record:
        _ensure_finite_scalar(record["afir_force_norm"], "afir_force_norm", record_index)
    if "total_force_norm" in record:
        _ensure_finite_scalar(record["total_force_norm"], "total_force_norm", record_index)
    if "max_force" in record:
        _ensure_finite_scalar(record["max_force"], "max_force", record_index)

    if "backend_forces_hartree_per_bohr" in record:
        backend_forces = np.asarray(record["backend_forces_hartree_per_bohr"], dtype=float)
        if backend_forces.shape != (len(normalized_symbols), 3):
            raise ValueError(
                f"Trace record {record_index} field 'backend_forces_hartree_per_bohr' must have "
                f"shape ({len(normalized_symbols)}, 3); got {backend_forces.shape!r}"
            )
        if not np.all(np.isfinite(backend_forces)):
            raise ValueError(
                f"Trace record {record_index} field 'backend_forces_hartree_per_bohr' must "
                "contain finite values"
            )
    if "afir_forces_hartree_per_bohr" in record:
        afir_forces = np.asarray(record["afir_forces_hartree_per_bohr"], dtype=float)
        if afir_forces.shape != (len(normalized_symbols), 3):
            raise ValueError(
                f"Trace record {record_index} field 'afir_forces_hartree_per_bohr' must have "
                f"shape ({len(normalized_symbols)}, 3); got {afir_forces.shape!r}"
            )
        if not np.all(np.isfinite(afir_forces)):
            raise ValueError(
                f"Trace record {record_index} field 'afir_forces_hartree_per_bohr' must contain "
                "finite values"
            )
    if "total_forces_hartree_per_bohr" in record:
        total_forces = np.asarray(record["total_forces_hartree_per_bohr"], dtype=float)
        if total_forces.shape != (len(normalized_symbols), 3):
            raise ValueError(
                f"Trace record {record_index} field 'total_forces_hartree_per_bohr' must have "
                f"shape ({len(normalized_symbols)}, 3); got {total_forces.shape!r}"
            )
        if not np.all(np.isfinite(total_forces)):
            raise ValueError(
                f"Trace record {record_index} field 'total_forces_hartree_per_bohr' must contain "
                "finite values"
            )

    current_bonds = _validate_bond_list(record["current_bonds"], "current_bonds", len(normalized_symbols), record_index)
    formed_bonds = _validate_bond_list(record["formed_bonds"], "formed_bonds", len(normalized_symbols), record_index)
    broken_bonds = _validate_bond_list(record["broken_bonds"], "broken_bonds", len(normalized_symbols), record_index)
    bond_change_count = int(float(record["bond_change_count"]))
    expected_change_count = len(formed_bonds) + len(broken_bonds)
    if bond_change_count != expected_change_count:
        raise ValueError(
            f"Trace record {record_index} field 'bond_change_count' must equal the number of "
            f"formed and broken bonds ({expected_change_count}); got {bond_change_count}"
        )

    normalized = {
        "step_index": step_index,
        "symbols": normalized_symbols,
        "coordinates_angstrom": coordinates.tolist(),
        "backend_energy_hartree": float(record["backend_energy_hartree"]),
        "afir_energy_hartree": float(record["afir_energy_hartree"]),
        "total_energy_hartree": float(record["total_energy_hartree"]),
        "backend_force_norm": float(record["backend_force_norm"]) if "backend_force_norm" in record else None,
        "afir_force_norm": float(record["afir_force_norm"]) if "afir_force_norm" in record else None,
        "total_force_norm": float(record["total_force_norm"]) if "total_force_norm" in record else None,
        "max_force": float(record["max_force"]) if "max_force" in record else None,
        "current_bonds": [list(pair) for pair in current_bonds],
        "formed_bonds": [list(pair) for pair in formed_bonds],
        "broken_bonds": [list(pair) for pair in broken_bonds],
        "bond_change_count": bond_change_count,
        "min_interfragment_distance_angstrom": (
            None if min_distance is None else float(min_distance)
        ),
    }

    if "backend_forces_hartree_per_bohr" in record:
        normalized["backend_forces_hartree_per_bohr"] = backend_forces.tolist()
    if "afir_forces_hartree_per_bohr" in record:
        normalized["afir_forces_hartree_per_bohr"] = afir_forces.tolist()
    if "total_forces_hartree_per_bohr" in record:
        normalized["total_forces_hartree_per_bohr"] = total_forces.tolist()

    return normalized


class ReactionTraceRecorder:
    """Append geomeTRIC evaluation data to JSONL and XYZ trace files.

    The recorder is stateful: it tracks the next step index, reconstructs the
    previous bond set when appending to an existing trace, and keeps the JSONL
    trace and XYZ snapshots synchronized on disk.
    """

    def __init__(self, job_directory, trace_name="reaction_trace", mode="write"):
        """Initialize the trace recorder.

        ``mode='write'`` starts a fresh trace tree, while ``mode='append'``
        resumes an existing trace and continues numbering from the last step
        that was already recorded.
        """
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
        """Write one trace record and its corresponding XYZ snapshot.

        The recorded JSON object includes the raw energies/forces, bond-change
        counts, and interfragment distance used by the downstream analysis
        tools. The companion XYZ snapshot preserves the optimized geometry for
        the same step index.
        """
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
    """Load and validate JSONL trace records from ``trace_file``.

    ``trace_file`` may be either the JSONL file itself or the containing
    ``reaction_trace/`` directory. Each line is validated and normalized
    before being returned, so callers can rely on the trace schema.
    """
    trace_path = Path(trace_file)
    if trace_path.is_dir():
        trace_path = trace_path / "trace.jsonl"
    if not trace_path.exists():
        return []
    records = []
    with trace_path.open(encoding="utf-8") as fp:
        for line_number, line in enumerate(fp, start=1):
            line = line.strip()
            if not line:
                continue
            try:
                record = json.loads(line)
            except json.JSONDecodeError as exc:
                raise ValueError(
                    f"Could not parse JSON trace record at line {line_number} in {trace_path}: {exc}"
                ) from exc
            records.append(validate_trace_record(record, line_number))
    return records
