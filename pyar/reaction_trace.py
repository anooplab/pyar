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
    "bias_energy_hartree",
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
    record = dict(record)
    if "bias_energy_hartree" not in record and "afir_energy_hartree" in record:
        record["bias_energy_hartree"] = record["afir_energy_hartree"]
    if "bias_force_norm" not in record and "afir_force_norm" in record:
        record["bias_force_norm"] = record["afir_force_norm"]
    if "bias_forces_hartree_per_bohr" not in record and "afir_forces_hartree_per_bohr" in record:
        record["bias_forces_hartree_per_bohr"] = record["afir_forces_hartree_per_bohr"]

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
    if "bias_energy_hartree" not in record:
        raise ValueError(f"Trace record {record_index} is missing required key 'bias_energy_hartree'")
    _ensure_finite_scalar(record["bias_energy_hartree"], "bias_energy_hartree", record_index)
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
    if "bias_force_norm" in record:
        _ensure_finite_scalar(record["bias_force_norm"], "bias_force_norm", record_index)
    if "total_force_norm" in record:
        _ensure_finite_scalar(record["total_force_norm"], "total_force_norm", record_index)
    if "max_force" in record:
        _ensure_finite_scalar(record["max_force"], "max_force", record_index)
    if "backend_max_force" in record:
        _ensure_finite_scalar(record["backend_max_force"], "backend_max_force", record_index)
    if "collective_coordinate_bohr" in record:
        _ensure_finite_scalar(record["collective_coordinate_bohr"], "collective_coordinate_bohr", record_index)
    if "softmin_beta" in record:
        beta = _ensure_finite_scalar(record["softmin_beta"], "softmin_beta", record_index)
        if beta <= 0.0:
            raise ValueError(f"Trace record {record_index} field 'softmin_beta' must be positive")
    if "contact_diagnostics" in record and not isinstance(record["contact_diagnostics"], dict):
        raise ValueError(f"Trace record {record_index} field 'contact_diagnostics' must be an object")
    if "bias_parameters" in record and not isinstance(record["bias_parameters"], dict):
        raise ValueError(f"Trace record {record_index} field 'bias_parameters' must be an object")

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
    if "bias_forces_hartree_per_bohr" in record:
        bias_forces = np.asarray(record["bias_forces_hartree_per_bohr"], dtype=float)
        if bias_forces.shape != (len(normalized_symbols), 3):
            raise ValueError(
                f"Trace record {record_index} field 'bias_forces_hartree_per_bohr' must have "
                f"shape ({len(normalized_symbols)}, 3); got {bias_forces.shape!r}"
            )
        if not np.all(np.isfinite(bias_forces)):
            raise ValueError(
                f"Trace record {record_index} field 'bias_forces_hartree_per_bohr' must contain "
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
        "bias_energy_hartree": float(record["bias_energy_hartree"]),
        "afir_energy_hartree": float(record["bias_energy_hartree"]),
        "total_energy_hartree": float(record["total_energy_hartree"]),
        "backend_force_norm": float(record["backend_force_norm"]) if "backend_force_norm" in record else None,
        "bias_force_norm": float(record["bias_force_norm"]) if "bias_force_norm" in record else None,
        "afir_force_norm": float(record["bias_force_norm"]) if "bias_force_norm" in record else None,
        "total_force_norm": float(record["total_force_norm"]) if "total_force_norm" in record else None,
        "max_force": float(record["max_force"]) if "max_force" in record else None,
        "backend_max_force": (
            float(record["backend_max_force"]) if "backend_max_force" in record else None
        ),
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
    if "bias_forces_hartree_per_bohr" in record:
        normalized["bias_forces_hartree_per_bohr"] = bias_forces.tolist()
        normalized["afir_forces_hartree_per_bohr"] = bias_forces.tolist()
    if "total_forces_hartree_per_bohr" in record:
        normalized["total_forces_hartree_per_bohr"] = total_forces.tolist()
    if "collective_coordinate_bohr" in record:
        normalized["collective_coordinate_bohr"] = float(record["collective_coordinate_bohr"])
    if "softmin_beta" in record:
        normalized["softmin_beta"] = float(record["softmin_beta"])
    if "contact_diagnostics" in record:
        normalized["contact_diagnostics"] = record["contact_diagnostics"]
    if "bias_controller" in record:
        if not isinstance(record["bias_controller"], dict):
            raise ValueError(f"Trace record {record_index} field 'bias_controller' must be an object")
        normalized["bias_controller"] = record["bias_controller"]
    if "bias_parameters" in record:
        normalized["bias_parameters"] = record["bias_parameters"]
    for key in (
        "accepted_step", "segment_index", "release_state", "release_reason",
        "release_evidence", "forming_pairs", "forming_pair_distance",
        "normalized_distance", "topology_change_state", "persistence_counter",
        "alpha", "alpha_critical", "alpha_target", "bond_order_available",
        "bond_order_scheme", "bond_order_delta", "bond_order_stabilized",
        "comparable_bond_orders",
    ):
        if key in record:
            normalized[key] = record[key]

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
        bias_energy_hartree=None,
        total_energy_hartree,
        backend_forces_hartree_per_bohr,
        bias_forces_hartree_per_bohr=None,
        total_forces_hartree_per_bohr,
        backend_force_norm,
        bias_force_norm=None,
        total_force_norm,
        max_force,
        backend_max_force=None,
        fragment_indices=None,
        collective_coordinate_bohr=None,
        contact_diagnostics=None,
        softmin_beta=None,
        bias_controller=None,
        bias_parameters=None,
        afir_energy_hartree=None,
        afir_forces_hartree_per_bohr=None,
        afir_force_norm=None,
        accepted_step=None,
        segment_index=None,
        release_state=None,
        release_reason=None,
        release_evidence=None,
        forming_pairs=None,
        forming_pair_distance=None,
        normalized_distance=None,
        topology_change_state=None,
        persistence_counter=None,
        alpha=None,
        alpha_critical=None,
        alpha_target=None,
        bond_order_available=None,
        bond_order_scheme=None,
        bond_order_delta=None,
        bond_order_stabilized=None,
        comparable_bond_orders=None,
    ):
        """Write one trace record and its corresponding XYZ snapshot.

        The recorded JSON object includes the raw energies/forces, bond-change
        counts, and interfragment distance used by the downstream analysis
        tools. The companion XYZ snapshot preserves the optimized geometry for
        the same step index.
        """
        coordinates = np.asarray(coordinates_angstrom, dtype=float)
        backend_forces = np.asarray(backend_forces_hartree_per_bohr, dtype=float)
        if bias_energy_hartree is None:
            bias_energy_hartree = afir_energy_hartree
        if bias_forces_hartree_per_bohr is None:
            bias_forces_hartree_per_bohr = afir_forces_hartree_per_bohr
        if bias_force_norm is None:
            bias_force_norm = afir_force_norm
        if bias_energy_hartree is None or bias_forces_hartree_per_bohr is None or bias_force_norm is None:
            raise ValueError("reaction trace requires bias energy, forces, and force norm")
        bias_forces = np.asarray(bias_forces_hartree_per_bohr, dtype=float)
        total_forces = np.asarray(total_forces_hartree_per_bohr, dtype=float)
        current_bonds = infer_bonds(symbols, coordinates, self.previous_bonds)
        formed_bonds, broken_bonds = bond_changes(self.previous_bonds, current_bonds)
        min_distance = min_interfragment_distance(coordinates, fragment_indices)

        record = {
            "step_index": self.step_index,
            "symbols": list(symbols),
            "coordinates_angstrom": coordinates.tolist(),
            "backend_energy_hartree": float(backend_energy_hartree),
            "bias_energy_hartree": float(bias_energy_hartree),
            "afir_energy_hartree": float(bias_energy_hartree),
            "total_energy_hartree": float(total_energy_hartree),
            "backend_force_norm": float(backend_force_norm),
            "bias_force_norm": float(bias_force_norm),
            "afir_force_norm": float(bias_force_norm),
            "total_force_norm": float(total_force_norm),
            "max_force": float(max_force),
            "backend_forces_hartree_per_bohr": backend_forces.tolist(),
            "bias_forces_hartree_per_bohr": bias_forces.tolist(),
            "afir_forces_hartree_per_bohr": bias_forces.tolist(),
            "total_forces_hartree_per_bohr": total_forces.tolist(),
            "current_bonds": [list(pair) for pair in sorted(current_bonds)],
            "formed_bonds": [list(pair) for pair in formed_bonds],
            "broken_bonds": [list(pair) for pair in broken_bonds],
            "bond_change_count": len(formed_bonds) + len(broken_bonds),
            "min_interfragment_distance_angstrom": min_distance,
        }
        if backend_max_force is not None:
            record["backend_max_force"] = float(backend_max_force)
        if collective_coordinate_bohr is not None:
            record["collective_coordinate_bohr"] = float(collective_coordinate_bohr)
        if contact_diagnostics is not None:
            record["contact_diagnostics"] = contact_diagnostics
        if softmin_beta is not None:
            record["softmin_beta"] = float(softmin_beta)
        if bias_controller is not None:
            record["bias_controller"] = bias_controller
        if bias_parameters is not None:
            record["bias_parameters"] = bias_parameters
        for key, value in (
            ("accepted_step", accepted_step),
            ("segment_index", segment_index),
            ("release_state", release_state),
            ("release_reason", release_reason),
            ("release_evidence", release_evidence),
            ("forming_pairs", forming_pairs),
            ("forming_pair_distance", forming_pair_distance),
            ("normalized_distance", normalized_distance),
            ("topology_change_state", topology_change_state),
            ("persistence_counter", persistence_counter),
            ("alpha", alpha),
            ("alpha_critical", alpha_critical),
            ("alpha_target", alpha_target),
            ("bond_order_available", bond_order_available),
            ("bond_order_scheme", bond_order_scheme),
            ("bond_order_delta", bond_order_delta),
            ("bond_order_stabilized", bond_order_stabilized),
            ("comparable_bond_orders", comparable_bond_orders),
        ):
            if value is not None:
                record[key] = value

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
