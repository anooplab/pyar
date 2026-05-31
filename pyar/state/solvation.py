"""Versioned, inspectable restart state for solvation workflows.

The solvation workflow persists a compact JSON state file plus geometry
snapshots so interrupted runs can resume without redoing completed cycles or
losing the selected seed geometries.
"""

from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path

import numpy as np

from pyar import __version__
from pyar.core.molecule import Molecule

STATE_VERSION = 1
STATE_FILENAME = "state.json"


class SolvationStateError(RuntimeError):
    """Raised when a solvation run cannot be safely resumed."""


def _json_value(value):
    """Return ``value`` in a JSON-serializable representation."""
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, dict):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_value(item) for item in value]
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    return str(value)


def _molecule_signature(molecule):
    """Return stable molecule metadata used to validate restarts."""
    return {
        "name": molecule.name,
        "atoms": list(molecule.atoms_list),
        "coordinates": np.asarray(molecule.coordinates, dtype=float).tolist(),
        "charge": molecule.charge,
        "multiplicity": molecule.multiplicity,
        "scftype": molecule.scftype,
        "fragment_definition": getattr(molecule, "fragments", []),
    }


class SolvationRunState:
    """Persist completed solvation cycles independently of optimizer objects.

    The state tracks the requested solvation calculation, the current seed
    geometries, the completed cycle history, and the terminal output paths.
    It is designed to be human-inspectable and safe to validate before a
    resumed run mutates the working directory.
    """

    def __init__(self, root_directory, data):
        """Bind the run state to the solvation directory tree on disk."""
        self.root_directory = Path(root_directory).resolve()
        self.solvation_directory = self.root_directory / "solvation"
        self.state_directory = self.solvation_directory / "state"
        self.geometry_directory = self.state_directory / "geometries"
        self.state_file = self.solvation_directory / STATE_FILENAME
        self.data = data

    @classmethod
    def create(cls, root_directory, request, current_seeds, sampling=None):
        """Create and persist state for a new solvation run.

        The initial state stores the request payload, snapshots of the input
        seed geometries, and the cycle counter used to drive the growth loop.
        """
        data = {
            "version": STATE_VERSION,
            "workflow": "solvation",
            "package_version": __version__,
            "status": "running",
            "request": _json_value(request),
            "next_cycle": 2,
            "completed_cycles": [],
            "current_seeds": [],
            "final_seeds": [],
            "next_snapshot": 0,
        }
        if sampling is not None:
            data["sampling"] = _json_value(sampling)
        state = cls(root_directory, data)
        state.geometry_directory.mkdir(parents=True, exist_ok=True)
        state._replace_current_seeds(current_seeds)
        state.save()
        return state

    @classmethod
    def load(cls, root_directory, expected_request):
        """Load a running solvation state file and validate its request.

        A matching request is required for restart. If the state file belongs
        to a different calculation, the load is rejected rather than silently
        mutating an unrelated run directory.
        """
        state_file = Path(root_directory).resolve() / "solvation" / STATE_FILENAME
        if not state_file.is_file():
            return None
        try:
            with state_file.open() as fp:
                data = json.load(fp)
        except (OSError, ValueError) as exc:
            raise SolvationStateError(
                f"Could not read solvation state file {state_file}: {exc}"
            ) from exc
        if data.get("version") != STATE_VERSION or data.get("workflow") != "solvation":
            raise SolvationStateError(f"Unsupported solvation state format in {state_file}")
        state = cls(root_directory, data)
        state.validate_progress()
        state.validate_request(expected_request)
        if data.get("status") != "running":
            raise SolvationStateError(
                f"Solvation state in {state_file} is already {data.get('status')!r}; "
                "start a new calculation in a new directory."
            )
        return state

    def validate_request(self, expected_request):
        """Reject resume attempts that change the solvation calculation."""
        if self.data.get("request") != _json_value(expected_request):
            raise SolvationStateError(
                "Existing solvation state does not match this invocation; "
                "resume with the original inputs and settings or use a new directory."
            )

    def validate_progress(self):
        """Reject incomplete or inconsistent persisted cycle progress."""
        completed = self.data.get("completed_cycles")
        if not isinstance(completed, list):
            raise SolvationStateError("Solvation state has invalid cycle progress data.")
        next_cycle = self.data.get("next_cycle")
        if not isinstance(next_cycle, int) or next_cycle < 2:
            raise SolvationStateError("Solvation state has an invalid next cycle marker.")
        aggregate_size = int(self.data.get("request", {}).get("aggregate_size", 0))
        if next_cycle > aggregate_size + 2:
            raise SolvationStateError("Solvation state records a cycle beyond the requested size.")

    def current_seed_molecules(self):
        """Load seed geometries to be used for the next solvation cycle."""
        return [self._load_molecule(reference) for reference in self.data.get("current_seeds", [])]

    def complete_cycle(self, cycle_number, selected_molecules):
        """Persist a completed cycle and update seeds for the next cycle.

        The selected molecules are snapshotted to disk so the next cycle can
        be resumed from an inspectable geometry set rather than from in-memory
        optimizer objects.
        """
        expected_cycle = self.data.get("next_cycle")
        if cycle_number != expected_cycle:
            raise SolvationStateError(
                "Solvation cycle completion is out of sequence; "
                "the saved state cannot be updated safely."
            )
        selected_refs = self._snapshot_molecules(selected_molecules)
        self.data["completed_cycles"].append(
            {
                "cycle": int(cycle_number),
                "status": "completed",
                "selected_results": [str(ref["path"]) for ref in selected_refs],
            }
        )
        self.data["current_seeds"] = selected_refs
        self.data["next_cycle"] = int(cycle_number) + 1
        self.save()

    def finish(self, final_status="completed"):
        """Persist terminal workflow state and final selected output paths."""
        aggregate_size = int(self.data.get("request", {}).get("aggregate_size", 0))
        if final_status == "completed" and self.data.get("next_cycle", 0) <= aggregate_size + 1:
            raise SolvationStateError(
                "Cannot complete solvation state while cycles remain unfinished."
            )
        self.data["status"] = final_status
        self.data["final_seeds"] = [str(ref["path"]) for ref in self.data.get("current_seeds", [])]
        self.save()

    def save(self):
        """Atomically write ``solvation/state.json``."""
        self.solvation_directory.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=self.solvation_directory,
            prefix=".state-",
            suffix=".json",
            delete=False,
        ) as fp:
            json.dump(self.data, fp, indent=2, sort_keys=True)
            fp.write("\n")
            temporary_path = Path(fp.name)
        os.replace(temporary_path, self.state_file)

    def _replace_current_seeds(self, molecules):
        """Replace the current seed list with snapshots of ``molecules``."""
        self.data["current_seeds"] = self._snapshot_molecules(molecules)

    def _snapshot_molecules(self, molecules):
        """Return persistent snapshot metadata for a sequence of molecules."""
        return [self._snapshot_molecule(molecule) for molecule in molecules]

    def _snapshot_molecule(self, molecule):
        """Write one molecule snapshot and return its persisted metadata."""
        sequence = self.data["next_snapshot"]
        self.data["next_snapshot"] += 1
        filename = f"{sequence:06d}_{_safe_name(molecule.name)}.xyz"
        output_path = self.geometry_directory / filename
        with tempfile.NamedTemporaryFile(
            dir=self.geometry_directory,
            prefix=".geometry-",
            suffix=".xyz",
            delete=False,
        ) as fp:
            temporary_path = Path(fp.name)
        molecule.mol_to_xyz(str(temporary_path))
        os.replace(temporary_path, output_path)
        return {
            "path": str(output_path.relative_to(self.solvation_directory)),
            "name": molecule.name,
            "title": getattr(molecule, "title", "Title"),
            "fragments": _json_value(getattr(molecule, "fragments", [])),
            "charge": _json_value(molecule.charge),
            "multiplicity": _json_value(molecule.multiplicity),
            "scftype": molecule.scftype,
            "energy": _json_value(getattr(molecule, "energy", None)),
        }

    def _load_molecule(self, reference):
        """Load one snapshot from disk and reconstruct the molecule object."""
        source = self.solvation_directory / reference["path"]
        try:
            loaded = Molecule.from_xyz(str(source))
        except (OSError, ValueError, SystemExit) as exc:
            raise SolvationStateError(
                f"Could not restore solvation geometry snapshot {source}: {exc}"
            ) from exc
        return Molecule(
            loaded.atoms_list,
            loaded.coordinates,
            name=reference["name"],
            title=reference["title"],
            fragments=reference.get("fragments", []),
            charge=reference.get("charge", 0),
            multiplicity=reference.get("multiplicity", 1),
            scftype=reference.get("scftype", "rhf"),
            energy=reference.get("energy"),
        )


def _safe_name(value):
    """Return a filesystem-safe fragment of a molecule name."""
    return "".join(ch if ch.isalnum() or ch in "_.-" else "_" for ch in str(value)).strip("_") or "molecule"
