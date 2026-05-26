"""Versioned, inspectable restart state for reaction workflows."""

from __future__ import annotations

import json
import os
import pickle
import re
import tempfile
from pathlib import Path

import numpy as np

from pyar import __version__
from pyar.core.molecule import Molecule

STATE_VERSION = 2
STATE_FILENAME = "state.json"


class ReactionStateError(RuntimeError):
    """Raised when a reaction run cannot be created or safely resumed."""


def _json_value(value):
    """Return ``value`` in a JSON-serializable representation."""
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_value(item) for item in value]
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    return str(value)


def _safe_name(value):
    """Return a filesystem-safe fragment of a molecule or job name."""
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_") or "molecule"


class ReactionRunState:
    """Persist reaction workflow progress independently of optimizer objects."""

    def __init__(self, root_directory, data):
        self.root_directory = Path(root_directory).resolve()
        self.reaction_directory = self.root_directory / "reaction"
        self.state_directory = self.reaction_directory / "state"
        self.geometry_directory = self.state_directory / "geometries"
        self.state_file = self.reaction_directory / STATE_FILENAME
        self.data = data

    @classmethod
    def create(cls, root_directory, request, orientations, reactants):
        """Create and persist state for a fresh reaction calculation."""
        data = {
            "version": STATE_VERSION,
            "workflow": "reaction",
            "package_version": __version__,
            "status": "running",
            "request": _json_value(request),
            "gamma_schedule": [float(value) for value in request["gamma_schedule"]],
            "current_cycle": 0,
            "completed_cycles": [],
            "completed_jobs": [],
            "products": [],
            "pending_orientations": [],
            "current_survivors": [],
            "input_reactants": [],
            "next_snapshot": 0,
        }
        state = cls(root_directory, data)
        state.geometry_directory.mkdir(parents=True, exist_ok=True)
        state.data["input_reactants"] = [
            state._snapshot_molecule(molecule, "input")
            for molecule in reactants
        ]
        state._replace_pending_orientations(orientations)
        state.save()
        return state

    @classmethod
    def load(cls, root_directory, expected_request):
        """Load a current-format reaction state file and validate its request."""
        state_file = Path(root_directory).resolve() / "reaction" / STATE_FILENAME
        if not state_file.is_file():
            return None
        try:
            with state_file.open() as fp:
                data = json.load(fp)
        except (OSError, ValueError) as exc:
            raise ReactionStateError(f"Could not read reaction state file {state_file}: {exc}") from exc
        if data.get("version") != STATE_VERSION:
            raise ReactionStateError(
                f"Unsupported reaction state version {data.get('version')!r} in {state_file}"
            )
        state = cls(root_directory, data)
        state.validate_request(expected_request)
        if data.get("status") != "running":
            raise ReactionStateError(
                f"Reaction state in {state_file} is already {data.get('status')!r}; "
                "start a new calculation in a new directory."
            )
        return state

    @classmethod
    def migrate_legacy(cls, root_directory, checkpoint, request):
        """Convert an unambiguous legacy ``jobs.pkl`` checkpoint to JSON state."""
        requested_schedule = [float(value) for value in request["gamma_schedule"]]
        legacy_label_map = {}
        for gamma in requested_schedule:
            legacy_label = f"{int(gamma):04d}"
            if legacy_label in legacy_label_map:
                raise ReactionStateError(
                    "Legacy jobs.pkl cannot be safely resumed: multiple requested "
                    "gamma values share one old checkpoint label."
                )
            legacy_label_map[legacy_label] = gamma

        remaining_schedule = []
        for key in checkpoint:
            label = str(key)
            if label in legacy_label_map:
                remaining_schedule.append(legacy_label_map[label])
                continue
            try:
                numeric_key = float(key)
            except (TypeError, ValueError) as exc:
                raise ReactionStateError(
                    f"Legacy jobs.pkl contains an unreadable gamma key: {key!r}"
                ) from exc
            if numeric_key not in requested_schedule:
                raise ReactionStateError(
                    "Legacy jobs.pkl does not match the requested gamma schedule."
                )
            remaining_schedule.append(numeric_key)

        if not remaining_schedule:
            raise ReactionStateError("Legacy jobs.pkl does not contain pending reaction cycles.")
        start_index = requested_schedule.index(remaining_schedule[0])
        if remaining_schedule != requested_schedule[start_index:]:
            raise ReactionStateError(
                "Legacy jobs.pkl does not describe a valid remaining gamma schedule."
            )

        data = {
            "version": STATE_VERSION,
            "workflow": "reaction",
            "package_version": __version__,
            "status": "running",
            "request": _json_value(request),
            "gamma_schedule": requested_schedule,
            "current_cycle": start_index,
            "completed_cycles": [
                {"gamma": gamma, "status": "imported_completed"}
                for gamma in requested_schedule[:start_index]
            ],
            "completed_jobs": [],
            "products": [],
            "pending_orientations": [],
            "current_survivors": [],
            "input_reactants": [],
            "next_snapshot": 0,
            "legacy_import": "jobs.pkl",
        }
        state = cls(root_directory, data)
        state.geometry_directory.mkdir(parents=True, exist_ok=True)
        state._replace_pending_orientations(checkpoint[next(iter(checkpoint))])
        state.save()
        return state

    def validate_request(self, expected_request):
        """Fail clearly when a resume invocation changes scientific settings."""
        if self.data.get("request") != _json_value(expected_request):
            raise ReactionStateError(
                "Existing reaction state does not match this invocation; "
                "resume with the original inputs and settings or use a new directory."
            )

    def remaining_gamma_schedule(self):
        """Return the numeric gamma values still to be processed."""
        return self.data["gamma_schedule"][self.data["current_cycle"]:]

    def pending_molecules(self):
        """Load the geometries still pending in the current gamma cycle."""
        return [self._load_molecule(reference) for reference in self.data.get("pending_orientations", [])]

    def current_survivor_molecules(self):
        """Load successful current-cycle candidates retained before interruption."""
        return [self._load_molecule(reference) for reference in self.data.get("current_survivors", [])]

    def saved_product_identities(self):
        """Return saved product identifiers for restart-time deduplication."""
        return {
            product["job_name"]: (product["inchi"], product["smiles"])
            for product in self.data.get("products", [])
        }

    def record_job(self, job_name, gamma, status, remaining_orientations, current_survivors):
        """Record processed work plus pending and retained current-cycle candidates."""
        self.data["completed_jobs"].append(
            {"job_name": job_name, "gamma": float(gamma), "status": str(status)}
        )
        self._replace_pending_orientations(remaining_orientations)
        self._replace_current_survivors(current_survivors)
        self.save()

    def record_product(self, job_name, gamma, inchi, smiles, path):
        """Record and immediately persist one newly discovered product."""
        self.data["products"].append(
            {
                "job_name": job_name,
                "gamma": float(gamma),
                "inchi": inchi,
                "smiles": smiles,
                "path": str(Path(path).resolve().relative_to(self.reaction_directory)),
            }
        )
        self.save()

    def complete_cycle(self, gamma, next_orientations):
        """Mark a gamma cycle complete and snapshot candidates for the next cycle."""
        self.data["completed_cycles"].append({"gamma": float(gamma), "status": "completed"})
        self.data["current_cycle"] += 1
        self._replace_pending_orientations(next_orientations)
        self.data["current_survivors"] = []
        self.save()

    def finish(self, status="completed"):
        """Persist terminal workflow state while retaining the run record."""
        self.data["status"] = status
        self.data["pending_orientations"] = []
        self.save()

    def save(self):
        """Atomically write ``reaction/state.json``."""
        self.reaction_directory.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=self.reaction_directory,
            prefix=".state-",
            suffix=".json",
            delete=False,
        ) as fp:
            json.dump(self.data, fp, indent=2, sort_keys=True)
            fp.write("\n")
            temporary_path = Path(fp.name)
        os.replace(temporary_path, self.state_file)

    def _replace_pending_orientations(self, orientations):
        self.data["pending_orientations"] = [
            self._snapshot_molecule(molecule, "pending")
            for molecule in orientations
        ]

    def _replace_current_survivors(self, orientations):
        self.data["current_survivors"] = [
            self._snapshot_molecule(molecule, "survivor")
            for molecule in orientations
        ]

    def _snapshot_molecule(self, molecule, role):
        sequence = self.data["next_snapshot"]
        self.data["next_snapshot"] += 1
        filename = f"{sequence:06d}_{role}_{_safe_name(molecule.name)}.xyz"
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
            "path": str(output_path.relative_to(self.reaction_directory)),
            "name": molecule.name,
            "title": molecule.title,
            "fragments": _json_value(getattr(molecule, "fragments", [])),
            "charge": _json_value(molecule.charge),
            "multiplicity": _json_value(molecule.multiplicity),
            "scftype": molecule.scftype,
            "energy": _json_value(molecule.energy),
        }

    def _load_molecule(self, reference):
        source = self.reaction_directory / reference["path"]
        try:
            loaded = Molecule.from_xyz(str(source))
        except (OSError, ValueError, SystemExit) as exc:
            raise ReactionStateError(
                f"Could not restore reaction geometry snapshot {source}: {exc}"
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


def read_legacy_checkpoint(root_directory):
    """Read a legacy ``jobs.pkl`` only for one-time migration to JSON state."""
    checkpoint_file = Path(root_directory).resolve() / "jobs.pkl"
    if not checkpoint_file.is_file():
        return None
    try:
        with checkpoint_file.open("rb") as fp:
            checkpoint = pickle.load(fp)
    except Exception as exc:
        raise ReactionStateError(
            f"Could not read legacy reaction checkpoint {checkpoint_file}: {exc}"
        ) from exc
    if not isinstance(checkpoint, dict) or not checkpoint:
        raise ReactionStateError(
            f"Legacy reaction checkpoint {checkpoint_file} has no pending state."
        )
    return checkpoint
