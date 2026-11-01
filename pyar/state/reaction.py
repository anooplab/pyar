"""Versioned, inspectable restart state for reaction workflows."""

from __future__ import annotations

import json
import os
import re
import tempfile
from pathlib import Path

import numpy as np

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


class ReactionStateBase:
    """Persist reaction workflow progress independently of optimizer objects."""

    def __init__(self, root_directory, data):
        self.root_directory = Path(root_directory).resolve()
        self.reaction_directory = self.root_directory / "reaction"
        self.state_directory = self.reaction_directory / "state"
        self.geometry_directory = self.state_directory / "geometries"
        self.state_file = self.reaction_directory / STATE_FILENAME
        self.data = data





    def validate_request(self, expected_request):
        """Fail clearly when a resume invocation changes scientific settings."""
        saved = self.data.get("request")
        expected = _json_value(expected_request)
        if saved != expected:
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

    def record_job(self, job_name, gamma, status, remaining_orientations, current_survivors, identity=None, failure=None):
        """Record processed work plus pending and retained current-cycle candidates."""
        record = {"job_name": job_name, "gamma": float(gamma), "status": str(status)}
        if identity:
            record["product_identity"] = _json_value(identity)
        if failure:
            record["failure"] = str(failure)
        self.data["completed_jobs"].append(record)
        self._replace_pending_orientations(remaining_orientations)
        self._replace_current_survivors(current_survivors)
        self.save()


    def complete_cycle(self, gamma, next_orientations):
        """Mark a gamma cycle complete and snapshot candidates for the next cycle."""
        self.data["completed_cycles"].append({"gamma": float(gamma), "status": "completed"})
        self.data["current_cycle"] += 1
        self._replace_pending_orientations(next_orientations)
        self.data["current_survivors"] = []
        self.save()


    def fail(self, error):
        """Persist a failed run without discarding pending orientation evidence."""
        self.data["status"] = "failed_backend"
        self.data["failure"] = str(error)
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
