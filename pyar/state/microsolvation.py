"""Versioned, checksummed restart state for solute-centred microsolvation."""

from __future__ import annotations

import hashlib
import json
import os
import tempfile
from pathlib import Path

from pyar.state.grow import write_json


class MicrosolvationStateError(RuntimeError):
    """A microsolvation run cannot be safely resumed."""


class MicrosolvationRunState:
    VERSION = 1

    def __init__(self, directory, data):
        self.directory = Path(directory).resolve()
        self.data = data
        self.path = self.directory / "state.json"

    @classmethod
    def load(cls, directory, request):
        directory = Path(directory).resolve()
        state_path = directory / "state.json"
        if not state_path.exists():
            if directory.exists() and any(directory.iterdir()):
                raise MicrosolvationStateError(
                    f"Non-empty microsolvation directory without valid state: {directory}"
                )
            return None
        try:
            data = json.loads(state_path.read_text(encoding="utf-8"))
        except (OSError, ValueError) as exc:
            raise MicrosolvationStateError(f"Cannot read microsolvation state {state_path}: {exc}") from exc
        if data.get("version") != cls.VERSION or data.get("workflow") != "microsolvation":
            raise MicrosolvationStateError(f"Unsupported microsolvation state format: {state_path}")
        if data.get("request") != request:
            raise MicrosolvationStateError(
                "Microsolvation request differs from saved state; use the original settings or a fresh output directory"
            )
        state = cls(directory, data)
        state.validate_progress()
        return state

    @classmethod
    def create(cls, directory, request, initial_molecule):
        directory = Path(directory).resolve()
        directory.mkdir(parents=True, exist_ok=True)
        (directory / "step_000" / "selected").mkdir(parents=True, exist_ok=True)
        data = {
            "version": cls.VERSION,
            "workflow": "microsolvation",
            "request": request,
            "status": "running",
            "next_step": 1,
            "completed_steps": [],
            "solute_atom_indices": list(range(len(initial_molecule))),
            "solvent_fragments": [],
            "coverage_history": [],
            "current_seeds": [],
        }
        state = cls(directory, data)
        data["current_seeds"] = state.snapshot_pool([initial_molecule], directory / "step_000" / "selected")
        write_json(directory / "request.json", request)
        state.save()
        return state

    def validate_progress(self):
        steps = self.data.get("completed_steps")
        if not isinstance(steps, list):
            raise MicrosolvationStateError("Microsolvation progress is invalid")
        step_ids = [item.get("step") for item in steps]
        if step_ids != list(range(1, len(steps) + 1)):
            raise MicrosolvationStateError("Microsolvation completed steps are not sequential")
        if self.data.get("next_step") != len(steps) + 1:
            raise MicrosolvationStateError("Microsolvation next-step marker is inconsistent")
        if self.data.get("status") not in {"running", "completed", "no_candidates", "stopped"}:
            raise MicrosolvationStateError("Microsolvation status is invalid")
        if steps and [item["path"] for item in self.data.get("current_seeds", [])] != steps[-1].get("selected_paths"):
            raise MicrosolvationStateError("Current pool does not match the last completed microsolvation step")
        for reference in self.data.get("current_seeds", []) + self.data.get("final_seeds", []):
            path = (self.directory / reference["path"]).resolve()
            if not path.is_relative_to(self.directory) or not path.is_file():
                raise MicrosolvationStateError(f"Missing microsolvation snapshot: {path}")
            if _sha256(path) != reference.get("sha256"):
                raise MicrosolvationStateError(f"Microsolvation snapshot was modified: {path}")
        for step in steps:
            for reference in step.get("selected_references", []):
                path = (self.directory / reference["path"]).resolve()
                if not path.is_relative_to(self.directory) or not path.is_file():
                    raise MicrosolvationStateError(f"Missing completed microsolvation snapshot: {path}")
                if _sha256(path) != reference.get("sha256"):
                    raise MicrosolvationStateError(f"Microsolvation snapshot was modified: {path}")

    def snapshot_pool(self, molecules, destination):
        destination = Path(destination)
        destination.mkdir(parents=True, exist_ok=True)
        references = []
        for index, molecule in enumerate(molecules):
            target = destination / f"selected_{index:03d}.xyz"
            with tempfile.NamedTemporaryFile(dir=destination, suffix=".xyz", delete=False) as handle:
                temporary = Path(handle.name)
            try:
                molecule.mol_to_xyz(str(temporary))
                os.replace(temporary, target)
            finally:
                if temporary.exists():
                    temporary.unlink()
            relative = str(target.relative_to(self.directory))
            references.append({
                "path": relative,
                "sha256": _sha256(target),
                "name": molecule.name,
                "title": molecule.title,
                "charge": molecule.charge,
                "multiplicity": molecule.multiplicity,
                "scftype": molecule.scftype,
                "energy": molecule.energy,
                "solute_atom_indices": list(molecule.solute_atom_indices),
                "solvent_fragments": [list(part) for part in molecule.solvent_fragments],
                "coverage": dict(getattr(molecule, "microsolvation_coverage", {})),
                "target_atom": getattr(molecule, "microsolvation_target_atom", None),
            })
        return references

    def restore_pool(self):
        from pyar.core.molecule import Molecule
        result = []
        for reference in self.data.get("current_seeds", []):
            path = (self.directory / reference["path"]).resolve()
            if not path.is_relative_to(self.directory) or _sha256(path) != reference.get("sha256"):
                raise MicrosolvationStateError(f"Microsolvation snapshot is missing or modified: {path}")
            loaded = Molecule.from_xyz(str(path))
            molecule = Molecule(
                loaded.atoms_list, loaded.coordinates, name=reference["name"], title=reference["title"],
                charge=reference["charge"], multiplicity=reference["multiplicity"],
                scftype=reference["scftype"], energy=reference.get("energy"),
            )
            molecule.solute_atom_indices = tuple(reference["solute_atom_indices"])
            molecule.solvent_fragments = tuple(tuple(part) for part in reference["solvent_fragments"])
            molecule.fragments = [list(molecule.solute_atom_indices), *[list(x) for x in molecule.solvent_fragments]]
            molecule.fragments_coordinates = molecule.split_coordinates()
            molecule.fragments_atoms_list = molecule.split_atoms_lists()
            molecule.microsolvation_coverage = dict(reference.get("coverage", {}))
            molecule.microsolvation_target_atom = reference.get("target_atom")
            result.append(molecule)
        return result

    def complete_step(self, step, molecules, diagnostics):
        if step != self.data["next_step"]:
            raise MicrosolvationStateError("Microsolvation step completion is out of sequence")
        references = self.snapshot_pool(molecules, self.directory / f"step_{step:03d}" / "selected")
        coverage_records = [dict(ref["coverage"], name=ref["name"]) for ref in references]
        step_record = {
            "step": step,
            "solvent_count": step,
            "selected_count": len(references),
            "selected_paths": [ref["path"] for ref in references],
            "selected_references": references,
            "coverage": coverage_records,
            "selection": diagnostics,
        }
        write_json(self.directory / f"step_{step:03d}" / "coverage.json", step_record)
        self.data["completed_steps"].append(step_record)
        self.data["coverage_history"].append(coverage_records)
        self.data["current_seeds"] = references
        self.data["next_step"] = step + 1
        self.save()

    def finish(self, status):
        self.data["status"] = status
        self.save()

    def save(self):
        write_json(self.path, self.data)


def _sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()
