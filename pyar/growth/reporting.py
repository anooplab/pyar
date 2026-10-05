"""Persistent geometry snapshots for sequential growth."""
import os
import hashlib
import tempfile
from pathlib import Path
from pyar.state.aggregate import _json_value
from pyar.core.molecule import Molecule


def snapshot_pool(molecules, directory, run_directory):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    refs = []
    for index, molecule in enumerate(molecules):
        path = directory / f"selected_{index:03d}.xyz"
        with tempfile.NamedTemporaryFile(dir=directory, suffix=".xyz", delete=False) as fp:
            temporary = fp.name
        try:
            molecule.mol_to_xyz(temporary)
            os.replace(temporary, path)
        finally:
            if os.path.exists(temporary):
                os.unlink(temporary)
        refs.append({"path": str(path.relative_to(run_directory)), "name": molecule.name,
                     "title": molecule.title, "charge": molecule.charge,
                     "multiplicity": molecule.multiplicity, "scftype": molecule.scftype,
                     "fragments": _json_value(molecule.fragments), "energy": molecule.energy,
                     "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                     "growth_metadata": {key: getattr(molecule, key) for key in
                                         ("connectivity_policy_hint", "growth_kind", "source_kind")
                                         if hasattr(molecule, key)}})
    return refs


def restore_pool(references, run_directory):
    pool = []
    for ref in references:
        path = (Path(run_directory) / ref["path"]).resolve()
        if not path.is_relative_to(Path(run_directory).resolve()):
            raise ValueError("Growth snapshot path is outside the run directory")
        if hashlib.sha256(path.read_bytes()).hexdigest() != ref["sha256"]:
            raise ValueError(f"Growth snapshot was modified: {path}")
        loaded = Molecule.from_xyz(str(path))
        pool.append(Molecule(loaded.atoms_list, loaded.coordinates,
                             **{key: ref[key] for key in ("name", "title", "charge", "multiplicity",
                                                         "scftype", "fragments", "energy")}))
        for key, value in ref.get("growth_metadata", {}).items():
            setattr(pool[-1], key, value)
    return pool
