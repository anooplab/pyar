"""RDKit-backed conformational-search workflow for PyAR."""

from __future__ import annotations

import csv
import json
import logging
import math
import os
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from pyar import __version__
from pyar.backends import require_executable
from pyar.backends.subprocess_utils import run_command
from pyar.core.molecule import Molecule
from pyar.optional_dependencies import optional_dependency_error
from pyar.optimiser import _status_label, is_usable, optimise
from pyar.workflow_results import ConformerResult
from pyar.workflows._growth import _working_directory, workflow_run_directory, workflow_state_path

conformer_logger = logging.getLogger("pyar.workflows.conformer")

STATE_VERSION = 2


class ConformerWorkflowError(RuntimeError):
    """Raised when conformational search cannot complete safely."""


@dataclass
class ConformerRecord:
    """One generated conformer and its ranking metadata."""

    seed: int
    source_conf_id: int
    rdkit_energy: float
    rdkit_status: str
    force_field: str
    rdkit_molecule: object | None = None
    rdkit_conf_id: int | None = None
    rank: int | None = None
    name: str | None = None
    molecule: Molecule | None = None
    rdkit_path: str | None = None
    generation_stage: str = "etkdg"
    generation_round: int = 0
    parent_name: str | None = None
    torsion_moves: str | None = None
    refinement_reason: str | None = None
    refinement_diversity: float | None = None
    contact_count: int | None = None
    radius_gyration: float | None = None
    backend_status: str | None = None
    backend_energy: float | None = None
    backend_path: str | None = None
    selected_path: str | None = None


def _json_value(value):
    """Return a JSON-serializable representation for workflow state."""
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


def _write_json_atomic(path, data):
    """Atomically write JSON data."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        dir=path.parent,
        prefix=".state-",
        suffix=".json",
        delete=False,
    ) as fp:
        json.dump(_json_value(data), fp, indent=2, sort_keys=True)
        fp.write("\n")
        temporary_path = Path(fp.name)
    os.replace(temporary_path, path)


def _rdkit_modules():
    """Import RDKit lazily and raise PyAR's optional-extra hint if missing."""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError as exc:
        raise optional_dependency_error(
            "rdkit",
            feature="RDKit conformer search",
            extra="conformer",
        ) from exc
    return Chem, AllChem


def _resolve_input_format(input_spec, input_format):
    """Return the concrete input format for one conformer-search input."""
    normalized = str(input_format).lower()
    if normalized != "auto":
        return normalized

    path = Path(input_spec)
    suffix = path.suffix.lower()
    if path.exists() and suffix == ".xyz":
        return "xyz"
    if path.exists() and suffix in {".sdf", ".sd"}:
        return "sdf"
    if path.exists() and suffix == ".mol":
        return "mol"
    return "smiles"


def _first_molecule_from_supplier(supplier, source):
    """Return the first valid RDKit molecule from a supplier."""
    for molecule in supplier:
        if molecule is not None:
            return molecule
    raise ConformerWorkflowError(f"No valid molecule could be read from {source}")


def _load_sdf_molecule(path, Chem):
    supplier = Chem.SDMolSupplier(str(path), removeHs=False)
    return _first_molecule_from_supplier(supplier, path)


def _load_xyz_with_openbabel(input_path, run_directory, Chem):
    """Convert an XYZ file to SDF with OpenBabel and load it into RDKit."""
    source = Path(input_path)
    if not source.is_file():
        raise ConformerWorkflowError(f"XYZ input file {input_path!r} does not exist")
    executable = require_executable("obabel", "OpenBabel")
    input_directory = run_directory / "input"
    input_directory.mkdir(parents=True, exist_ok=True)
    converted_sdf = input_directory / "input_from_xyz.sdf"
    log_path = input_directory / "openbabel_xyz_to_sdf.log"
    exit_status = run_command(
        [executable, "-ixyz", source, "-osdf", "-O", converted_sdf],
        stdout_path=log_path,
        stderr_path=log_path,
    )
    if exit_status != 0 or not converted_sdf.is_file():
        raise ConformerWorkflowError(
            "OpenBabel could not infer an RDKit-readable bond graph from the XYZ input. "
            f"See {log_path}."
        )
    return _load_sdf_molecule(converted_sdf, Chem)


def _load_rdkit_molecule(input_spec, input_format, Chem, run_directory):
    """Load the requested molecule into RDKit and add explicit hydrogens."""
    resolved_format = _resolve_input_format(input_spec, input_format)
    if resolved_format == "smiles":
        molecule = Chem.MolFromSmiles(str(input_spec))
    elif resolved_format == "sdf":
        molecule = _load_sdf_molecule(Path(input_spec), Chem)
    elif resolved_format == "mol":
        molecule = Chem.MolFromMolFile(str(input_spec), removeHs=False)
    elif resolved_format == "xyz":
        molecule = _load_xyz_with_openbabel(input_spec, run_directory, Chem)
    else:
        raise ConformerWorkflowError(f"Unsupported conformer input format: {input_format!r}")

    if molecule is None:
        raise ConformerWorkflowError(f"RDKit could not read input {input_spec!r}")
    if hasattr(Chem, "SanitizeMol"):
        Chem.SanitizeMol(molecule)
    molecule = Chem.AddHs(molecule, addCoords=True)
    return molecule, resolved_format


def _force_field_for_molecule(rdkit_molecule, AllChem, force_field):
    """Return the concrete RDKit force field to use."""
    requested = str(force_field).lower()
    has_mmff = bool(AllChem.MMFFHasAllMoleculeParams(rdkit_molecule))
    has_uff = bool(AllChem.UFFHasAllMoleculeParams(rdkit_molecule))
    if requested == "auto":
        if has_mmff:
            return "mmff"
        if has_uff:
            return "uff"
        raise ConformerWorkflowError("RDKit has neither MMFF nor UFF parameters for this molecule")
    if requested == "mmff" and not has_mmff:
        raise ConformerWorkflowError("RDKit MMFF parameters are unavailable for this molecule")
    if requested == "uff" and not has_uff:
        raise ConformerWorkflowError("RDKit UFF parameters are unavailable for this molecule")
    if requested not in {"mmff", "uff"}:
        raise ConformerWorkflowError(f"Unsupported RDKit force field: {force_field!r}")
    return requested


def _embed_parameters(AllChem, seed, prune_rms_threshold, num_threads, use_random_coords):
    """Build RDKit ETKDG parameters with deterministic defaults."""
    if hasattr(AllChem, "ETKDGv3"):
        parameters = AllChem.ETKDGv3()
    else:
        parameters = AllChem.ETKDG()
    parameters.randomSeed = int(seed)
    parameters.pruneRmsThresh = float(prune_rms_threshold)
    parameters.numThreads = int(num_threads)
    parameters.useRandomCoords = bool(use_random_coords)
    return parameters


def _embed_and_rank_conformers(
    rdkit_molecule,
    AllChem,
    *,
    seed,
    num_conformers,
    prune_rms_threshold,
    force_field,
    num_threads,
    use_random_coords,
    max_iterations,
):
    """Generate, minimize, and rank RDKit conformers."""
    parameters = _embed_parameters(
        AllChem,
        seed,
        prune_rms_threshold,
        num_threads,
        use_random_coords,
    )
    conformer_ids = list(
        AllChem.EmbedMultipleConfs(
            rdkit_molecule,
            numConfs=int(num_conformers),
            params=parameters,
        )
    )
    if not conformer_ids:
        return []

    selected_force_field = _force_field_for_molecule(rdkit_molecule, AllChem, force_field)
    if selected_force_field == "mmff":
        results = AllChem.MMFFOptimizeMoleculeConfs(
            rdkit_molecule,
            numThreads=int(num_threads),
            maxIters=int(max_iterations),
        )
    else:
        results = AllChem.UFFOptimizeMoleculeConfs(
            rdkit_molecule,
            numThreads=int(num_threads),
            maxIters=int(max_iterations),
        )

    records = []
    for conformer_id, result in zip(conformer_ids, results):
        not_converged, energy = result
        try:
            energy = float(energy)
        except (TypeError, ValueError):
            continue
        if not math.isfinite(energy):
            continue
        records.append(
            ConformerRecord(
                seed=int(seed),
                source_conf_id=int(conformer_id),
                rdkit_energy=energy,
                rdkit_status="not_converged" if int(not_converged) else "converged",
                force_field=selected_force_field,
                rdkit_molecule=rdkit_molecule,
                rdkit_conf_id=int(conformer_id),
                generation_stage="etkdg",
                generation_round=0,
            )
        )
    return sorted(records, key=lambda record: (record.rdkit_energy, record.source_conf_id))


def _parent_record_name(record):
    """Return the stable name assigned to an embedded RDKit conformer."""
    if record.name:
        return record.name
    return f"seed_{record.seed:06d}_conf_{record.source_conf_id:04d}"


def _record_conformer_id(record):
    """Return the RDKit conformer identifier stored in a workflow record."""
    if record.rdkit_conf_id is not None:
        return int(record.rdkit_conf_id)
    return int(record.source_conf_id)


def _rotatable_torsions(rdkit_molecule, Chem):
    """Return rotatable-bond torsion quartets for torsion-kick sampling."""
    pattern = Chem.MolFromSmarts("[!$(*#*)&!D1]-!@[!$(*#*)&!D1]")
    if pattern is None:
        return []
    torsions = []
    seen_bonds = set()
    for begin_index, end_index in rdkit_molecule.GetSubstructMatches(pattern):
        bond = rdkit_molecule.GetBondBetweenAtoms(int(begin_index), int(end_index))
        if bond is None:
            continue
        bond_key = tuple(sorted((int(begin_index), int(end_index))))
        if bond_key in seen_bonds:
            continue
        seen_bonds.add(bond_key)
        begin_atom = rdkit_molecule.GetAtomWithIdx(int(begin_index))
        end_atom = rdkit_molecule.GetAtomWithIdx(int(end_index))
        begin_neighbors = [
            atom.GetIdx()
            for atom in begin_atom.GetNeighbors()
            if atom.GetIdx() != int(end_index) and atom.GetAtomicNum() > 1
        ]
        end_neighbors = [
            atom.GetIdx()
            for atom in end_atom.GetNeighbors()
            if atom.GetIdx() != int(begin_index) and atom.GetAtomicNum() > 1
        ]
        if not begin_neighbors or not end_neighbors:
            continue
        torsions.append((int(begin_neighbors[0]), int(begin_index), int(end_index), int(end_neighbors[0])))
    return torsions


def _torsion_importance_scores(record, Chem):
    """Return torsion scores biased toward contacts that can be folded shut."""
    rdkit_molecule = record.rdkit_molecule
    if rdkit_molecule is None:
        return [], []
    torsions = _rotatable_torsions(rdkit_molecule, Chem)
    if not torsions:
        return [], []
    try:
        conformer = rdkit_molecule.GetConformer(_record_conformer_id(record))
    except (AttributeError, RuntimeError, ValueError):
        try:
            conformer = rdkit_molecule.GetConformer()
        except (AttributeError, RuntimeError, ValueError):
            return torsions, [0.0 for _ in torsions]

    try:
        topology = np.asarray(Chem.GetDistanceMatrix(rdkit_molecule), dtype=float)
    except Exception:
        topology = None

    coordinates = []
    heavy_indices = []
    for atom in rdkit_molecule.GetAtoms():
        atom_index = atom.GetIdx()
        position = conformer.GetAtomPosition(atom_index)
        coordinates.append([float(position.x), float(position.y), float(position.z)])
        if atom.GetAtomicNum() > 1:
            heavy_indices.append(atom_index)
    if not coordinates:
        return torsions, [0.0 for _ in torsions]

    coordinates = np.asarray(coordinates, dtype=float)
    heavy_indices = heavy_indices or list(range(len(coordinates)))
    scores = [0.0 for _ in torsions]
    torsion_bonds = [tuple(sorted((torsion[1], torsion[2]))) for torsion in torsions]

    for left_pos, left_index in enumerate(heavy_indices[:-1]):
        for right_index in heavy_indices[left_pos + 1 :]:
            if topology is not None and topology[left_index, right_index] <= 3:
                continue
            distance = float(np.linalg.norm(coordinates[left_index] - coordinates[right_index]))
            if distance > 3.5:
                continue
            if distance <= 0.0:
                continue
            try:
                path = Chem.GetShortestPath(rdkit_molecule, int(left_index), int(right_index))
            except Exception:
                path = ()
            if len(path) < 2:
                continue
            pair_weight = (3.5 - distance) + 0.1 * max(0.0, float(len(path)) - 2.0)
            path_bonds = {tuple(sorted((int(path[i]), int(path[i + 1])))) for i in range(len(path) - 1)}
            for torsion_index, bond in enumerate(torsion_bonds):
                if bond in path_bonds:
                    scores[torsion_index] += pair_weight

    return torsions, scores


def _rank_torsions_by_fold(record, Chem):
    """Return torsions ordered from most fold-relevant to least relevant."""
    torsions, scores = _torsion_importance_scores(record, Chem)
    if not torsions:
        return [], []
    if not any(score > 0.0 for score in scores):
        return torsions, list(range(len(torsions)))
    order = sorted(range(len(torsions)), key=lambda index: (-scores[index], index))
    return torsions, order


def _rdkit_minimize_conformer(rdkit_molecule, AllChem, force_field, conf_id, max_iterations):
    """Minimize one RDKit conformer and return status plus energy."""
    if force_field == "mmff":
        properties = AllChem.MMFFGetMoleculeProperties(rdkit_molecule)
        if properties is None:
            return None, None
        status = AllChem.MMFFOptimizeMolecule(
            rdkit_molecule,
            mmffVariant="MMFF94",
            confId=int(conf_id),
            maxIters=int(max_iterations),
        )
        field = AllChem.MMFFGetMoleculeForceField(
            rdkit_molecule,
            properties,
            confId=int(conf_id),
        )
    else:
        status = AllChem.UFFOptimizeMolecule(
            rdkit_molecule,
            confId=int(conf_id),
            maxIters=int(max_iterations),
        )
        field = AllChem.UFFGetMoleculeForceField(rdkit_molecule, confId=int(conf_id))
    if field is None:
        return None, None
    energy = float(field.CalcEnergy())
    if not math.isfinite(energy):
        return None, None
    return int(status), energy


def _apply_torsion_moves(work_molecule, conf_id, torsions, angle_choices, move_indices, angle_indices):
    """Apply torsion moves in-place and return a compact move description."""
    moves = []
    conformer = work_molecule.GetConformer(conf_id)
    from rdkit.Chem import rdMolTransforms

    for torsion_index, angle_index in zip(move_indices, angle_indices):
        atom_i, atom_j, atom_k, atom_l = torsions[int(torsion_index)]
        delta = float(angle_choices[int(angle_index)])
        current = float(rdMolTransforms.GetDihedralDeg(conformer, atom_i, atom_j, atom_k, atom_l))
        rdMolTransforms.SetDihedralDeg(conformer, atom_i, atom_j, atom_k, atom_l, current + delta)
        moves.append(f"{atom_i}-{atom_j}-{atom_k}-{atom_l}:{delta:.1f}")
    return moves


def _torsion_move_indices(num_torsions, *, parent_index, round_index, kick_index, max_bonds):
    """Return a deterministic torsion subset for grid-style exploration."""
    move_count = min(int(max_bonds), int(num_torsions))
    if move_count < 1:
        return []
    start = (int(parent_index) + int(round_index) + int(kick_index)) % int(num_torsions)
    return [int((start + 2 * offset) % int(num_torsions)) for offset in range(move_count)]


def _generate_torsion_random_records(
    records,
    Chem,
    AllChem,
    *,
    enabled,
    kicks_per_conformer,
    max_bonds,
    seed,
    round_index,
    max_iterations,
):
    """Generate locally minimized torsion-kicked conformers from RDKit records."""
    if not enabled or int(kicks_per_conformer) < 1 or int(max_bonds) < 1:
        return []
    angle_choices = np.asarray([60.0, 120.0, 180.0, 240.0, 300.0], dtype=float)
    kicked_records = []
    next_ids_by_seed = {}
    for parent_index, parent in enumerate(sorted(records, key=_record_sort_key)):
        rdkit_molecule = parent.rdkit_molecule
        if rdkit_molecule is None:
            continue
        torsions = _rotatable_torsions(rdkit_molecule, Chem)
        if not torsions:
            continue
        next_ids_by_seed.setdefault(parent.seed, max(record.source_conf_id for record in records if record.seed == parent.seed) + 1)
        for kick_index in range(int(kicks_per_conformer)):
            try:
                work_molecule = Chem.Mol(rdkit_molecule)
                parent_conformer = rdkit_molecule.GetConformer(_record_conformer_id(parent))
                kicked_conformer = Chem.Conformer(parent_conformer)
                if hasattr(work_molecule, "RemoveAllConformers"):
                    work_molecule.RemoveAllConformers()
                conf_id = int(work_molecule.AddConformer(kicked_conformer, assignId=True))
            except (AttributeError, RuntimeError, ValueError):
                continue

            rng_seed = (
                int(seed) * 1_000_003
                + int(round_index) * 100_003
                + int(parent.seed) * 10_007
                + int(parent.source_conf_id) * 101
                + int(kick_index)
                + int(parent_index)
            )
            rng = np.random.default_rng(rng_seed)
            move_count = min(int(max_bonds), len(torsions))
            torsion_indices = rng.choice(len(torsions), size=move_count, replace=False)
            angle_indices = rng.choice(len(angle_choices), size=move_count, replace=True)
            try:
                moves = _apply_torsion_moves(
                    work_molecule,
                    conf_id,
                    torsions,
                    angle_choices,
                    np.atleast_1d(torsion_indices),
                    np.atleast_1d(angle_indices),
                )
            except Exception:
                continue

            status, energy = _rdkit_minimize_conformer(
                work_molecule,
                AllChem,
                parent.force_field,
                conf_id,
                max_iterations,
            )
            if status is None or energy is None:
                continue
            source_conf_id = next_ids_by_seed[parent.seed]
            next_ids_by_seed[parent.seed] += 1
            kicked_records.append(
                ConformerRecord(
                    seed=parent.seed,
                    source_conf_id=int(source_conf_id),
                    rdkit_energy=float(energy),
                    rdkit_status="not_converged" if status else "converged",
                    force_field=parent.force_field,
                    rdkit_molecule=work_molecule,
                    rdkit_conf_id=int(conf_id),
                    generation_stage="torsion_kick",
                    generation_round=int(round_index),
                    parent_name=_parent_record_name(parent),
                    torsion_moves=";".join(moves),
                )
            )
    return sorted(kicked_records, key=_record_sort_key)


def _generate_torsion_kick_records(*args, **kwargs):
    """Backward-compatible alias for the stochastic torsion generator."""
    return _generate_torsion_random_records(*args, **kwargs)


def _generate_torsion_guided_records(
    records,
    Chem,
    AllChem,
    *,
    enabled,
    kicks_per_conformer,
    max_bonds,
    round_index,
    max_iterations,
):
    """Generate deterministic torsion trials guided by fold-contact scores."""
    if not enabled or int(kicks_per_conformer) < 1 or int(max_bonds) < 1:
        return []
    angle_choices = np.asarray([-180.0, -120.0, -60.0, 60.0, 120.0, 180.0], dtype=float)
    kicked_records = []
    next_ids_by_seed = {}
    for parent_index, parent in enumerate(sorted(records, key=_record_sort_key)):
        rdkit_molecule = parent.rdkit_molecule
        if rdkit_molecule is None:
            continue
        torsions, torsion_order = _rank_torsions_by_fold(parent, Chem)
        if not torsions:
            continue
        if not torsion_order:
            torsion_order = list(range(len(torsions)))
        next_ids_by_seed.setdefault(parent.seed, max(record.source_conf_id for record in records if record.seed == parent.seed) + 1)
        move_count = min(int(max_bonds), len(torsion_order))
        if move_count < 1:
            continue
        for kick_index in range(int(kicks_per_conformer)):
            try:
                work_molecule = Chem.Mol(rdkit_molecule)
                parent_conformer = rdkit_molecule.GetConformer(_record_conformer_id(parent))
                kicked_conformer = Chem.Conformer(parent_conformer)
                if hasattr(work_molecule, "RemoveAllConformers"):
                    work_molecule.RemoveAllConformers()
                conf_id = int(work_molecule.AddConformer(kicked_conformer, assignId=True))
            except (AttributeError, RuntimeError, ValueError):
                continue

            start = (int(parent_index) + int(round_index) + int(kick_index)) % len(torsion_order)
            move_indices = [int(torsion_order[(start + offset) % len(torsion_order)]) for offset in range(move_count)]
            angle_indices = [
                (int(parent_index) + int(round_index) + int(kick_index) + offset) % len(angle_choices)
                for offset in range(move_count)
            ]
            try:
                moves = _apply_torsion_moves(
                    work_molecule,
                    conf_id,
                    torsions,
                    angle_choices,
                    move_indices,
                    angle_indices,
                )
            except Exception:
                continue

            status, energy = _rdkit_minimize_conformer(
                work_molecule,
                AllChem,
                parent.force_field,
                conf_id,
                max_iterations,
            )
            if status is None or energy is None:
                continue
            source_conf_id = next_ids_by_seed[parent.seed]
            next_ids_by_seed[parent.seed] += 1
            kicked_records.append(
                ConformerRecord(
                    seed=parent.seed,
                    source_conf_id=int(source_conf_id),
                    rdkit_energy=float(energy),
                    rdkit_status="not_converged" if status else "converged",
                    force_field=parent.force_field,
                    rdkit_molecule=work_molecule,
                    rdkit_conf_id=int(conf_id),
                    generation_stage="torsion_evolve",
                    generation_round=int(round_index),
                    parent_name=_parent_record_name(parent),
                    torsion_moves=";".join(moves),
                )
            )
    return sorted(kicked_records, key=_record_sort_key)


def _generate_torsion_grid_records(
    records,
    Chem,
    AllChem,
    *,
    enabled,
    kicks_per_conformer,
    max_bonds,
    round_index,
    max_iterations,
):
    """Generate deterministic torsion-grid conformers from RDKit records."""
    if not enabled or int(kicks_per_conformer) < 1 or int(max_bonds) < 1:
        return []
    angle_choices = np.asarray([60.0, 120.0, 180.0, 240.0, 300.0], dtype=float)
    kicked_records = []
    next_ids_by_seed = {}
    for parent_index, parent in enumerate(sorted(records, key=_record_sort_key)):
        rdkit_molecule = parent.rdkit_molecule
        if rdkit_molecule is None:
            continue
        torsions = _rotatable_torsions(rdkit_molecule, Chem)
        if not torsions:
            continue
        next_ids_by_seed.setdefault(parent.seed, max(record.source_conf_id for record in records if record.seed == parent.seed) + 1)
        for kick_index in range(int(kicks_per_conformer)):
            try:
                work_molecule = Chem.Mol(rdkit_molecule)
                parent_conformer = rdkit_molecule.GetConformer(_record_conformer_id(parent))
                kicked_conformer = Chem.Conformer(parent_conformer)
                if hasattr(work_molecule, "RemoveAllConformers"):
                    work_molecule.RemoveAllConformers()
                conf_id = int(work_molecule.AddConformer(kicked_conformer, assignId=True))
            except (AttributeError, RuntimeError, ValueError):
                continue

            move_indices = _torsion_move_indices(
                len(torsions),
                parent_index=parent_index,
                round_index=round_index,
                kick_index=kick_index,
                max_bonds=max_bonds,
            )
            if not move_indices:
                continue
            angle_indices = [
                (parent_index + round_index + kick_index + offset) % len(angle_choices)
                for offset in range(len(move_indices))
            ]
            try:
                moves = _apply_torsion_moves(
                    work_molecule,
                    conf_id,
                    torsions,
                    angle_choices,
                    move_indices,
                    angle_indices,
                )
            except Exception:
                continue

            status, energy = _rdkit_minimize_conformer(
                work_molecule,
                AllChem,
                parent.force_field,
                conf_id,
                max_iterations,
            )
            if status is None or energy is None:
                continue
            source_conf_id = next_ids_by_seed[parent.seed]
            next_ids_by_seed[parent.seed] += 1
            kicked_records.append(
                ConformerRecord(
                    seed=parent.seed,
                    source_conf_id=int(source_conf_id),
                    rdkit_energy=float(energy),
                    rdkit_status="not_converged" if status else "converged",
                    force_field=parent.force_field,
                    rdkit_molecule=work_molecule,
                    rdkit_conf_id=int(conf_id),
                    generation_stage="torsion_kick",
                    generation_round=int(round_index),
                    parent_name=_parent_record_name(parent),
                    torsion_moves=";".join(moves),
                )
            )
    return sorted(kicked_records, key=_record_sort_key)


def _generate_torsion_mc_records(
    records,
    Chem,
    AllChem,
    *,
    enabled,
    kicks_per_conformer,
    max_bonds,
    seed,
    round_index,
    max_iterations,
    temperature=2.0,
    tabu_rms=0.5,
):
    """Generate Monte Carlo torsion-walk conformers with tabu-style rejection."""
    if not enabled or int(kicks_per_conformer) < 1 or int(max_bonds) < 1:
        return []
    angle_choices = np.asarray([60.0, 120.0, 180.0, 240.0, 300.0], dtype=float)
    kicked_records = []
    tabu_records = list(records)
    next_ids_by_seed = {}
    for parent_index, parent in enumerate(sorted(records, key=_record_sort_key)):
        rdkit_molecule = parent.rdkit_molecule
        if rdkit_molecule is None:
            continue
        torsions = _rotatable_torsions(rdkit_molecule, Chem)
        if not torsions:
            continue
        next_ids_by_seed.setdefault(parent.seed, max(record.source_conf_id for record in records if record.seed == parent.seed) + 1)
        current_record = parent
        current_energy = float(parent.rdkit_energy)
        for step_index in range(int(kicks_per_conformer)):
            try:
                work_molecule = Chem.Mol(current_record.rdkit_molecule)
                parent_conformer = current_record.rdkit_molecule.GetConformer(_record_conformer_id(current_record))
                kicked_conformer = Chem.Conformer(parent_conformer)
                if hasattr(work_molecule, "RemoveAllConformers"):
                    work_molecule.RemoveAllConformers()
                conf_id = int(work_molecule.AddConformer(kicked_conformer, assignId=True))
            except (AttributeError, RuntimeError, ValueError):
                continue

            rng_seed = (
                int(seed) * 1_000_003
                + int(round_index) * 100_003
                + int(parent.seed) * 10_007
                + int(current_record.source_conf_id) * 101
                + int(step_index)
                + int(parent_index)
            )
            rng = np.random.default_rng(rng_seed)
            move_count = min(int(max_bonds), len(torsions))
            torsion_indices = rng.choice(len(torsions), size=move_count, replace=False)
            angle_indices = rng.choice(len(angle_choices), size=move_count, replace=True)
            try:
                moves = _apply_torsion_moves(
                    work_molecule,
                    conf_id,
                    torsions,
                    angle_choices,
                    np.atleast_1d(torsion_indices),
                    np.atleast_1d(angle_indices),
                )
            except Exception:
                continue

            status, energy = _rdkit_minimize_conformer(
                work_molecule,
                AllChem,
                parent.force_field,
                conf_id,
                max_iterations,
            )
            if status is None or energy is None:
                continue

            candidate = ConformerRecord(
                seed=parent.seed,
                source_conf_id=int(next_ids_by_seed[parent.seed]),
                rdkit_energy=float(energy),
                rdkit_status="not_converged" if status else "converged",
                force_field=parent.force_field,
                rdkit_molecule=work_molecule,
                rdkit_conf_id=int(conf_id),
                generation_stage="torsion_mc",
                generation_round=int(round_index),
                parent_name=_parent_record_name(parent),
                torsion_moves=";".join(moves),
            )

            if _is_duplicate_record(candidate, tabu_records, tabu_rms):
                continue
            next_ids_by_seed[parent.seed] += 1
            tabu_records.append(candidate)

            delta = float(energy - current_energy)
            accept = delta <= 0.0
            if not accept:
                try:
                    accept = float(rng.random()) < math.exp(-delta / float(temperature))
                except OverflowError:
                    accept = False
            if accept:
                candidate.generation_stage = "torsion_mc"
                kicked_records.append(candidate)
                current_record = candidate
                current_energy = float(energy)
    return sorted(kicked_records, key=_record_sort_key)


def _generate_torsion_evolution_records(
    records,
    Chem,
    AllChem,
    *,
    enabled,
    generations,
    kicks_per_conformer,
    max_bonds,
    max_iterations,
    diversity_fraction=0.2,
    compactness_fraction=0.2,
    parent_limit=None,
    tabu_rms=0.5,
):
    """Generate population-based torsion mutants with elitist selection."""
    if not enabled or int(generations) < 1 or int(kicks_per_conformer) < 1 or int(max_bonds) < 1:
        return []
    population = sorted(records, key=_record_sort_key)
    if not population:
        return []
    limit = len(population) if parent_limit is None else int(parent_limit)
    limit = max(1, min(limit, len(population)))
    offspring_records = []
    tabu_records = list(population)
    for generation_index in range(1, int(generations) + 1):
        parents = _select_torsion_parent_records(
            population,
            limit,
            diversity_fraction,
            compactness_fraction,
        )
        if not parents:
            break
        generation_records = _generate_torsion_guided_records(
            parents,
            Chem,
            AllChem,
            enabled=True,
            kicks_per_conformer=kicks_per_conformer,
            max_bonds=max_bonds,
            round_index=generation_index,
            max_iterations=max_iterations,
        )
        if not generation_records:
            break
        novel_records = []
        for candidate in generation_records:
            if _is_duplicate_record(candidate, tabu_records, tabu_rms):
                continue
            candidate.generation_round = int(generation_index)
            novel_records.append(candidate)
            tabu_records.append(candidate)
        if not novel_records:
            break
        offspring_records.extend(novel_records)
        population = _deduplicate_records(population + novel_records, tabu_rms)
        population = _select_torsion_parent_records(
            population,
            limit,
            diversity_fraction,
            compactness_fraction,
        )
    return sorted(offspring_records, key=_record_sort_key)


def _formal_charge(rdkit_molecule):
    total = 0
    for atom in rdkit_molecule.GetAtoms():
        if hasattr(atom, "GetFormalCharge"):
            total += int(atom.GetFormalCharge())
    return total


def _molecule_from_rdkit_conformer(
    rdkit_molecule,
    conformer_id,
    *,
    name,
    seed,
    charge,
    multiplicity,
    scftype,
    energy,
):
    """Convert one RDKit conformer to a PyAR Molecule."""
    conformer = rdkit_molecule.GetConformer(int(conformer_id))
    atoms = [atom.GetSymbol() for atom in rdkit_molecule.GetAtoms()]
    coordinates = []
    for atom_index in range(len(atoms)):
        position = conformer.GetAtomPosition(atom_index)
        coordinates.append([float(position.x), float(position.y), float(position.z)])
    return Molecule(
        atoms,
        coordinates,
        name=name,
        title=f"RDKit seed {seed} conformer {conformer_id}",
        charge=_formal_charge(rdkit_molecule) if charge is None else int(charge),
        multiplicity=int(multiplicity),
        scftype=scftype,
        energy=float(energy),
    )


def _prepare_run_directory(root_directory):
    """Create the conformer run directory, rejecting completed prior runs."""
    run_directory = Path(workflow_run_directory(root_directory, "conformers"))
    state_file = Path(workflow_state_path(str(run_directory)))
    if state_file.exists():
        try:
            with state_file.open(encoding="utf-8") as fp:
                existing_state = json.load(fp)
        except (OSError, ValueError):
            existing_state = {"status": "unknown"}
        raise ConformerWorkflowError(
            f"Conformer state in {state_file} is already "
            f"{existing_state.get('status', 'unknown')!r}; start a new calculation "
            "in a new directory."
        )
    for child in ("input", "rdkit", "selected"):
        (run_directory / child).mkdir(parents=True, exist_ok=True)
    return run_directory


def _write_rdkit_conformers(
    records,
    rdkit_molecule,
    rdkit_directory,
    *,
    charge,
    multiplicity,
    scftype,
):
    """Write ranked RDKit conformers and attach PyAR molecules to records."""
    for rank, record in enumerate(records):
        if record.generation_stage == "torsion_kick":
            name = f"seed_{record.seed:06d}_kick_{record.source_conf_id:04d}"
        else:
            name = _parent_record_name(record)
        source_molecule = record.rdkit_molecule if record.rdkit_molecule is not None else rdkit_molecule
        conformer_id = record.source_conf_id if record.rdkit_conf_id is None else record.rdkit_conf_id
        molecule = _molecule_from_rdkit_conformer(
            source_molecule,
            conformer_id,
            name=name,
            seed=record.seed,
            charge=charge,
            multiplicity=multiplicity,
            scftype=scftype,
            energy=record.rdkit_energy,
        )
        path = rdkit_directory / f"{name}.xyz"
        molecule.mol_to_xyz(str(path))
        record.rank = rank
        record.name = name
        record.molecule = molecule
        record.rdkit_path = str(path)


def _record_sort_key(record):
    return (record.rdkit_energy, record.seed, record.source_conf_id)


def _diversity_coordinates(record):
    """Return heavy-atom coordinates for conformer diversity scoring."""
    molecule = record.molecule
    if molecule is not None and molecule.coordinates is not None:
        coordinates = np.asarray(molecule.coordinates, dtype=float)
        if coordinates.ndim != 2 or coordinates.shape[1] != 3:
            return None
        heavy_indices = [
            index
            for index, atom in enumerate(molecule.atoms_list)
            if str(atom).upper() != "H"
        ]
        if not heavy_indices:
            heavy_indices = list(range(len(molecule.atoms_list)))
        return coordinates[np.asarray(heavy_indices, dtype=int)]

    rdkit_molecule = record.rdkit_molecule
    if rdkit_molecule is None:
        return None
    try:
        conformer_id = record.source_conf_id if record.rdkit_conf_id is None else record.rdkit_conf_id
        conformer = rdkit_molecule.GetConformer(int(conformer_id))
    except (RuntimeError, ValueError):
        try:
            conformer = rdkit_molecule.GetConformer()
        except (RuntimeError, ValueError):
            return None
    coordinates = []
    heavy_indices = []
    for atom in rdkit_molecule.GetAtoms():
        atom_index = atom.GetIdx()
        position = conformer.GetAtomPosition(atom_index)
        coordinates.append([float(position.x), float(position.y), float(position.z)])
        if atom.GetAtomicNum() > 1:
            heavy_indices.append(atom_index)
    if not heavy_indices:
        heavy_indices = list(range(len(coordinates)))
    if not coordinates:
        return None
    return np.asarray(coordinates, dtype=float)[np.asarray(heavy_indices, dtype=int)]


def _deduplicate_records(records, rms_threshold, *, sort_key=_record_sort_key):
    """Return low-energy conformers after heavy-atom RMSD deduplication."""
    threshold = float(rms_threshold)
    ordered_records = list(records) if sort_key is None else sorted(records, key=sort_key)
    if threshold <= 0.0:
        return ordered_records
    selected = []
    for record in ordered_records:
        record_coordinates = _diversity_coordinates(record)
        if record_coordinates is None or len(record_coordinates) < 3:
            selected.append(record)
            continue
        is_duplicate = False
        for accepted in selected:
            accepted_coordinates = _diversity_coordinates(accepted)
            if accepted_coordinates is None or len(accepted_coordinates) < 3:
                continue
            if _kabsch_rmsd(record_coordinates, accepted_coordinates) < threshold:
                is_duplicate = True
                break
        if not is_duplicate:
            selected.append(record)
    return selected


def _select_unique_records(records, limit, rms_threshold):
    """Return up to limit low-energy records after geometry deduplication."""
    return _deduplicate_records(records, rms_threshold, sort_key=None)[: int(limit)]


def _is_duplicate_record(record, references, rms_threshold):
    """Return True when record is within RMS threshold of any reference record."""
    threshold = float(rms_threshold)
    if threshold <= 0.0:
        return False
    record_coordinates = _diversity_coordinates(record)
    if record_coordinates is None or len(record_coordinates) < 3:
        return False
    for reference in references:
        reference_coordinates = _diversity_coordinates(reference)
        if reference_coordinates is None or len(reference_coordinates) < 3:
            continue
        if _kabsch_rmsd(record_coordinates, reference_coordinates) < threshold:
            return True
    return False


def _select_torsion_parent_records(records, limit, diversity_fraction, compactness_fraction):
    """Return a bounded low-energy and compact/diverse pool for another torsion round."""
    records = sorted(records, key=_record_sort_key)
    if int(limit) >= len(records):
        return records
    return _select_refinement_records(
        records,
        int(limit),
        min(int(limit), max(1, int(limit) // 2)),
        diversity_fraction,
        compactness_fraction,
    )


def _compactness_metrics(record):
    """Return intramolecular contact and compactness scores for a conformer."""
    molecule = record.molecule
    if molecule is None or molecule.coordinates is None:
        return 0, float("inf")

    coordinates = np.asarray(molecule.coordinates, dtype=float)
    if coordinates.ndim != 2 or coordinates.shape[1] != 3:
        return 0, float("inf")

    heavy_indices = [
        index
        for index, atom in enumerate(molecule.atoms_list)
        if str(atom).upper() != "H"
    ]
    if len(heavy_indices) < 2:
        heavy_indices = list(range(len(molecule.atoms_list)))
    if len(heavy_indices) < 2:
        return 0, float("inf")

    selected = coordinates[np.asarray(heavy_indices, dtype=int)]
    center = np.mean(selected, axis=0)
    rgyr = float(np.sqrt(np.mean(np.sum((selected - center) ** 2, axis=1))))

    contact_count = 0
    rdkit_molecule = record.rdkit_molecule
    if rdkit_molecule is not None:
        try:
            from rdkit import Chem

            topology = np.asarray(Chem.GetDistanceMatrix(rdkit_molecule), dtype=float)
            for left in range(len(heavy_indices) - 1):
                left_index = heavy_indices[left]
                for right in range(left + 1, len(heavy_indices)):
                    right_index = heavy_indices[right]
                    if topology[left_index, right_index] <= 2:
                        continue
                    if np.linalg.norm(coordinates[left_index] - coordinates[right_index]) <= 3.5:
                        contact_count += 1
        except Exception:
            contact_count = 0

    record.contact_count = int(contact_count)
    record.radius_gyration = float(rgyr)
    return contact_count, rgyr


def _fold_score(record):
    """Return a tuple that prefers contact-rich, compact conformers."""
    contact_count, rgyr = _compactness_metrics(record)
    return contact_count, -rgyr


def _open_score(record):
    """Return a tuple that prefers extended conformers for basin coverage."""
    contact_count, rgyr = _compactness_metrics(record)
    return rgyr, -contact_count


def _kabsch_rmsd(reference, mobile):
    """Return RMSD after optimal translation and proper rotation."""
    reference = np.asarray(reference, dtype=float)
    mobile = np.asarray(mobile, dtype=float)
    if reference.shape != mobile.shape:
        return float("inf")
    reference_centered = reference - np.mean(reference, axis=0)
    mobile_centered = mobile - np.mean(mobile, axis=0)
    covariance = mobile_centered.T @ reference_centered
    left_vectors, _, right_vectors = np.linalg.svd(covariance)
    correction = np.eye(3)
    if np.linalg.det(left_vectors @ right_vectors) < 0:
        correction[-1, -1] = -1.0
    rotation = left_vectors @ correction @ right_vectors
    aligned = mobile_centered @ rotation
    difference = aligned - reference_centered
    return float(np.sqrt(np.mean(np.sum(difference * difference, axis=1))))


def _conformer_rmsd(left, right):
    """Return a heavy-atom RMSD between two conformer records."""
    left_coordinates = _diversity_coordinates(left)
    right_coordinates = _diversity_coordinates(right)
    if left_coordinates is None or right_coordinates is None:
        return 0.0
    return _kabsch_rmsd(left_coordinates, right_coordinates)


def _append_selected_record(selected, selected_ids, record, reason, diversity=None):
    """Append record if not selected and annotate why it was retained."""
    if record is None or id(record) in selected_ids:
        return False
    record.refinement_reason = reason
    record.refinement_diversity = diversity
    _compactness_metrics(record)
    selected.append(record)
    selected_ids.add(id(record))
    return True


def _select_diverse_record(candidates, selected):
    """Return the candidate farthest from the current selected set."""
    best_candidate = None
    best_key = None
    for candidate in candidates:
        if not selected:
            nearest_selected = 0.0
        else:
            nearest_selected = min(
                _conformer_rmsd(candidate, chosen)
                for chosen in selected
            )
        key = (nearest_selected, -candidate.rdkit_energy, -candidate.seed, -candidate.source_conf_id)
        if best_key is None or key > best_key:
            best_key = key
            best_candidate = candidate
    return best_candidate, None if best_key is None else float(best_key[0])


def _selection_quota(limit, fraction):
    """Return a small protected quota for a selection fraction."""
    if limit <= 0 or float(fraction) <= 0.0:
        return 0
    return max(1, int(math.ceil(limit * float(fraction))))


def _select_refinement_records(records, limit, top_n, diversity_fraction, compactness_fraction):
    """Select backend candidates from energy, diversity, compact, and open basins."""
    limit = min(int(limit), len(records))
    if limit <= 0:
        return []
    diversity_fraction = float(diversity_fraction)
    compactness_fraction = float(compactness_fraction)

    # compactness_fraction is a protected folded/contact-rich quota, not a
    # global preference. Mirror it with an open-structure quota so compact
    # structures cannot crowd out extended conformer families.
    compact_count = min(limit, _selection_quota(limit, compactness_fraction))
    open_count = min(limit, _selection_quota(limit, compactness_fraction))
    diversity_count = min(limit, _selection_quota(limit, diversity_fraction))
    outlier_count = 1 if limit >= 10 else 0
    protected_count = compact_count + open_count + diversity_count + outlier_count
    energy_count = max(1, min(int(top_n), max(1, limit - protected_count)))
    energy_count = min(limit, energy_count)

    ordered_records = sorted(records, key=_record_sort_key)
    selected = []
    selected_ids = set()
    for record in ordered_records[:energy_count]:
        _append_selected_record(selected, selected_ids, record, "energy", 0.0)

    compact_candidates = sorted(
        (candidate for candidate in ordered_records if id(candidate) not in selected_ids),
        key=lambda record: (
            -_fold_score(record)[0],
            -_fold_score(record)[1],
            record.rdkit_energy,
            record.seed,
            record.source_conf_id,
        ),
    )
    for record in compact_candidates[: max(0, min(compact_count, limit - len(selected)))]:
        _append_selected_record(selected, selected_ids, record, "contact_compact", None)

    open_candidates = sorted(
        (candidate for candidate in ordered_records if id(candidate) not in selected_ids),
        key=lambda record: (
            -_open_score(record)[0],
            -_open_score(record)[1],
            record.rdkit_energy,
            record.seed,
            record.source_conf_id,
        ),
    )
    for record in open_candidates[: max(0, min(open_count, limit - len(selected)))]:
        _append_selected_record(selected, selected_ids, record, "open", None)

    for _ in range(max(0, min(diversity_count, limit - len(selected)))):
        candidates = [candidate for candidate in ordered_records if id(candidate) not in selected_ids]
        candidate, diversity = _select_diverse_record(candidates, selected)
        if candidate is None:
            break
        _append_selected_record(selected, selected_ids, candidate, "diversity", diversity)

    outlier_candidates = sorted(
        (candidate for candidate in ordered_records if id(candidate) not in selected_ids),
        key=lambda record: (
            -record.rdkit_energy,
            -_open_score(record)[0],
            -record.seed,
            -record.source_conf_id,
        ),
    )
    for record in outlier_candidates[: max(0, min(outlier_count, limit - len(selected)))]:
        _append_selected_record(selected, selected_ids, record, "outlier", None)

    while len(selected) < limit:
        candidates = [candidate for candidate in ordered_records if id(candidate) not in selected_ids]
        candidate, diversity = _select_diverse_record(candidates, selected)
        if candidate is None:
            break
        _append_selected_record(selected, selected_ids, candidate, "diversity", diversity)

    return sorted(selected, key=_record_sort_key)


def _refine_with_backend(records, run_directory, qc_params):
    """Optimize selected RDKit conformers with an existing PyAR backend."""
    backend_directory = run_directory / "backend"
    backend_directory.mkdir(parents=True, exist_ok=True)
    refined = []
    with _working_directory(backend_directory):
        for record in records:
            molecule = record.molecule.copy()
            status = optimise(molecule, qc_params)
            record.backend_status = _status_label(status)
            backend_result = backend_directory / f"job_{molecule.name}" / f"result_{molecule.name}.xyz"
            if backend_result.is_file():
                record.backend_path = str(backend_result)
            if is_usable(status) and molecule.coordinates is not None:
                record.molecule = molecule
                if molecule.energy is not None:
                    record.backend_energy = float(molecule.energy)
                refined.append(record)
    return sorted(
        refined,
        key=lambda record: (
            float("inf") if record.backend_energy is None else record.backend_energy,
            record.rdkit_energy,
            record.seed,
            record.source_conf_id,
        ),
    )


def _write_selected_conformers(records, selected_directory):
    """Write the final selected conformers in ranking order."""
    selected_paths = []
    for selected_rank, record in enumerate(records):
        path = selected_directory / f"conf_{selected_rank:04d}.xyz"
        record.molecule.mol_to_xyz(str(path))
        record.selected_path = str(path)
        selected_paths.append(str(path))
    return selected_paths


def _write_summary(records, summary_path):
    """Write a CSV summary for all accepted conformers."""
    columns = [
        "rank",
        "seed",
        "name",
        "source_conf_id",
        "force_field",
        "rdkit_status",
        "rdkit_energy",
        "rdkit_path",
        "generation_stage",
        "generation_round",
        "parent_name",
        "torsion_moves",
        "refinement_reason",
        "refinement_diversity",
        "contact_count",
        "radius_gyration",
        "backend_status",
        "backend_energy",
        "backend_path",
        "selected_path",
    ]
    with summary_path.open("w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(fp, fieldnames=columns)
        writer.writeheader()
        for record in records:
            writer.writerow(
                {
                    "rank": record.rank,
                    "seed": record.seed,
                    "name": record.name,
                    "source_conf_id": record.source_conf_id,
                    "force_field": record.force_field,
                    "rdkit_status": record.rdkit_status,
                    "rdkit_energy": record.rdkit_energy,
                    "rdkit_path": record.rdkit_path,
                    "generation_stage": record.generation_stage,
                    "generation_round": record.generation_round,
                    "parent_name": record.parent_name,
                    "torsion_moves": record.torsion_moves,
                    "refinement_reason": record.refinement_reason,
                    "refinement_diversity": record.refinement_diversity,
                    "contact_count": record.contact_count,
                    "radius_gyration": record.radius_gyration,
                    "backend_status": record.backend_status,
                    "backend_energy": record.backend_energy,
                    "backend_path": record.backend_path,
                    "selected_path": record.selected_path,
                }
            )


def _build_state(
    *,
    status,
    request,
    selected_paths,
    records,
    source_format,
    backend_requested,
):
    """Return the persisted conformer state payload."""
    return {
        "version": STATE_VERSION,
        "workflow": "conformer",
        "package_version": __version__,
        "status": status,
        "request": request,
        "source_format": source_format,
        "backend_requested": bool(backend_requested),
        "generated_conformers": len(records),
        "selected_results": selected_paths,
        "records": [
            {
                "rank": record.rank,
                "seed": record.seed,
                "name": record.name,
                "source_conf_id": record.source_conf_id,
                "force_field": record.force_field,
                "rdkit_status": record.rdkit_status,
                "rdkit_energy": record.rdkit_energy,
                "rdkit_path": record.rdkit_path,
                "generation_stage": record.generation_stage,
                "generation_round": record.generation_round,
                "parent_name": record.parent_name,
                "torsion_moves": record.torsion_moves,
                "refinement_reason": record.refinement_reason,
                "refinement_diversity": record.refinement_diversity,
                "contact_count": record.contact_count,
                "radius_gyration": record.radius_gyration,
                "backend_status": record.backend_status,
                "backend_energy": record.backend_energy,
                "backend_path": record.backend_path,
                "selected_path": record.selected_path,
            }
            for record in records
        ],
    }


def conformer_search(
    input_spec,
    *,
    input_format="auto",
    num_conformers=100,
    top_n=10,
    backend_top_n=None,
    num_seeds=3,
    diversity_fraction=0.2,
    compactness_fraction=0.2,
    rms_threshold=0.25,
    use_random_coords=True,
    torsion_kicks=True,
    torsion_rounds=2,
    torsion_mode="evolve",
    torsion_kicks_per_conformer=8,
    torsion_max_bonds=3,
    torsion_dedup_rms=0.5,
    force_field="auto",
    seed=1,
    num_threads=0,
    max_iterations=200,
    charge=None,
    multiplicity=1,
    scftype="rhf",
    qc_params=None,
    root_directory=None,
):
    """Generate and optionally refine conformers for one molecule."""
    if int(num_conformers) < 1:
        raise ConformerWorkflowError("--num-conformers must be at least 1")
    if int(top_n) < 1:
        raise ConformerWorkflowError("--top-n must be at least 1")
    if int(num_seeds) < 1:
        raise ConformerWorkflowError("--num-seeds must be at least 1")
    if float(rms_threshold) < 0.0:
        raise ConformerWorkflowError("--rms-threshold must be non-negative")
    if not isinstance(use_random_coords, bool):
        use_random_coords = bool(use_random_coords)
    if not isinstance(torsion_kicks, bool):
        torsion_kicks = bool(torsion_kicks)
    torsion_mode = str(torsion_mode).lower()
    if torsion_mode not in {"evolve", "mc", "grid", "random"}:
        raise ConformerWorkflowError("--torsion-mode must be 'evolve', 'mc', 'grid', or 'random'")
    if int(torsion_rounds) < 0:
        raise ConformerWorkflowError("--torsion-rounds must be non-negative")
    if int(torsion_kicks_per_conformer) < 0:
        raise ConformerWorkflowError("--torsion-kicks-per-conformer must be non-negative")
    if int(torsion_max_bonds) < 1:
        raise ConformerWorkflowError("--torsion-max-bonds must be at least 1")
    if float(torsion_dedup_rms) < 0.0:
        raise ConformerWorkflowError("--torsion-dedup-rms must be non-negative")
    if backend_top_n is not None and int(backend_top_n) < 1:
        raise ConformerWorkflowError("--backend-top-n must be at least 1")
    if not 0.0 <= float(diversity_fraction) <= 1.0:
        raise ConformerWorkflowError("--diversity-fraction must be between 0 and 1")
    if not 0.0 <= float(compactness_fraction) <= 1.0:
        raise ConformerWorkflowError("--compactness-fraction must be between 0 and 1")

    root_directory = os.getcwd() if root_directory is None else root_directory
    Chem, AllChem = _rdkit_modules()
    run_directory = _prepare_run_directory(root_directory)
    state_path = Path(workflow_state_path(str(run_directory)))
    seed_values = [int(seed) + index for index in range(int(num_seeds))]
    request = {
        "input": str(input_spec),
        "input_format": input_format,
        "num_conformers": int(num_conformers),
        "top_n": int(top_n),
        "backend_top_n": None if backend_top_n is None else int(backend_top_n),
        "num_seeds": int(num_seeds),
        "diversity_fraction": float(diversity_fraction),
        "compactness_fraction": float(compactness_fraction),
        "rms_threshold": float(rms_threshold),
        "use_random_coords": bool(use_random_coords),
        "torsion_kicks": bool(torsion_kicks),
        "torsion_rounds": int(torsion_rounds),
        "torsion_mode": torsion_mode,
        "torsion_kicks_per_conformer": int(torsion_kicks_per_conformer),
        "torsion_max_bonds": int(torsion_max_bonds),
        "torsion_dedup_rms": float(torsion_dedup_rms),
        "force_field": force_field,
        "seed": int(seed),
        "seed_values": seed_values,
        "num_threads": int(num_threads),
        "max_iterations": int(max_iterations),
        "charge": charge,
        "multiplicity": int(multiplicity),
        "scftype": scftype,
        "backend_parameters": dict(qc_params or {}),
    }

    rdkit_molecule, source_format = _load_rdkit_molecule(
        input_spec,
        input_format,
        Chem,
        run_directory,
    )
    records = []
    for seed_value in seed_values:
        try:
            seed_molecule = Chem.Mol(rdkit_molecule)
        except (AttributeError, TypeError):
            seed_molecule = rdkit_molecule
        seed_records = _embed_and_rank_conformers(
            seed_molecule,
            AllChem,
            seed=seed_value,
            num_conformers=num_conformers,
            prune_rms_threshold=rms_threshold,
            force_field=force_field,
            num_threads=num_threads,
            use_random_coords=use_random_coords,
            max_iterations=max_iterations,
        )
        records.extend(seed_records)
    records = sorted(records, key=lambda record: (record.rdkit_energy, record.seed, record.source_conf_id))
    if records and torsion_kicks and int(torsion_rounds) > 0:
        torsion_parent_limit = max(int(top_n) * 10, 50)
        if torsion_mode == "evolve":
            kicked_records = _generate_torsion_evolution_records(
                records,
                Chem,
                AllChem,
                enabled=torsion_kicks,
                generations=torsion_rounds,
                kicks_per_conformer=torsion_kicks_per_conformer,
                max_bonds=torsion_max_bonds,
                max_iterations=max_iterations,
                diversity_fraction=diversity_fraction,
                compactness_fraction=compactness_fraction,
                parent_limit=torsion_parent_limit,
                tabu_rms=torsion_dedup_rms,
            )
            if kicked_records:
                records = _deduplicate_records(records + kicked_records, torsion_dedup_rms)
        else:
            for round_index in range(1, int(torsion_rounds) + 1):
                parent_records = _select_torsion_parent_records(
                    records,
                    torsion_parent_limit,
                    diversity_fraction,
                    compactness_fraction,
                )
                kicked_records = (
                    _generate_torsion_mc_records(
                        parent_records,
                        Chem,
                        AllChem,
                        enabled=torsion_kicks,
                        kicks_per_conformer=torsion_kicks_per_conformer,
                        max_bonds=torsion_max_bonds,
                        seed=seed,
                        round_index=round_index,
                        max_iterations=max_iterations,
                        tabu_rms=torsion_dedup_rms,
                    )
                    if torsion_mode == "mc"
                    else _generate_torsion_grid_records(
                        parent_records,
                        Chem,
                        AllChem,
                        enabled=torsion_kicks,
                        kicks_per_conformer=torsion_kicks_per_conformer,
                        max_bonds=torsion_max_bonds,
                        round_index=round_index,
                        max_iterations=max_iterations,
                    )
                    if torsion_mode == "grid"
                    else _generate_torsion_random_records(
                        parent_records,
                        Chem,
                        AllChem,
                        enabled=torsion_kicks,
                        kicks_per_conformer=torsion_kicks_per_conformer,
                        max_bonds=torsion_max_bonds,
                        seed=seed,
                        round_index=round_index,
                        max_iterations=max_iterations,
                    )
                )
                if not kicked_records:
                    break
                records = _deduplicate_records(records + kicked_records, torsion_dedup_rms)
    records = sorted(records, key=lambda record: (record.rdkit_energy, record.seed, record.source_conf_id))
    if not records:
        state = _build_state(
            status="failed_no_conformers",
            request=request,
            selected_paths=[],
            records=[],
            source_format=source_format,
            backend_requested=bool(qc_params and qc_params.get("software")),
        )
        _write_json_atomic(state_path, state)
        raise ConformerWorkflowError("RDKit did not generate any conformers")

    _write_rdkit_conformers(
        records,
        rdkit_molecule,
        run_directory / "rdkit",
        charge=charge,
        multiplicity=multiplicity,
        scftype=scftype,
    )

    if backend_top_n is None:
        refinement_limit = max(int(top_n) * 5, int(top_n) + 20)
    else:
        refinement_limit = int(backend_top_n)
    refinement_limit = min(refinement_limit, len(records))
    backend_requested = bool(qc_params and qc_params.get("software"))
    if backend_requested:
        retained_records = _select_refinement_records(
            records,
            refinement_limit,
            top_n,
            diversity_fraction,
            compactness_fraction,
        )
        refined_records = _refine_with_backend(retained_records, run_directory, qc_params)
        selected_records = _select_unique_records(
            refined_records,
            min(int(top_n), len(refined_records)),
            max(float(torsion_dedup_rms), 0.75),
        )
        status = "completed" if selected_records else "completed_no_backend_success"
    else:
        retained_records = _select_refinement_records(
            records,
            min(int(top_n), len(records)),
            top_n,
            diversity_fraction,
            compactness_fraction,
        )
        selected_records = _select_unique_records(
            retained_records,
            min(int(top_n), len(retained_records)),
            float(torsion_dedup_rms),
        )
        status = "completed"

    selected_paths = _write_selected_conformers(selected_records, run_directory / "selected")
    _write_summary(records, run_directory / "summary.csv")
    state = _build_state(
        status=status,
        request=request,
        selected_paths=selected_paths,
        records=records,
        source_format=source_format,
        backend_requested=backend_requested,
    )
    _write_json_atomic(state_path, state)

    conformer_logger.info(
        "Conformer workflow completed: generated=%d selected=%d backend=%s",
        len(records),
        len(selected_records),
        "yes" if backend_requested else "no",
    )
    return ConformerResult(
        workflow="conformer",
        status=status,
        run_directory=str(run_directory),
        state_path=str(state_path),
        selected_paths=tuple(selected_paths),
        metadata={
            "source_format": source_format,
            "generated_conformers": len(records),
            "selected_conformers": len(selected_records),
            "summary_path": str(run_directory / "summary.csv"),
            "backend_requested": backend_requested,
            "seed_values": seed_values,
            "backend_top_n": refinement_limit,
            "diversity_fraction": float(diversity_fraction),
            "compactness_fraction": float(compactness_fraction),
            "use_random_coords": bool(use_random_coords),
            "torsion_kicks": bool(torsion_kicks),
            "torsion_rounds": int(torsion_rounds),
            "torsion_mode": torsion_mode,
            "torsion_kicks_per_conformer": int(torsion_kicks_per_conformer),
            "torsion_max_bonds": int(torsion_max_bonds),
            "torsion_dedup_rms": float(torsion_dedup_rms),
        },
    )
