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
    rank: int | None = None
    name: str | None = None
    molecule: Molecule | None = None
    rdkit_path: str | None = None
    refinement_reason: str | None = None
    refinement_diversity: float | None = None
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
            )
        )
    return sorted(records, key=lambda record: (record.rdkit_energy, record.source_conf_id))


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
        name = f"seed_{record.seed:06d}_conf_{record.source_conf_id:04d}"
        source_molecule = record.rdkit_molecule if record.rdkit_molecule is not None else rdkit_molecule
        molecule = _molecule_from_rdkit_conformer(
            source_molecule,
            record.source_conf_id,
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
    if molecule is None or molecule.coordinates is None:
        return None
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

    return contact_count, rgyr


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


def _select_refinement_records(records, limit, top_n, diversity_fraction, compactness_fraction):
    """Select backend candidates by RDKit energy plus geometric diversity."""
    limit = min(int(limit), len(records))
    if limit <= 0:
        return []
    diversity_fraction = float(diversity_fraction)
    compactness_fraction = float(compactness_fraction)
    energy_count = math.ceil(limit * (1.0 - diversity_fraction))
    energy_count = max(1, int(top_n), min(limit, energy_count))
    energy_count = min(limit, energy_count)

    ordered_records = sorted(records, key=_record_sort_key)
    selected = list(ordered_records[:energy_count])
    selected_ids = {id(record) for record in selected}
    for record in selected:
        record.refinement_reason = "energy"
        record.refinement_diversity = 0.0

    compact_count = max(0, min(limit - len(selected), math.ceil((limit - len(selected)) * compactness_fraction)))
    compact_candidates = []
    if compact_count > 0:
        compact_candidates = sorted(
            (candidate for candidate in ordered_records if id(candidate) not in selected_ids),
            key=lambda record: (
                -_compactness_metrics(record)[0],
                _compactness_metrics(record)[1],
                record.rdkit_energy,
                record.seed,
                record.source_conf_id,
            ),
        )[:compact_count]
        for record in compact_candidates:
            record.refinement_reason = "compact"
            record.refinement_diversity = None
            selected.append(record)
            selected_ids.add(id(record))

    while len(selected) < limit:
        best_candidate = None
        best_key = None
        for candidate in ordered_records:
            if id(candidate) in selected_ids:
                continue
            nearest_selected = min(
                _conformer_rmsd(candidate, chosen)
                for chosen in selected
            )
            key = (nearest_selected, -candidate.rdkit_energy, -candidate.seed, -candidate.source_conf_id)
            if best_key is None or key > best_key:
                best_key = key
                best_candidate = candidate
        if best_candidate is None:
            break
        best_candidate.refinement_reason = "diversity"
        best_candidate.refinement_diversity = float(best_key[0])
        selected.append(best_candidate)
        selected_ids.add(id(best_candidate))

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
        "refinement_reason",
        "refinement_diversity",
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
                    "refinement_reason": record.refinement_reason,
                    "refinement_diversity": record.refinement_diversity,
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
                "refinement_reason": record.refinement_reason,
                "refinement_diversity": record.refinement_diversity,
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
    num_seeds=1,
    diversity_fraction=0.5,
    compactness_fraction=0.5,
    rms_threshold=0.5,
    use_random_coords=False,
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

    refinement_limit = int(top_n) if backend_top_n is None else max(int(top_n), int(backend_top_n))
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
        selected_records = refined_records[: min(int(top_n), len(refined_records))]
        status = "completed" if selected_records else "completed_no_backend_success"
    else:
        retained_records = records[: min(refinement_limit, len(records))]
        selected_records = retained_records[: min(int(top_n), len(retained_records))]
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
        },
    )
