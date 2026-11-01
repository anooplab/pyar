"""Controlled backend fixtures and provenance-preserving golden normalization."""

from collections import Counter
import hashlib
import importlib
import json
from pathlib import Path
import random
import re
import sys
import uuid

import numpy as np

BASELINE = "9f1e567ce22b5d382d1facc9ff1da624de5d55a9"
VOLATILE_KEYS = {
    "created_at",
    "updated_at",
    "timestamp",
    "created_epoch",
    "deadline_epoch",
    "elapsed_seconds",
    "started_at",
    "finished_at",
    "acquired_utc",
    "released_utc",
    "recovered_utc",
    "started_at_epoch",
}


def digest(value):
    return hashlib.sha256(value).hexdigest()


class DeterministicBackend:
    """Records explicit chemistry inputs; does not establish scientific truth."""

    def __init__(self, monkeypatch):
        self.calls = []
        self.scans = {}
        for module in (
            "pyar.optimiser",
            "pyar.growth.service",
            "pyar.workflows.conformer",
            "pyar.workflows.reaction",
            "pyar.workflows.scan_bond",
        ):
            monkeypatch.setattr(
                importlib.import_module(module), "optimise", self.optimize
            )
        monkeypatch.setattr(
            "pyar.energy_gradient_providers.get_energy_gradient_provider",
            lambda *args, **kwargs: self,
        )
        monkeypatch.setattr("pyar.backends.bond_scan.run_bond_scan", self.scan)
        monkeypatch.setattr(
            "pyar.backends.bond_scan.load_bond_scan_result",
            lambda molecule, request, directory: self.scans.get(str(directory)),
        )

    def record(self, kind, molecule, settings):
        self.calls.append(
            dict(
                operation=kind,
                name=molecule.name,
                symbols=list(molecule.atoms_list),
                charge=molecule.charge,
                multiplicity=molecule.multiplicity,
                settings=dict(settings),
                coordinates_sha256=digest(
                    np.asarray(molecule.coordinates, dtype="<f8").tobytes()
                ),
            )
        )

    def optimize(self, molecule, settings, *, directory=None):
        self.record("optimize", molecule, settings)
        if settings.get("microsolvation_evaluation_budget_path"):
            from pyar.evaluation_budget import reserve_backend_evaluation

            reserve_backend_evaluation(
                settings["microsolvation_evaluation_budget_path"]
            )
        directory = Path(directory or Path.cwd()) / f"job_{molecule.name}"
        directory.mkdir(parents=True, exist_ok=True)
        molecule.mol_to_xyz(directory / f"trial_{molecule.name}.xyz")
        molecule.energy = -float(len(molecule))
        molecule.optimized_coordinates = molecule.coordinates.copy()
        molecule.mol_to_xyz(directory / f"result_{molecule.name}.xyz")
        return True

    def evaluate(self, molecule, coordinates_bohr):
        from pyar.energy_gradient_providers import EnergyGradientResult

        self.record("evaluate", molecule, {})
        return EnergyGradientResult(
            -float(len(molecule)), np.zeros_like(coordinates_bohr)
        )

    def scan(self, molecule, request, directory, settings):
        from pyar.backends.orca_scan import OrcaBondScanResult
        from pyar.neb import _write_xyz_trajectory

        self.record("scan", molecule, settings)
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)
        grid = np.linspace(
            request.start_distance_angstrom,
            request.end_distance_angstrom,
            request.n_points,
        )
        delta = (
            molecule.coordinates[request.atom_j] - molecule.coordinates[request.atom_i]
        )
        direction = delta / np.linalg.norm(delta)
        frames, profile = [], []
        for distance in grid:
            xyz = molecule.coordinates.copy()
            xyz[request.atom_j] = xyz[request.atom_i] + distance * direction
            frames.append((list(molecule.atoms_list), xyz, ""))
            profile.append(
                dict(target_distance_angstrom=float(distance), energy_hartree=-1.0)
            )
        trajectory, final = (
            directory / "scan_trajectory.xyz",
            directory / "final_scan.xyz",
        )
        _write_xyz_trajectory(
            trajectory,
            molecule.atoms_list,
            [f[1] for f in frames],
            [-1.0] * len(frames),
        )
        _write_xyz_trajectory(final, molecule.atoms_list, [frames[-1][1]], [-1.0])
        request_path, output_path = (
            directory / "scan_request.json",
            directory / "scan_summary.json",
        )
        request_path.write_text(json.dumps(vars(request), sort_keys=True))
        output_path.write_text('{"status": "success"}')
        result = OrcaBondScanResult(
            True,
            trajectory,
            final,
            frames[-1][1],
            request_path,
            output_path,
            "success",
            frames,
            profile,
        )
        self.scans[str(directory)] = result
        return result


def controlled_conformers(monkeypatch):
    """Fixed embedding evidence keeps the golden independent of RDKit minimizers."""
    from pyar.workflows import conformer

    def embed(molecule, AllChem, *, seed, **options):
        from rdkit import Chem

        geometry = Chem.Conformer(molecule.GetNumAtoms())
        xyz = [
            (0, 0, 0),
            (0.629, 0.629, 0.629),
            (-0.629, -0.629, 0.629),
            (-0.629, 0.629, -0.629),
            (0.629, -0.629, -0.629),
        ]
        assert molecule.GetNumAtoms() == len(xyz)
        for index, position in enumerate(xyz):
            geometry.SetAtomPosition(index, position)
        index = molecule.AddConformer(geometry, assignId=True)
        return [
            conformer.ConformerRecord(
                seed=seed,
                source_conf_id=index,
                rdkit_energy=1.0,
                rdkit_status="converged",
                force_field="uff",
                rdkit_molecule=molecule,
                rdkit_conf_id=index,
            )
        ]

    monkeypatch.setattr(conformer, "_embed_and_rank_conformers", embed)


def deterministic_metadata(monkeypatch):
    """Control UUID generation without weakening reference-identity comparisons."""
    original = uuid.uuid4
    aliases = [
        (module, name)
        for module_name, module in tuple(sys.modules.items())
        if module_name.startswith("pyar") and module is not None
        for name, value in tuple(vars(module).items())
        if value is original or getattr(value, "_pyar_golden_uuid", False) is True
    ]
    sequence = iter(range(1, 100000))

    def next_uuid():
        return uuid.UUID(int=next(sequence), version=4)

    next_uuid._pyar_golden_uuid = True
    monkeypatch.setattr(uuid, "uuid4", next_uuid)
    for module, name in aliases:
        monkeypatch.setattr(module, name, next_uuid)


def run_scenario(kind, root, monkeypatch, backend=None):
    original_uuid = uuid.uuid4
    try:
        return _run_scenario(kind, root, monkeypatch, backend)
    finally:
        # Modules first imported while metadata is controlled must not retain
        # the test generator after monkeypatch restores the uuid module.
        for name, module in tuple(sys.modules.items()):
            if name.startswith("pyar") and module is not None:
                for key, value in tuple(vars(module).items()):
                    if getattr(value, "_pyar_golden_uuid", False) is True:
                        setattr(module, key, original_uuid)


def _run_scenario(kind, root, monkeypatch, backend=None):
    from pyar.core.molecule import Molecule

    deterministic_metadata(monkeypatch)
    root.mkdir(parents=True)
    monkeypatch.chdir(root)
    random.seed(43)
    np.random.seed(43)
    backend = backend or DeterministicBackend(monkeypatch)
    qc = dict(software="xtb", opt_threshold="normal", nprocs=1)

    def atom(symbol, name):
        return Molecule([symbol], [[0, 0, 0]], name=name)

    if kind == "conformer":
        from pyar.workflows.conformer import conformer_search

        controlled_conformers(monkeypatch)
        result = conformer_search(
            "C",
            num_conformers=2,
            top_n=1,
            num_seeds=1,
            torsion_kicks=False,
            seed=43,
            num_threads=1,
            qc_params=qc,
            root_directory=root,
        )
    elif kind == "aggregate":
        from pyar.workflows.aggregate import aggregate

        result = aggregate(
            [atom("He", "A"), atom("Ne", "B")],
            [1, 1],
            2,
            qc,
            1,
            0,
            1,
            None,
            connectivity_policy="off",
            selection_feature="distance-histogram",
        )
    elif kind == "grow":
        from pyar.growth.request import GrowRequest
        from pyar.workflows.grow import grow

        result = grow(
            GrowRequest(
                atom("He", "seed"),
                atom("Ne", "monomer"),
                1,
                number_of_orientations=2,
                backend_parameters=qc,
                maximum_number_of_seeds=1,
                connectivity_policy="off",
                selection_feature="distance-histogram",
            ),
            output=root / "grow",
        )
    elif kind == "microsolvation":
        from pyar.microsolvation.request import MicrosolvationRequest
        from pyar.workflows.microsolvation import microsolvate

        solute = Molecule(
            ["C", "H", "H", "H", "H"],
            [
                [0, 0, 0],
                [0.629, 0.629, 0.629],
                [-0.629, -0.629, 0.629],
                [-0.629, 0.629, -0.629],
                [0.629, -0.629, -0.629],
            ],
            name="methane",
        )
        water = Molecule(
            ["O", "H", "H"],
            [[0, 0, 0], [0.9572, 0, 0], [-0.239, 0.927, 0]],
            name="water",
        )
        request = MicrosolvationRequest(
            solute, water, 1, 2, qc, 1, None, surface_points_per_atom=16
        )
        result = microsolvate(request, output=root / "microsolvation")
    elif kind == "microsolvation_resampling":
        from pyar.microsolvation.request import MicrosolvationRequest
        from pyar.workflows.microsolvation import microsolvate

        request = MicrosolvationRequest(
            atom("He", "solute"),
            atom("He", "solvent"),
            1,
            2,
            dict(qc, geometry_optimizer="geometric"),
            2,
            surface_points_per_atom=8,
            objective="cluster",
            placement_score="energy",
            resampling_moves=2,
            resampling_evaluation_budget=20,
            resampling_time_limit=30,
            resampling_seed=12,
        )
        result = microsolvate(request, output=root / "microsolvation")
    elif kind == "reaction":
        from pyar.workflows.reaction import react

        result = react(atom("He", "A"), atom("Ne", "B"), 100.0, 100.0, 2, qc, None, 2.3)
    else:
        from pyar.backends import write_xyz
        from pyar.workflows.scan_bond import run_scan_bond

        first, second = root / "first.xyz", root / "second.xyz"
        if kind == "scan_inter":
            write_xyz(["H"], [[0, 0, 0]], first)
            write_xyz(["H"], [[0, 0, 0]], second)
            result = run_scan_bond(
                first.name,
                second.name,
                (0, 0),
                1,
                qc,
                root / "scan",
                scan_end=1.0,
                scan_points=3,
                through="scan",
            )
        else:
            write_xyz(["H", "H"], [[0, 0, 0], [0.74, 0, 0]], first)
            result = run_scan_bond(
                first.name,
                None,
                (0, 1),
                1,
                qc,
                root / "scan",
                scan_end=3.0,
                scan_points=3,
                through="scan",
            )
    value = result.to_dict() if hasattr(result, "to_dict") else result
    return snapshot(root, value, backend.calls)


def snapshot(root, result, calls):
    """Full JSON/CSV and exact selected-XYZ hashes; normalize metadata only.

    JSON artifact-reference hashes are recomputed over the normalized referred
    document, recursively, so timestamps/absolute-root changes do not hide a
    changed scientific field. XYZ/binary hashes are never normalized.
    """
    files = sorted(path for path in root.rglob("*") if path.is_file())
    replacements = {str(root): "<ROOT>"}
    identities = Counter()
    for path in files:
        if re.fullmatch(r"[0-9a-f]{32}\.json", path.name):
            record = json.loads(path.read_text())
            identity = f"{record['workflow_role']}-{record['molecule_name']}"
            identities[identity] += 1
            replacements[path.stem] = f"{identity}-{identities[identity]}"
    raw = {path: path.read_bytes() for path in files}
    hashes = {}

    def normalize(value, key=None):
        if key in VOLATILE_KEYS:
            return "<TIME>"
        if isinstance(value, dict):
            if value.get("format") == "pyar-run-lock":
                value = dict(value, pid="<PID>", host="<HOST>")
            return {k: normalize(v, k) for k, v in value.items()}
        if isinstance(value, (list, tuple)):
            return [normalize(v) for v in value]
        if isinstance(value, str):
            value = re.sub(
                r"[0-9a-f]{64}", lambda match: hashes.get(match[0], match[0]), value
            )
            for before, after in replacements.items():
                value = value.replace(before, after)
        return value

    documents = {}
    for path, data in raw.items():
        if path.suffix in {".json", ".lock"}:
            documents[path] = json.loads(data)
    # Recompute path-bearing request fingerprints rather than ignoring them.
    from pyar.recovery.models import OptimizationRequest

    request_hashes = {}

    def request_fingerprints(value):
        if isinstance(value, dict):
            if {
                "request_fingerprint",
                "scientific_fingerprint",
                "settings",
                "symbols",
            } <= value.keys():

                def request(settings):
                    return OptimizationRequest(
                        value["name"],
                        tuple(value["symbols"]),
                        value["coordinates_angstrom"],
                        value["backend"],
                        value["optimizer"],
                        value["target"],
                        json.dumps(settings),
                        json.dumps(value["asset_hashes"]),
                        value["workflow_role"],
                    )

                original, normalized = (
                    request(value["settings"]),
                    request(normalize(value["settings"])),
                )
                for field in ("request_fingerprint", "scientific_fingerprint"):
                    assert value[field] == getattr(original, field), (
                        "Invalid original request identity"
                    )
                    request_hashes[value[field]] = getattr(normalized, field)
            for child in value.values():
                request_fingerprints(child)
        elif isinstance(value, list):
            for child in value:
                request_fingerprints(child)

    for value in documents.values():
        request_fingerprints(value)
    hashes.update(request_hashes)
    for _ in range(len(documents) + 1):
        updated = {
            digest(raw[path]): digest(
                json.dumps(
                    normalize(value), sort_keys=True, separators=(",", ":")
                ).encode()
            )
            for path, value in documents.items()
        }
        updated.update(request_hashes)
        if updated == hashes:
            break
        hashes = updated
    inventory = {}
    for path, data in raw.items():
        name = normalize(str(path.relative_to(root)))
        if path in documents:
            inventory[name] = normalize(documents[path])
        elif path.suffix == ".jsonl":
            inventory[name] = [
                normalize(json.loads(line)) for line in data.decode().splitlines()
            ]
        elif path.suffix == ".xyz":
            inventory[name] = dict(sha256=digest(data), text=data.decode())
        else:
            inventory[name] = normalize(data.decode())
    return dict(
        result=normalize(result), backend_calls=normalize(calls), files=inventory
    )
