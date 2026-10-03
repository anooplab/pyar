"""Constrained relaxed scans using registered energy/gradient backends."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import tempfile

import numpy as np

from pyar.backends.orca_scan import OrcaBondScanRequest as BondScanRequest
from pyar.backends.orca_scan import OrcaBondScanResult as BondScanResult
from pyar.backends.orca_scan import parse_xyz_trajectory


def _hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_bond_scan_result(molecule, request, directory):
    """Recover only complete, hash-verified generic scan artifacts."""
    directory = Path(directory)
    summary = directory / "scan_summary.json"
    trajectory = directory / "scan_trajectory.xyz"
    final = directory / "final_scan.xyz"
    failure = BondScanResult(False, None, None, None, directory / "scan_request.json",
                             summary, "scan_incomplete")
    try:
        scan_request = json.loads((directory / "scan_request.json").read_text())
        if scan_request["scan"] != vars(request):
            return failure
        if scan_request["symbols"] != list(molecule.atoms_list):
            return failure
        if scan_request["coordinates"] != np.asarray(molecule.coordinates).tolist():
            return failure
        if any(scan_request["qc_params"].get(key) != getattr(molecule, key)
               for key in ("charge", "multiplicity", "scftype")):
            return failure
        state = json.loads(summary.read_text())
        if state["status"] != "success" or state["qc_params"] != scan_request["qc_params"]:
            return failure
        if set(state["artifacts"]) != {"scan_trajectory.xyz", "final_scan.xyz", "scan_request.json"}:
            return failure
        if any(_hash(directory / name) != digest for name, digest in state["artifacts"].items()):
            return failure
        frames = parse_xyz_trajectory(trajectory, molecule.number_of_atoms, molecule.atoms_list)
        profile = state["profile"]
        grid = np.linspace(request.start_distance_angstrom, request.end_distance_angstrom, request.n_points)
        if len(frames) != request.n_points or len(profile) != request.n_points:
            return failure
        for (_, xyz, comment), point, target in zip(frames, profile, grid):
            stored_energy = float(comment.split("energy_hartree=", 1)[1].split()[0])
            if not np.isclose(stored_energy, point["energy_hartree"], atol=1e-12, rtol=0):
                return failure
            distance = np.linalg.norm(xyz[request.atom_i] - xyz[request.atom_j])
            if (not np.isfinite(point["energy_hartree"])
                    or not np.isclose(point["target_distance_angstrom"], target, atol=1e-8, rtol=0)
                    or not np.isclose(distance, target, atol=1e-3, rtol=0)):
                return failure
        stored = parse_xyz_trajectory(final, molecule.number_of_atoms, molecule.atoms_list)
        if len(stored) != 1 or not np.array_equal(stored[0][1], frames[-1][1]):
            return failure
        return BondScanResult(True, trajectory, final, frames[-1][1],
                              directory / "scan_request.json", summary, "success", frames, profile)
    except (OSError, ValueError, KeyError, TypeError, IndexError):
        return failure


def run_bond_scan(molecule, request, directory, qc_params):
    """Run with calculator state and optimizer side files isolated from the caller."""
    directory = Path(directory).resolve()
    directory.mkdir(parents=True, exist_ok=True)
    previous = Path.cwd()
    try:
        os.chdir(directory)
        return _run_bond_scan(molecule, request, directory, qc_params)
    finally:
        os.chdir(previous)


def _run_bond_scan(molecule, request, directory, qc_params):
    """Optimize each distance with geomeTRIC, seeding it from the previous point.

    Constrained-coordinate gradients are handled by geomeTRIC; Cartesian forces
    along the fixed bond are not incorrectly required to vanish.
    """
    from ase.units import Bohr
    from geometric.errors import GeomOptNotConvergedError
    from geometric.internal import DelocalizedInternalCoordinates, Distance
    from geometric.optimize import OPT_STATE, Optimizer
    from geometric.params import OptParams

    from pyar.backends.geometric import PyarGeometricCalculator
    from pyar.neb import _engine, _write_xyz_trajectory

    if (isinstance(request.n_points, bool) or not isinstance(request.n_points, int)
            or request.n_points < 2):
        raise ValueError("scan requires at least two integer points")
    if (request.atom_i == request.atom_j
            or any(isinstance(index, bool) or not isinstance(index, int)
                   or index < 0 or index >= molecule.number_of_atoms
                   for index in (request.atom_i, request.atom_j))):
        raise ValueError("scan requires two distinct valid atom indices")
    if (not np.isfinite([request.start_distance_angstrom, request.end_distance_angstrom]).all()
            or min(request.start_distance_angstrom, request.end_distance_angstrom) <= 0
            or request.start_distance_angstrom == request.end_distance_angstrom):
        raise ValueError("scan endpoints must be distinct finite positive distances")
    directory = Path(directory).resolve()
    directory.mkdir(parents=True, exist_ok=True)
    params = dict(qc_params, charge=molecule.charge, multiplicity=molecule.multiplicity,
                  scftype=molecule.scftype, gamma=0.0)
    calculator = PyarGeometricCalculator(params)
    request_path = directory / "scan_request.json"
    request_path.write_text(json.dumps({"scan": vars(request), "qc_params": params,
                                        "symbols": list(molecule.atoms_list),
                                        "coordinates": np.asarray(molecule.coordinates).tolist()}, indent=2))
    summary_path = directory / "scan_summary.json"
    summary_path.write_text(json.dumps({"status": "running"}))
    frames, energies, profile = [], [], []
    xyz = np.asarray(molecule.coordinates, dtype=float).copy()
    threshold = {"loose": "GAU_LOOSE", "normal": "GAU", "tight": "GAU_TIGHT"}[
        params.get("opt_threshold", "normal")]
    for index, target in enumerate(np.linspace(request.start_distance_angstrom,
                                               request.end_distance_angstrom, request.n_points)):
        # Shift the complete second fragment before optimization, preserving its
        # internal geometry and avoiding a large initial constraint violation.
        vector = xyz[request.atom_j] - xyz[request.atom_i]
        distance = np.linalg.norm(vector)
        if not np.isfinite(distance) or distance <= 0:
            raise ValueError("scan target atoms must have a finite nonzero separation")
        group = next((fragment for fragment in (molecule.fragments or [])
                      if request.atom_j in fragment and request.atom_i not in fragment), [request.atom_j])
        xyz[group] += (float(target) - distance) * vector / distance
        gmolecule, engine = _engine(molecule.atoms_list, xyz, calculator)
        internal = DelocalizedInternalCoordinates(
            gmolecule, build=True, connect=False, addcart=False,
            # Orthogonalize the constrained subspace so unconstrained bonds
            # retain their correct minimization directions.
            constraints=[Distance(request.atom_i, request.atom_j)],
            cvals=[float(target) / Bohr], conmethod=1,
        )
        scratch = tempfile.mkdtemp(prefix=f"point_{index:03d}_", dir=directory)
        optimizer = Optimizer(xyz.flatten() / Bohr, gmolecule, internal, engine, scratch,
                              OptParams(maxiter=params.get("opt_cycles", 100),
                                        convergence_set=threshold, enforce=0.1, frequency=False),
                              print_info=False)
        try:
            progress = optimizer.optimizeGeometry()
        except GeomOptNotConvergedError:
            progress = optimizer.progress
        xyz = np.asarray(progress.xyzs[-1]).copy()
        energy = float(progress.qm_energies[-1])
        valid = (optimizer.state == OPT_STATE.CONVERGED and np.isfinite(xyz).all()
                 and np.isfinite(energy)
                 and np.isclose(np.linalg.norm(xyz[request.atom_i] - xyz[request.atom_j]),
                                target, atol=1e-3, rtol=0))
        if not valid:
            summary_path.write_text(json.dumps({"status": "point_not_converged", "point": index}))
            return BondScanResult(False, None, None, None, request_path, summary_path, "point_not_converged")
        frames.append((list(molecule.atoms_list), xyz.copy(), ""))
        energies.append(energy)
        profile.append({"target_distance_angstrom": float(target), "energy_hartree": energy})
    trajectory = directory / "scan_trajectory.xyz"
    _write_xyz_trajectory(trajectory, molecule.atoms_list, [frame[1] for frame in frames], energies)
    final = directory / "final_scan.xyz"
    _write_xyz_trajectory(final, molecule.atoms_list, [xyz], [energies[-1]])
    summary_path.write_text(json.dumps({"status": "success", "profile": profile,
                                       "optimizer": "geometric", "qc_params": params,
                                       "artifacts": {path.name: _hash(path) for path in (trajectory, final, request_path)}}, indent=2))
    return load_bond_scan_result(molecule, request, directory)
