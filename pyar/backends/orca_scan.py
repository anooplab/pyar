"""ORCA relaxed bond-scan input, execution, and trajectory recovery."""

from __future__ import annotations

from dataclasses import dataclass
import math
import os
from pathlib import Path

import numpy as np

from pyar.backends import require_executable, write_xyz
from pyar.backends.subprocess_utils import run_command


@dataclass(frozen=True)
class OrcaBondScanRequest:
    atom_i: int
    atom_j: int
    start_distance_angstrom: float
    end_distance_angstrom: float
    n_points: int


@dataclass
class OrcaBondScanResult:
    success: bool
    trajectory_path: Path | None
    final_geometry_path: Path | None
    final_coordinates: np.ndarray | None
    input_path: Path
    output_path: Path
    status: str
    frames: list | None = None
    profile: list | None = None


def parse_xyz_trajectory(path, expected_atoms, expected_symbols=None):
    """Parse a concatenated XYZ trajectory and return all frames."""
    lines = Path(path).read_text().splitlines()
    frames = []
    cursor = 0
    while cursor < len(lines):
        # ORCA's ``scan.allxyz`` separates frames with a standalone ``>``.
        while cursor < len(lines) and lines[cursor].strip() in {"", ">"}:
            cursor += 1
        if cursor >= len(lines):
            break
        try:
            count = int(lines[cursor].strip())
        except ValueError as exc:
            raise ValueError(f"Invalid XYZ atom count in {path} at line {cursor + 1}") from exc
        if count != expected_atoms:
            raise ValueError(f"XYZ frame has {count} atoms; expected {expected_atoms}")
        if cursor + count + 2 > len(lines):
            raise ValueError(f"Truncated XYZ frame in {path}")
        symbols = []
        coordinates = []
        for line in lines[cursor + 2:cursor + count + 2]:
            fields = line.split()
            if len(fields) < 4:
                raise ValueError(f"Malformed XYZ coordinate line in {path}")
            symbols.append(fields[0].capitalize())
            try:
                coordinate = [float(value) for value in fields[1:4]]
            except ValueError as exc:
                raise ValueError(f"Malformed XYZ coordinates in {path}") from exc
            if not all(math.isfinite(value) for value in coordinate):
                raise ValueError(f"Non-finite XYZ coordinates in {path}")
            coordinates.append(coordinate)
        if expected_symbols is not None and symbols != list(expected_symbols):
            raise ValueError("XYZ trajectory atom symbols do not match the scan input")
        frames.append((symbols, np.asarray(coordinates, dtype=float), lines[cursor + 1]))
        cursor += count + 2
    if not frames:
        raise ValueError(f"XYZ trajectory contains no frames: {path}")
    return frames


def parse_orca_scan_profile(path):
    """Parse ORCA's two-column relaxed-scan distance/actual-energy table."""
    profile = []
    for line_number, line in enumerate(Path(path).read_text().splitlines(), start=1):
        fields = line.split()
        if not fields or fields[0].startswith(("#", "!")):
            continue
        if len(fields) != 2:
            raise ValueError(f"Malformed ORCA scan profile row at line {line_number}")
        try:
            distance, energy = map(float, fields)
        except ValueError as exc:
            raise ValueError(f"Malformed ORCA scan profile row at line {line_number}") from exc
        if not math.isfinite(distance) or distance <= 0 or not math.isfinite(energy):
            raise ValueError(f"Non-finite or invalid ORCA scan profile value at line {line_number}")
        profile.append({"target_distance_angstrom": distance, "energy_hartree": energy})
    if not profile:
        raise ValueError(f"ORCA scan profile contains no points: {path}")
    return profile


def _orca_keyword(qc_params, scftype=None):
    threshold = {"loose": "LooseOpt", "tight": "TightOpt"}.get(
        qc_params.get("opt_threshold"), "Opt"
    )
    keyword = f"! {qc_params['method']} {qc_params['basis']} {threshold} RI def2/J D3BJ KDIIS"
    # ``merged_with`` represents unrestricted molecules as ``uhf`` even when
    # the requested ORCA calculation is DFT, where ORCA's keyword is ``UKS``.
    if str(scftype or qc_params.get("scftype", "rhf")).lower() in {"uhf", "uks"}:
        keyword += " UKS"
    return keyword


def _write_input(path, atoms, coordinates, charge, multiplicity, scftype, qc_params, scan=None):
    lines = [_orca_keyword(qc_params, scftype)]
    lines.append(f"%pal nprocs {int(qc_params.get('nprocs') or 1)} end")
    lines.append(f"%scf maxiter {int(qc_params.get('scf_cycles') or 1000)} end")
    if qc_params.get("opt_cycles") is not None or scan is not None:
        lines.append("%geom")
        if qc_params.get("opt_cycles") is not None:
            lines.append(f"    MaxIter {int(qc_params['opt_cycles'])}")
        if scan is not None:
            lines.extend([
                "    Scan",
                f"        B {scan.atom_i} {scan.atom_j} = "
                f"{scan.start_distance_angstrom:.10f}, {scan.end_distance_angstrom:.10f}, {scan.n_points}",
                "    end",
            ])
        lines.append("end")
    lines.append(f"*xyz {int(charge)} {int(multiplicity)}")
    for symbol, coordinate in zip(atoms, coordinates):
        lines.append(f"{symbol:<2} {coordinate[0]: .12f} {coordinate[1]: .12f} {coordinate[2]: .12f}")
    lines.append("*")
    Path(path).write_text("\n".join(lines) + "\n")


def load_orca_bond_scan_result(molecule, request, directory):
    """Recover and validate a completed scan from its retained ORCA files."""
    directory = Path(directory)
    input_path = directory / "scan.inp"
    output_path = directory / "scan.out"
    trajectory_path = directory / "scan.allxyz"
    if not output_path.exists() or "****ORCA TERMINATED NORMALLY****" not in output_path.read_text():
        return OrcaBondScanResult(False, None, None, None, input_path, output_path,
                                  "abnormal_termination")
    if not trajectory_path.exists():
        return OrcaBondScanResult(False, None, None, None, input_path, output_path,
                                  "trajectory_missing")
    try:
        frames = parse_xyz_trajectory(trajectory_path, molecule.number_of_atoms, molecule.atoms_list)
    except ValueError:
        return OrcaBondScanResult(False, trajectory_path, None, None, input_path, output_path,
                                  "trajectory_invalid")
    if len(frames) != request.n_points:
        return OrcaBondScanResult(False, trajectory_path, None, None, input_path, output_path,
                                  "trajectory_incomplete")
    profile_path = directory / "scan.relaxscanact.dat"
    if not profile_path.exists():
        return OrcaBondScanResult(False, trajectory_path, None, None, input_path, output_path,
                                  "scan_profile_missing")
    try:
        profile = parse_orca_scan_profile(profile_path)
    except (OSError, ValueError):
        return OrcaBondScanResult(False, trajectory_path, None, None, input_path, output_path,
                                  "scan_profile_invalid")
    if len(profile) != request.n_points:
        return OrcaBondScanResult(False, trajectory_path, None, None, input_path, output_path,
                                  "scan_profile_incomplete")
    expected_distances = np.linspace(
        request.start_distance_angstrom,
        request.end_distance_angstrom,
        request.n_points,
    )
    for expected_distance, point, (_, coordinates, _) in zip(expected_distances, profile, frames):
        if not math.isclose(
            point["target_distance_angstrom"], float(expected_distance),
            rel_tol=0.0, abs_tol=1e-3,
        ):
            return OrcaBondScanResult(False, trajectory_path, None, None, input_path, output_path,
                                      "scan_profile_grid_mismatch")
        distance = float(np.linalg.norm(coordinates[request.atom_i] - coordinates[request.atom_j]))
        if not math.isclose(distance, point["target_distance_angstrom"], rel_tol=0.0, abs_tol=1e-3):
            return OrcaBondScanResult(False, trajectory_path, None, None, input_path, output_path,
                                      "scan_profile_geometry_mismatch")
    final_coordinates = frames[-1][1]
    final_path = directory / "final_scan.xyz"
    write_xyz(molecule.atoms_list, final_coordinates, final_path, job_name="final_scan", precision=12)
    stable_trajectory = directory / "scan_trajectory.xyz"
    if stable_trajectory != trajectory_path:
        stable_trajectory.write_text(trajectory_path.read_text())
    return OrcaBondScanResult(True, stable_trajectory, final_path, final_coordinates,
                              input_path, output_path, "success", frames, profile)


def run_orca_bond_scan(molecule, request, directory, qc_params):
    """Run one ORCA relaxed scan and recover its final trajectory frame."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    input_path = directory / "scan.inp"
    output_path = directory / "scan.out"
    trajectory_path = directory / "scan.allxyz"
    _write_input(input_path, molecule.atoms_list, molecule.coordinates,
                 molecule.charge, molecule.multiplicity, molecule.scftype,
                 qc_params, request)
    try:
        executable = require_executable("orca", "ORCA")
    except FileNotFoundError:
        return OrcaBondScanResult(False, None, None, None, input_path, output_path, "orca_missing")
    old_cwd = Path.cwd()
    try:
        os.chdir(directory)
        status = run_command([executable, input_path.name], stdout_path=output_path.name, stderr_path=output_path.name)
    finally:
        os.chdir(old_cwd)
    if status != 0:
        return OrcaBondScanResult(False, None, None, None, input_path, output_path, "orca_failed")
    return load_orca_bond_scan_result(molecule, request, directory)
