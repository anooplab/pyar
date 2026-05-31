"""Energy and gradient providers for geomeTRIC-compatible backends."""

from __future__ import annotations

import subprocess as subp
import tempfile
import re
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Protocol, runtime_checkable

import numpy as np
from ase import Atoms
from ase.units import Bohr, Hartree

from pyar import backends
from pyar.backend_capabilities import (
    backend_supports_geometry_optimization,
    normalize_backend_name,
    get_backend_capabilities,
    register_backend_capabilities,
    supported_geometry_backends,
)
from pyar.data.units import bohr2angstrom
from pyar.backends import require_executable
from pyar.backends.xtb_utils import xtb_parallel_args


@dataclass(frozen=True)
class EnergyGradientResult:
    """Backend energy and Cartesian gradient in atomic units."""

    energy_hartree: float
    gradient_hartree_per_bohr: np.ndarray

    def __post_init__(self):
        """Validate and normalize the returned energy/gradient pair."""
        try:
            energy = float(self.energy_hartree)
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"energy_hartree must be a finite scalar; got {self.energy_hartree!r}"
            ) from exc
        if not np.isfinite(energy):
            raise ValueError(
                f"energy_hartree must be a finite scalar; got {self.energy_hartree!r}"
            )

        try:
            gradient = np.asarray(self.gradient_hartree_per_bohr, dtype=float)
        except (TypeError, ValueError) as exc:
            raise ValueError(
                "gradient_hartree_per_bohr must be a finite array with shape (natoms, 3); "
                f"got invalid data {self.gradient_hartree_per_bohr!r}"
            ) from exc
        if gradient.ndim != 2 or gradient.shape[1] != 3:
            raise ValueError(
                "gradient_hartree_per_bohr must be a finite array with shape (natoms, 3); "
                f"got shape {gradient.shape!r}"
            )
        if not np.all(np.isfinite(gradient)):
            raise ValueError(
                "gradient_hartree_per_bohr must contain only finite values; "
                f"got array with shape {gradient.shape!r}"
            )

        object.__setattr__(self, "energy_hartree", energy)
        object.__setattr__(self, "gradient_hartree_per_bohr", gradient)


@runtime_checkable
class EnergyGradientProvider(Protocol):
    """Provide a backend energy and Cartesian gradient."""

    def evaluate(self, molecule, coordinates_bohr: np.ndarray) -> EnergyGradientResult:
        """Evaluate the backend objective in Bohr/Hartree units."""


def _as_numpy_positions(coordinates_bohr):
    """Return Cartesian positions as a contiguous float array in Bohr."""
    return np.asarray(coordinates_bohr, dtype=float)


def _atomic_numbers(molecule):
    """Return atomic numbers for either an ASE or PyAR molecule object."""
    if hasattr(molecule, "get_atomic_numbers"):
        return molecule.get_atomic_numbers()
    return molecule.atomic_number


def _read_xtb_energy(stdout_lines):
    """Extract the final xTB total energy from output text."""
    energy = None
    for line in stdout_lines:
        if "TOTAL ENERGY" in line:
            try:
                energy = float(line.split()[3])
            except (IndexError, ValueError):
                continue
        elif "total E" in line:
            try:
                energy = float(line.split()[-1])
            except (IndexError, ValueError):
                continue

    if energy is None:
        raise ValueError("Could not read xTB total energy from output")
    return energy


def _read_xtb_gradient(gradient_file):
    """Read an xTB gradient file into a Cartesian gradient array."""
    gradients = []
    with open(gradient_file) as fp:
        for line in fp:
            parts = line.split()
            if len(parts) != 3:
                continue
            gradients.append([float(xx) for xx in parts])
    if not gradients:
        raise ValueError(f"Could not read xTB gradients from {gradient_file}")
    return np.asarray(gradients, dtype=float)


def _forces_ev_per_angstrom_to_gradient_hartree_per_bohr(forces_ev_per_angstrom):
    """Convert ASE forces to an energy gradient in Hartree/Bohr."""
    return -np.asarray(forces_ev_per_angstrom, dtype=float) * Bohr / Hartree


class XtbEnergyGradientProvider:
    """Evaluate xTB energy and gradient for the current coordinates."""

    def __init__(self, qc_params=None):
        self.qc_params = dict(qc_params or {})

    def evaluate(self, molecule, coordinates_bohr):
        xtb_executable = require_executable("xtb", "xTB")
        coordinates_bohr = _as_numpy_positions(coordinates_bohr)
        coordinates_angstrom = bohr2angstrom(coordinates_bohr)

        with tempfile.TemporaryDirectory() as tmpdir:
            xyz_path = Path(tmpdir) / "input.xyz"
            backends.write_xyz(
                molecule.get_chemical_symbols(),
                coordinates_angstrom,
                str(xyz_path),
                job_name="pyar_geometric_xtb",
            )
            command = [xtb_executable, str(xyz_path)]
            command.extend(xtb_parallel_args(self.qc_params))
            if self.qc_params.get("charge", 0) != 0:
                command.extend(["-chrg", str(self.qc_params["charge"])])
            multiplicity = int(self.qc_params.get("multiplicity", 1) or 1)
            if multiplicity != 1:
                command.extend(["-uhf", str(max(multiplicity - 1, 1))])
            scftype = self.qc_params.get("scftype", "rhf")
            if multiplicity == 1 and scftype != "rhf":
                command.append(f"-{scftype}")
            command.append("--grad")
            proc = subp.run(
                command,
                cwd=tmpdir,
                capture_output=True,
                text=True,
                check=False,
            )
            if proc.returncode != 0:
                raise RuntimeError(
                    "xTB energy/gradient evaluation failed: "
                    + (proc.stderr.strip() or proc.stdout.strip() or "unknown error")
                )

            gradient_path = Path(tmpdir) / "gradient"
            if not gradient_path.exists():
                raise RuntimeError("xTB gradient file was not written")

            energy_hartree = _read_xtb_energy(proc.stdout.splitlines())
            gradient_hartree_per_bohr = _read_xtb_gradient(gradient_path)
            return EnergyGradientResult(energy_hartree, gradient_hartree_per_bohr)


class Aimnet2EnergyGradientProvider:
    """Evaluate AIMNet2 energy and gradient for the current coordinates."""

    def __init__(self, qc_params=None):
        self.qc_params = dict(qc_params or {})

    def evaluate(self, molecule, coordinates_bohr):
        from pyar.AIMNet2.calculators.aimnet2ase import AIMNet2Calculator
        from pyar.backends.aimnet_2 import load_default_aimnet2_model

        coordinates_bohr = _as_numpy_positions(coordinates_bohr)
        coordinates_angstrom = bohr2angstrom(coordinates_bohr)
        calculator = AIMNet2Calculator(
            load_default_aimnet2_model(),
            charge=int(self.qc_params.get("charge", 0) or 0),
        )
        ase_atoms = Atoms(
            symbols=molecule.get_chemical_symbols(),
            positions=coordinates_angstrom,
        )
        ase_atoms.calc = calculator
        energy_ev = float(ase_atoms.get_potential_energy())
        forces_ev_per_angstrom = np.asarray(ase_atoms.get_forces(), dtype=float)
        return EnergyGradientResult(
            energy_ev / Hartree,
            _forces_ev_per_angstrom_to_gradient_hartree_per_bohr(forces_ev_per_angstrom),
        )


def _read_gaussian_energy(output_lines):
    """Extract the final SCF energy from Gaussian output."""
    for line in reversed(output_lines):
        if "SCF Done" in line:
            try:
                return float(line.split()[4])
            except (IndexError, ValueError) as exc:
                raise ValueError("Could not parse Gaussian energy from output") from exc
    raise ValueError("Could not read Gaussian energy from output")


def _read_gaussian_forces(output_lines, number_of_atoms):
    """Extract Gaussian forces from the gradient block and convert to gradients."""
    for index, line in enumerate(output_lines):
        if "Forces (Hartrees/Bohr)" in line:
            forces = []
            for row in output_lines[index + 1:]:
                parts = row.split()
                if len(parts) < 5 or not parts[0].isdigit():
                    continue
                forces.append([float(parts[2]), float(parts[3]), float(parts[4])])
                if len(forces) == number_of_atoms:
                    break
            if len(forces) != number_of_atoms:
                break
            return -np.asarray(forces, dtype=float)
    raise ValueError("Could not read Gaussian forces from output")


class GaussianEnergyGradientProvider:
    """Evaluate Gaussian energy and gradient for the current coordinates."""

    def __init__(self, qc_params=None):
        self.qc_params = dict(qc_params or {})

    def _build_keyword(self, molecule):
        """Build a Gaussian input route section for a force evaluation."""
        method = self.qc_params["method"]
        basis = self.qc_params["basis"]
        if basis.lower() == "def2-svp":
            basis = "def2SVP"
        keyword = f"%nprocshared={int(self.qc_params.get('nprocs', 1) or 1)}\n"
        keyword += f"%chk=trial_gaussian.chk\n"
        keyword += f"%mem=2GB\n"
        keyword += f"# {method} {basis} Force SCF=(MaxCycle={int(self.qc_params.get('scf_cycles', 350) or 350)})"
        if int(self.qc_params.get("multiplicity", 1) or 1) != 1:
            keyword += " guess=mix"
        return keyword

    def _write_input(self, molecule, coordinates_angstrom, inp_path):
        """Write a Gaussian force input file for the current geometry."""
        with open(inp_path, "w") as fp:
            fp.write(self._build_keyword(molecule) + "\n\n")
            fp.write("PyAR Gaussian energy/gradient provider\n\n")
            fp.write(f"{int(self.qc_params.get('charge', 0) or 0)} {int(self.qc_params.get('multiplicity', 1) or 1)}\n")
            for symbol, coord in zip(molecule.get_chemical_symbols(), coordinates_angstrom):
                fp.write(f" {symbol:>3}  {coord[0]:10.7f}  {coord[1]:10.7f} {coord[2]:10.7f}\n")
            fp.write("\n")

    def evaluate(self, molecule, coordinates_bohr):
        """Run one Gaussian force job and parse the resulting gradient."""
        gaussian_executable = require_executable("g16", "Gaussian")
        coordinates_bohr = _as_numpy_positions(coordinates_bohr)
        coordinates_angstrom = bohr2angstrom(coordinates_bohr)

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            inp_path = tmpdir / "input.com"

            self._write_input(molecule, coordinates_angstrom, inp_path)
            proc = subp.run(
                [gaussian_executable, str(inp_path)],
                cwd=tmpdir,
                capture_output=True,
                text=True,
                check=False,
            )
            if proc.returncode != 0:
                raise RuntimeError(
                    "Gaussian energy/gradient evaluation failed: "
                    + (proc.stderr.strip() or proc.stdout.strip() or "unknown error")
                )
            output_lines = (proc.stdout or "").splitlines()
            energy_hartree = _read_gaussian_energy(output_lines)
            gradient_hartree_per_bohr = _read_gaussian_forces(output_lines, len(molecule.get_chemical_symbols()))
            return EnergyGradientResult(energy_hartree, gradient_hartree_per_bohr)


def _read_orca_energy(output_text):
    """Extract the final ORCA single-point energy from the output text."""
    lines = output_text.splitlines()
    for line in reversed(lines):
        if "FINAL SINGLE POINT ENERGY" in line:
            try:
                return float(line.split()[-1])
            except (IndexError, ValueError) as exc:
                raise ValueError("Could not parse ORCA energy from ORCA output") from exc

    raise ValueError("Could not read ORCA single-point energy from ORCA output")


def _read_orca_engrad(engrad_path, number_of_atoms):
    """Read an ORCA `.engrad` file into an energy and Cartesian gradient."""
    text = Path(engrad_path).read_text()
    tokens = re.findall(r"[-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?", text)
    if len(tokens) < 2 + 3 * number_of_atoms:
        raise ValueError(f"Could not read ORCA gradients from {engrad_path}")

    atoms = int(float(tokens[0]))
    if atoms != number_of_atoms:
        raise ValueError(
            f"ORCA gradient file {engrad_path} reported {atoms} atoms, expected {number_of_atoms}"
        )

    energy_hartree = float(tokens[1])
    gradient_values = np.asarray([float(xx) for xx in tokens[2: 2 + 3 * atoms]], dtype=float)
    gradient_hartree_per_bohr = gradient_values.reshape(atoms, 3)
    return energy_hartree, gradient_hartree_per_bohr


class OrcaEnergyGradientProvider:
    """Evaluate ORCA energy and gradient for the current coordinates."""

    def __init__(self, qc_params=None):
        self.qc_params = dict(qc_params or {})

    def _build_keyword(self, molecule):
        """Build an ORCA input keyword block for a gradient evaluation."""
        method = self.qc_params["method"]
        basis = self.qc_params["basis"]

        keyword = f"! {method} {basis}"
        if any(x >= 21 for x in _atomic_numbers(molecule)):
            keyword += " def2-ECP"
        keyword += " RI def2/J D3BJ KDIIS ENGRAD"
        if int(self.qc_params.get("multiplicity", 1) or 1) != 1:
            keyword += " UKS"
        keyword += f"\n%pal nprocs {int(self.qc_params.get('nprocs', 1) or 1)} end\n"
        keyword += f"%scf maxiter {int(self.qc_params.get('scf_cycles', 350) or 350)} end\n"
        return keyword

    def _write_input(self, molecule, coordinates_angstrom, inp_path):
        """Write a minimal ORCA input file for one energy/gradient calculation."""
        with open(inp_path, "w") as fp:
            fp.write(self._build_keyword(molecule) + "\n")
            fp.write(
                f"*xyz {int(self.qc_params.get('charge', 0) or 0)} "
                f"{int(self.qc_params.get('multiplicity', 1) or 1)}\n"
            )
            for symbol, coord in zip(molecule.get_chemical_symbols(), coordinates_angstrom):
                fp.write(f" {symbol:>3}  {coord[0]:10.7f}  {coord[1]:10.7f} {coord[2]:10.7f}\n")
            fp.write("*\n")

    def evaluate(self, molecule, coordinates_bohr):
        """Run one ORCA single-point calculation and parse its gradient."""
        orca_executable = require_executable("orca", "ORCA")
        coordinates_bohr = _as_numpy_positions(coordinates_bohr)
        coordinates_angstrom = bohr2angstrom(coordinates_bohr)

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            inp_path = tmpdir / "input.inp"
            engrad_path = tmpdir / "input.engrad"

            self._write_input(molecule, coordinates_angstrom, inp_path)
            proc = subp.run(
                [orca_executable, str(inp_path)],
                cwd=tmpdir,
                capture_output=True,
                text=True,
                check=False,
            )
            if proc.returncode != 0:
                raise RuntimeError(
                    "ORCA energy/gradient evaluation failed: "
                    + (proc.stderr.strip() or proc.stdout.strip() or "unknown error")
                )
            if not engrad_path.exists():
                raise RuntimeError("ORCA gradient file was not written")

            try:
                energy_hartree = _read_orca_energy(proc.stdout or "")
            except ValueError:
                energy_hartree, _ = _read_orca_engrad(engrad_path, len(molecule.get_chemical_symbols()))

            _, gradient_hartree_per_bohr = _read_orca_engrad(
                engrad_path,
                len(molecule.get_chemical_symbols()),
            )
            return EnergyGradientResult(energy_hartree, gradient_hartree_per_bohr)


ENERGY_GRADIENT_PROVIDERS = {
    "xtb": XtbEnergyGradientProvider,
    "aimnet_2": Aimnet2EnergyGradientProvider,
    "gaussian": GaussianEnergyGradientProvider,
    "orca": OrcaEnergyGradientProvider,
}


def get_energy_gradient_provider(software, qc_params=None):
    """Return the registered energy/gradient provider for a backend."""
    canonical_software = normalize_backend_name(software)
    provider_factory = ENERGY_GRADIENT_PROVIDERS.get(canonical_software)
    if provider_factory is not None and backend_supports_geometry_optimization(canonical_software):
        if hasattr(provider_factory, "evaluate") and not isinstance(provider_factory, type):
            return provider_factory
        return provider_factory(qc_params)

    if get_backend_capabilities(canonical_software).energy_gradient:
        raise ValueError(
            f"Backend {software!r} is marked as energy-gradient capable but has no registered provider"
        )

    raise ValueError(
        "geomeTRIC optimization currently supports only "
        f"{', '.join(supported_geometry_backends())} as backend energy/gradient providers, "
        f"not {software!r}"
    )


def register_energy_gradient_provider(software, provider, capabilities=None):
    """Register or replace an energy/gradient provider for a backend."""
    if not callable(provider) and not hasattr(provider, "evaluate"):
        raise TypeError("provider must be callable or expose an evaluate() method")
    canonical_software = normalize_backend_name(software)
    ENERGY_GRADIENT_PROVIDERS[canonical_software] = provider
    if capabilities is None:
        capabilities = get_backend_capabilities(canonical_software)
        if not capabilities.energy_gradient:
            capabilities = replace(capabilities, energy_gradient=True)
    register_backend_capabilities(canonical_software, capabilities)
