"""geomeTRIC-backed optimization backend with optional AFIR bias.

This module provides a small bridge between PyAR's backend selection and
geomeTRIC's internal-coordinate optimizer.  The optimizer itself is backend
agnostic: it asks a selected PyAR backend for energy and gradients, then adds
the AFIR term when ``gamma`` is non-zero.
"""

from __future__ import annotations

import json
import logging
import subprocess as subp
import sys
from pathlib import Path

import numpy as np
from ase.calculators.calculator import Calculator, all_changes
from ase.units import Bohr, Hartree

from pyar.afir import restraints
from pyar.afir.utils import resolve_gamma
from pyar.energy_gradient_providers import EnergyGradientResult, get_energy_gradient_provider
from pyar.data.units import angstrom2bohr
from pyar.backends import SF, require_executable, write_xyz

geometric_logger = logging.getLogger("pyar.geometric")

_GEOMETRIC_STATE_FILE = "pyar_geometric_state.json"


def _find_geometric_executable():
    """Find geomeTRIC from the active Python environment or ``PATH``."""
    environment_script = Path(sys.executable).with_name("geometric-optimize")
    if environment_script.is_file():
        return str(environment_script)
    return require_executable("geometric-optimize", "geomeTRIC")


def _read_last_xyz(path):
    """Read the last XYZ frame from a geomeTRIC output file."""
    lines = Path(path).read_text().splitlines()
    if len(lines) < 2:
        raise ValueError(f"Invalid XYZ file: {path}")

    natoms = int(lines[0].strip())
    start = len(lines) - (natoms + 2)
    if start < 0:
        raise ValueError(f"Invalid XYZ file: {path}")

    frame = lines[start:]
    coords = np.loadtxt(frame[2:], usecols=(1, 2, 3), dtype=float)
    return coords


def _resolve_backend_evaluator(software, qc_params):
    """Return the registered backend energy/gradient provider."""
    return get_energy_gradient_provider(software, qc_params)


class PyarGeometricCalculator(Calculator):
    """ASE calculator that combines a PyAR backend with optional AFIR bias."""

    implemented_properties = ["energy", "forces"]

    def __init__(self, qc_params, fragment_indices=None, opt_target="minimum"):
        super().__init__()
        self.qc_params = dict(qc_params or {})
        self.software = self.qc_params.get("software")
        self.gamma = resolve_gamma(self.qc_params.get("gamma"), fallback=0.0)
        self.fragment_indices = fragment_indices
        self.opt_target = opt_target
        self._backend_evaluator = _resolve_backend_evaluator(self.software, self.qc_params)

    @staticmethod
    def _result_to_ase_units(result):
        """Convert backend Hartree/Bohr data to ASE energy/force units."""
        if not isinstance(result, EnergyGradientResult):
            raise TypeError("backend provider must return EnergyGradientResult")
        backend_energy_ev = float(result.energy_hartree) * Hartree
        backend_forces_ev_per_angstrom = (
            -np.asarray(result.gradient_hartree_per_bohr, dtype=float) * Hartree / Bohr
        )
        return backend_energy_ev, backend_forces_ev_per_angstrom

    def _afir_contribution(self):
        """Return AFIR energy and forces in ASE units."""
        if self.gamma == 0.0:
            natoms = len(self.atoms)
            return 0.0, np.zeros((natoms, 3), dtype=float)

        if not self.fragment_indices or len(self.fragment_indices) != 2:
            raise ValueError(
                "AFIR bias requires exactly two fragment index lists in the geometric optimizer"
            )

        coordinates_bohr = angstrom2bohr(np.asarray(self.atoms.get_positions(), dtype=float))
        afir_energy_hartree, afir_force_hartree_per_bohr = restraints.isotropic(
            self.fragment_indices,
            list(self.atoms.get_chemical_symbols()),
            coordinates_bohr,
            self.gamma,
        )
        afir_energy_ev = afir_energy_hartree * Hartree
        # ``restraints.isotropic`` returns -dE/dx, i.e. force, in atomic units.
        afir_forces_ev_per_angstrom = afir_force_hartree_per_bohr * Hartree / Bohr
        return afir_energy_ev, afir_forces_ev_per_angstrom

    def _write_state(self, energy, forces):
        """Persist the latest evaluation so the parent process can recover it."""
        state = {
            "software": self.software,
            "gamma": self.gamma,
            "opt_target": self.opt_target,
            "energy_ev": float(energy),
            "energy_hartree": float(energy / Hartree),
            "positions_angstrom": np.asarray(self.atoms.get_positions(), dtype=float).tolist(),
            "forces_ev_per_angstrom": np.asarray(forces, dtype=float).tolist(),
        }
        with open(_GEOMETRIC_STATE_FILE, "w") as fp:
            json.dump(state, fp, indent=2, sort_keys=True)

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        """Compute the backend objective plus the AFIR bias if enabled."""
        super().calculate(atoms, properties, system_changes)

        coordinates_bohr = angstrom2bohr(np.asarray(self.atoms.get_positions(), dtype=float))
        backend_result = self._backend_evaluator.evaluate(self.atoms, coordinates_bohr)
        backend_energy, backend_forces = self._result_to_ase_units(backend_result)
        afir_energy, afir_forces = self._afir_contribution()
        total_energy = backend_energy + afir_energy
        total_forces = np.asarray(backend_forces, dtype=float) + np.asarray(afir_forces, dtype=float)

        self.results["energy"] = float(total_energy)
        self.results["forces"] = total_forces
        self._write_state(total_energy, total_forces)


class Geometric(SF):
    """Run geomeTRIC with a PyAR backend calculator and optional AFIR bias."""

    def __init__(self, molecule, qc_params):
        super().__init__(molecule)
        self.qc_params = dict(qc_params or {})
        self.qc_params.update(
            charge=self.charge,
            multiplicity=self.multiplicity,
            scftype=self.scftype,
        )
        self.software = self.qc_params.get("software")
        self.gamma = resolve_gamma(self.qc_params.get("gamma"), fallback=0.0)
        self.opt_target = self.qc_params.get("opt_target", "minimum")
        self.fragment_indices = molecule.fragments
        self.geometric_executable = _find_geometric_executable()

    def _build_command(self):
        """Build the geomeTRIC command line."""
        if self.opt_target not in {"minimum", "ts"}:
            raise ValueError(f"Unsupported geomeTRIC optimization target: {self.opt_target!r}")
        if self.opt_target == "ts":
            raise NotImplementedError(
                "Transition-state optimization is reserved for a future reaction-product workflow"
            )

        if self.gamma != 0.0 and (not self.fragment_indices or len(self.fragment_indices) != 2):
            raise ValueError("AFIR geometry optimization requires exactly two fragments")

        ase_kwargs = {
            "qc_params": self.qc_params,
            "fragment_indices": self.fragment_indices,
            "opt_target": self.opt_target,
        }
        command = [
            self.geometric_executable,
            "--engine",
            "ase",
            "--ase-class",
            "pyar.backends.geometric.PyarGeometricCalculator",
            "--ase-kwargs",
            json.dumps(ase_kwargs),
            self.start_xyz_file,
        ]
        command.extend(["--coordsys", "tric"])
        if self.qc_params.get("opt_cycles") is not None:
            command.extend(["--maxiter", str(int(self.qc_params["opt_cycles"]))])
        convergence = self.qc_params.get("opt_threshold")
        if convergence is not None:
            geometric_threshold_map = {
                "loose": "GAU_LOOSE",
                "normal": "GAU",
                "tight": "GAU_TIGHT",
            }
            command.extend(["--converge", "set", geometric_threshold_map.get(convergence, "GAU")])
        return command

    def _read_final_energy(self):
        """Recover the latest energy from the calculator state file."""
        state_path = Path(_GEOMETRIC_STATE_FILE)
        if not state_path.exists():
            return None
        with state_path.open() as fp:
            state = json.load(fp)
        return state.get("energy_hartree")

    def _read_final_xyz(self):
        """Find the final geometry written by geomeTRIC."""
        stem = Path(self.start_xyz_file).stem
        candidates = [
            Path(f"{stem}_optim.xyz"),
            Path("opt.xyz"),
        ]
        for candidate in candidates:
            if candidate.exists():
                return _read_last_xyz(candidate)

        xyz_candidates = sorted(
            (
                path
                for path in Path(".").glob("*.xyz")
                if path.name != self.start_xyz_file
            ),
            key=lambda path: path.stat().st_mtime,
        )
        if xyz_candidates:
            return _read_last_xyz(xyz_candidates[-1])
        raise FileNotFoundError("geomeTRIC did not write a final XYZ file")

    def optimize(self):
        """Run the geomeTRIC optimization loop."""
        state_path = Path(_GEOMETRIC_STATE_FILE)
        state_path.unlink(missing_ok=True)
        for stale_path in (
            Path(f"{Path(self.start_xyz_file).stem}_optim.xyz"),
            Path("opt.xyz"),
        ):
            stale_path.unlink(missing_ok=True)

        command = self._build_command()
        geometric_logger.info(
            "geomeTRIC start: name=%s software=%s gamma=%s target=%s",
            self.job_name,
            self.software,
            self.gamma,
            self.opt_target,
        )

        with open("geometric.out", "w") as output_file_pointer:
            proc = subp.run(
                command,
                stdout=output_file_pointer,
                stderr=output_file_pointer,
                text=True,
                check=False,
            )

        if proc.returncode != 0:
            geometric_logger.error(
                "geomeTRIC failed: name=%s software=%s returncode=%s",
                self.job_name,
                self.software,
                proc.returncode,
            )
            return False

        try:
            self.optimized_coordinates = self._read_final_xyz()
            self.coordinates = self.optimized_coordinates
            self.energy = self._read_final_energy()
            if self.energy is None:
                raise FileNotFoundError(_GEOMETRIC_STATE_FILE)
        except Exception as exc:
            geometric_logger.error("geomeTRIC completed but final state was incomplete: %s", exc)
            return False

        write_xyz(
            self.atoms_list,
            self.optimized_coordinates,
            self.result_xyz_file,
            job_name=self.job_name,
            energy=self.energy,
        )
        geometric_logger.info(
            "geomeTRIC completed: name=%s energy=%15.6f",
            self.job_name,
            float(self.energy),
        )
        return True


def main():
    """Module entry point for direct execution."""
    pass
