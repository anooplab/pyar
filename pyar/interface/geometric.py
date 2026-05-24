"""geomeTRIC-backed geometry optimization with optional AFIR bias.

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
import tempfile
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
from ase.units import Bohr, Hartree

from pyar import interface
from pyar.afir import restraints
from pyar.afir.utils import resolve_gamma
from pyar.data.units import angstrom2bohr
from pyar.interface import SF, require_executable, write_xyz
from pyar.interface.xtb_utils import xtb_parallel_args

geometric_logger = logging.getLogger("pyar.geometric")

_GEOMETRIC_STATE_FILE = "pyar_geometric_state.json"


def _as_numpy_positions(atoms):
    """Return ASE positions as a contiguous float array."""
    return np.asarray(atoms.get_positions(), dtype=float)


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


def _backend_result_to_forces(energy_hartree, gradient_hartree_per_bohr):
    """Convert atomic-unit energy/gradient to ASE eV/forces units."""
    energy_ev = energy_hartree * Hartree
    forces_ev_per_angstrom = -gradient_hartree_per_bohr * Hartree / Bohr
    return energy_ev, forces_ev_per_angstrom


def _evaluate_xtb(atoms, qc_params):
    """Evaluate xTB energy and forces for the current coordinates."""
    xtb_executable = require_executable("xtb", "xTB")
    coordinates = _as_numpy_positions(atoms)

    with tempfile.TemporaryDirectory() as tmpdir:
        xyz_path = Path(tmpdir) / "input.xyz"
        interface.write_xyz(
            atoms.get_chemical_symbols(),
            coordinates,
            str(xyz_path),
            job_name="pyar_geometric_xtb",
        )
        command = [xtb_executable, str(xyz_path)]
        command.extend(xtb_parallel_args(qc_params))
        if qc_params.get("charge", 0) != 0:
            command.extend(["-chrg", str(qc_params["charge"])])
        multiplicity = int(qc_params.get("multiplicity", 1) or 1)
        if multiplicity != 1:
            command.extend(["-uhf", str(max(multiplicity - 1, 1))])
        scftype = qc_params.get("scftype", "rhf")
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
        return _backend_result_to_forces(energy_hartree, gradient_hartree_per_bohr)


def _evaluate_aimnet2(atoms, qc_params):
    """Evaluate AIMNet2 energy and forces for the current coordinates."""
    from pyar.AIMNet2.calculators.aimnet2ase import AIMNet2Calculator
    from pyar.interface.aimnet_2 import aimnet2

    calculator = AIMNet2Calculator(aimnet2, charge=int(qc_params.get("charge", 0) or 0))
    ase_atoms = Atoms(
        symbols=atoms.get_chemical_symbols(),
        positions=_as_numpy_positions(atoms),
    )
    ase_atoms.calc = calculator
    energy_ev = float(ase_atoms.get_potential_energy())
    forces_ev_per_angstrom = np.asarray(ase_atoms.get_forces(), dtype=float)
    return energy_ev, forces_ev_per_angstrom


def _resolve_backend_evaluator(software):
    """Return a callable that evaluates the selected backend."""
    if software == "xtb":
        return _evaluate_xtb
    if software == "aimnet_2":
        return _evaluate_aimnet2

    raise ValueError(
        "geomeTRIC optimization currently supports only 'xtb' and 'aimnet_2' "
        f"as backend energy/gradient providers, not {software!r}"
    )


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
        self._backend_evaluator = _resolve_backend_evaluator(self.software)

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

        backend_energy, backend_forces = self._backend_evaluator(self.atoms, self.qc_params)
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
            "pyar.interface.geometric.PyarGeometricCalculator",
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
