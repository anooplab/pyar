"""geomeTRIC-backed optimization backend with an optional reaction bias.

This module provides a small bridge between PyAR's backend selection and
geomeTRIC's internal-coordinate optimizer.  The optimizer itself is backend
agnostic: it asks a selected PyAR backend for energy and gradients, then adds
the selected bias term when ``gamma`` is non-zero.
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

from pyar.biases import afir as restraints, softmin
from pyar.biases.afir import resolve_gamma
from pyar.biases.collective_coordinates import evaluate_contact_coordinate
from pyar.energy_gradient_providers import EnergyGradientResult, get_energy_gradient_provider
from pyar.data.units import angstrom2bohr
from pyar.backends import SF, require_executable, write_xyz
from pyar.reaction_trace import ReactionTraceRecorder

geometric_logger = logging.getLogger("pyar.geometric")

_GEOMETRIC_STATE_FILE = "pyar_geometric_state.json"
_BIAS_POTENTIALS = {
    "afir",
    "softmin",
}


def _resolve_bias_potential(value):
    """Return the selected reaction-bias name."""
    name = "afir" if value is None else str(value).lower()
    if name not in _BIAS_POTENTIALS:
        choices = ", ".join(sorted(_BIAS_POTENTIALS))
        raise ValueError(f"Unsupported bias potential: {value!r}; choose one of {choices}") from None
    return name


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
    """ASE calculator that combines a PyAR backend with an optional reaction bias."""

    implemented_properties = ["energy", "forces"]

    def __init__(self, qc_params, fragment_indices=None, opt_target="minimum"):
        super().__init__()
        self.qc_params = dict(qc_params or {})
        self.software = self.qc_params.get("software")
        self.gamma = resolve_gamma(self.qc_params.get("gamma"), fallback=0.0)
        self.bias_potential = _resolve_bias_potential(self.qc_params.get("bias_potential"))
        self.softmin_beta = softmin.resolve_softmin_beta(self.qc_params.get("softmin_beta"))
        self.fragment_indices = fragment_indices
        self.opt_target = opt_target
        self._backend_evaluator = _resolve_backend_evaluator(self.software, self.qc_params)
        self.trace_enabled = bool(
            self.qc_params.get("trace_enabled") or self.qc_params.get("reaction_trace")
        )
        self._trace_recorder = None
        if self.trace_enabled:
            trace_name = self.qc_params.get("trace_name", "reaction_trace")
            trace_root = Path.cwd() / trace_name
            trace_mode = self.qc_params.get("trace_mode")
            if trace_mode is None:
                trace_mode = "append" if (trace_root / "trace.jsonl").exists() else "write"
            self._trace_recorder = ReactionTraceRecorder(
                Path.cwd(),
                trace_name=trace_name,
                mode=trace_mode,
            )

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

    def _bias_contribution(self):
        """Return selected bias energy and forces in atomic and ASE units."""
        if self.gamma == 0.0:
            natoms = len(self.atoms)
            zero_forces = np.zeros((natoms, 3), dtype=float)
            return 0.0, zero_forces, 0.0, zero_forces

        if not self.fragment_indices or len(self.fragment_indices) != 2:
            raise ValueError(
                "Reaction bias requires exactly two fragment index lists in the geometric optimizer"
            )

        coordinates_bohr = angstrom2bohr(np.asarray(self.atoms.get_positions(), dtype=float))
        evaluator = restraints.isotropic if self.bias_potential == "afir" else softmin.softmin
        bias_arguments = (
            self.fragment_indices,
            list(self.atoms.get_chemical_symbols()),
            coordinates_bohr,
            self.gamma,
        )
        if self.bias_potential == "softmin":
            bias_energy_hartree, bias_force_hartree_per_bohr = evaluator(
                *bias_arguments, beta=self.softmin_beta
            )
        else:
            bias_energy_hartree, bias_force_hartree_per_bohr = evaluator(*bias_arguments)
        bias_forces_hartree_per_bohr = np.asarray(bias_force_hartree_per_bohr, dtype=float)
        bias_energy_ev = bias_energy_hartree * Hartree
        bias_forces_ev_per_angstrom = bias_forces_hartree_per_bohr * Hartree / Bohr
        return (
            float(bias_energy_hartree),
            bias_forces_hartree_per_bohr,
            bias_energy_ev,
            bias_forces_ev_per_angstrom,
        )

    def _contact_coordinate_report(self):
        """Return the current bias coordinate and JSON-ready contact diagnostics."""
        if not self.fragment_indices or len(self.fragment_indices) != 2:
            return None, None
        coordinates_bohr = angstrom2bohr(np.asarray(self.atoms.get_positions(), dtype=float))
        coordinate_kwargs = {"kind": self.bias_potential}
        if self.bias_potential == "softmin":
            coordinate_kwargs["beta"] = self.softmin_beta
        q, _gradient, diagnostics = evaluate_contact_coordinate(
            self.fragment_indices,
            list(self.atoms.get_chemical_symbols()),
            coordinates_bohr,
            **coordinate_kwargs,
        )
        return float(q), diagnostics.as_dict()

    def _write_state(self, energy, forces, collective_coordinate_bohr=None, contact_diagnostics=None):
        """Persist the latest evaluation so the parent process can recover it."""
        state = {
            "software": self.software,
            "gamma": self.gamma,
            "bias_potential": self.bias_potential,
            "softmin_beta": self.softmin_beta,
            "opt_target": self.opt_target,
            "energy_ev": float(energy),
            "energy_hartree": float(energy / Hartree),
            "positions_angstrom": np.asarray(self.atoms.get_positions(), dtype=float).tolist(),
            "forces_ev_per_angstrom": np.asarray(forces, dtype=float).tolist(),
        }
        if collective_coordinate_bohr is not None:
            state["collective_coordinate_bohr"] = float(collective_coordinate_bohr)
        if contact_diagnostics is not None:
            state["contact_diagnostics"] = contact_diagnostics
        with open(_GEOMETRIC_STATE_FILE, "w") as fp:
            json.dump(state, fp, indent=2, sort_keys=True)

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        """Compute the backend objective plus the selected bias if enabled."""
        super().calculate(atoms, properties, system_changes)

        coordinates_bohr = angstrom2bohr(np.asarray(self.atoms.get_positions(), dtype=float))
        backend_result = self._backend_evaluator.evaluate(self.atoms, coordinates_bohr)
        backend_energy, backend_forces = self._result_to_ase_units(backend_result)
        backend_energy_hartree = float(backend_result.energy_hartree)
        # Providers return the energy gradient; traces and force diagnostics
        # must use the physical force on the same sign convention as ASE.
        backend_forces_hartree_per_bohr = -np.asarray(
            backend_result.gradient_hartree_per_bohr, dtype=float
        )

        (
            bias_energy_hartree,
            bias_forces_hartree_per_bohr,
            bias_energy,
            bias_forces,
        ) = self._bias_contribution()
        total_energy = backend_energy + bias_energy
        total_forces = np.asarray(backend_forces, dtype=float) + np.asarray(bias_forces, dtype=float)
        total_energy_hartree = backend_energy_hartree + bias_energy_hartree
        total_forces_hartree_per_bohr = backend_forces_hartree_per_bohr + bias_forces_hartree_per_bohr

        backend_force_norm = float(np.linalg.norm(backend_forces_hartree_per_bohr))
        bias_force_norm = float(np.linalg.norm(bias_forces_hartree_per_bohr))
        total_force_norm = float(np.linalg.norm(total_forces_hartree_per_bohr))
        max_force = float(np.max(np.linalg.norm(total_forces_hartree_per_bohr, axis=1)))
        collective_coordinate_bohr, contact_diagnostics = self._contact_coordinate_report()

        self.results["energy"] = float(total_energy)
        self.results["forces"] = total_forces
        self._write_state(
            total_energy,
            total_forces,
            collective_coordinate_bohr=collective_coordinate_bohr,
            contact_diagnostics=contact_diagnostics,
        )

        if self._trace_recorder is not None:
            self._trace_recorder.record(
                symbols=self.atoms.get_chemical_symbols(),
                coordinates_angstrom=np.asarray(self.atoms.get_positions(), dtype=float),
                backend_energy_hartree=backend_energy_hartree,
                bias_energy_hartree=bias_energy_hartree,
                total_energy_hartree=total_energy_hartree,
                backend_forces_hartree_per_bohr=backend_forces_hartree_per_bohr,
                bias_forces_hartree_per_bohr=bias_forces_hartree_per_bohr,
                total_forces_hartree_per_bohr=total_forces_hartree_per_bohr,
                backend_force_norm=backend_force_norm,
                bias_force_norm=bias_force_norm,
                total_force_norm=total_force_norm,
                max_force=max_force,
                fragment_indices=self.fragment_indices,
                collective_coordinate_bohr=collective_coordinate_bohr,
                contact_diagnostics=contact_diagnostics,
                softmin_beta=self.softmin_beta if self.bias_potential == "softmin" else None,
            )


class Geometric(SF):
    """Run geomeTRIC with a selected reaction bias."""

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
        self.bias_potential = _resolve_bias_potential(self.qc_params.get("bias_potential"))
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
            raise ValueError("Reaction-bias geometry optimization requires exactly two fragments")

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
            "geomeTRIC start: name=%s software=%s bias=%s gamma=%s target=%s",
            self.job_name,
            self.software,
            self.bias_potential,
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
