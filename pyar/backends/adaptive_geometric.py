"""Accepted-step geomeTRIC driver for piecewise-constant adaptive bias.

The ordinary ASE calculator interface has no acceptance notification. This
driver owns that boundary; trial steps and Hessian probes use the current
segment unchanged. Fixed-bias runs continue to use geomeTRIC's standard CLI.
"""

from __future__ import annotations

import json
import logging
import sys
from pathlib import Path

import numpy as np
from ase.units import Bohr
from geometric.ase_engine import EngineASE
from geometric.internal import DelocalizedInternalCoordinates
from geometric.molecule import Molecule
from geometric.errors import GeomOptNotConvergedError
from geometric.optimize import OPT_STATE, Optimizer
from geometric.params import OptParams

from pyar.backends.geometric import PyarGeometricCalculator, _CONTROLLER_STATE_FILE


class ReleaseCandidateReached(Exception):
    """Stop at a fully accepted geometry that meets release criteria."""


class AdaptiveOptimizer(Optimizer):
    """Refresh the objective only after geomeTRIC has accepted a trial."""

    def evaluateStep(self):
        trial_coordinates = self.X.copy()
        super().evaluateStep()
        if self.state == OPT_STATE.FAILED:
            return
        # geomeTRIC restores Xprev on rejection. Do not infer acceptance from
        # evaluation count: cached trials and rejected steps are normal.
        if not np.array_equal(self.X, trial_coordinates):
            return

        self.engine.update_atoms(self.X)
        changed = self.engine.calculator.accept_geometry(
            self.engine.ase_atoms, advance_bias=False
        )
        self.engine.clearCalcs()
        result = self.engine.calc(self.X, self.dirname)
        self.E = result["energy"]
        self.gradx = result["gradient"]
        self.G = self.IC.calcGrad(self.X, self.gradx).flatten()
        if changed:
            # Secant history from different bias strengths is not a Hessian
            # of the new objective. Start its model at this accepted geometry.
            self.H0 = self.IC.guess_hessian(self.X)
            self.H = self.H0.copy()
            self.X_hist = [self.X.copy()]
            self.Gx_hist = [self.gradx.copy()]
            if hasattr(self, "X_rj"):
                del self.X_rj
            if self.state == OPT_STATE.CONVERGED:
                rms_gradient, max_gradient = self.calcGradNorm()
                if (rms_gradient >= self.params.Convergence_grms or
                        max_gradient >= self.params.Convergence_gmax):
                    self.state = OPT_STATE.NEEDS_EVALUATION
        self.progress.qm_energies[-1] = self.E
        self.progress.qm_grads[-1] = self.gradx.copy()
        tracker = self.engine.calculator.release_tracker
        suppress_candidate = self.engine.calculator.qc_params.get(
            "adaptive_suppress_release_candidate", False
        )
        if (tracker is not None and tracker.evidence.state == "CANDIDATE"
                and not suppress_candidate):
            raise ReleaseCandidateReached

    def optimizeGeometry(self):
        return super().optimizeGeometry()


def run_adaptive_optimization(input_xyz, calculator_arguments):
    """Run minimum optimization with accepted-step updates and fresh Hessians.

    Explicit controller restart resumes the saved accepted geometry and bias
    history, but intentionally rebuilds geomeTRIC's trust/Hessian model.
    """
    molecule = Molecule(str(input_xyz))
    qc_params = calculator_arguments["qc_params"]
    if calculator_arguments.get("opt_target", "minimum") != "minimum":
        raise ValueError("Adaptive geomeTRIC currently supports minima only")
    if qc_params.get("bias_controller_restart"):
        checkpoint = json.loads(Path(_CONTROLLER_STATE_FILE).read_text())
        if molecule.elem != checkpoint["symbols"]:
            raise ValueError("Restart atom symbols do not match the input molecule")
        molecule.xyzs = [np.asarray(checkpoint["positions_angstrom"], dtype=float)]
        molecule.build_topology()
    calculator = PyarGeometricCalculator(**calculator_arguments)
    engine = EngineASE(molecule, calculator)
    # Use EngineASE's conversion constant, including on checkpoint reload.
    coordinates = molecule.xyzs[0].flatten() / Bohr
    threshold = {"loose": "GAU_LOOSE", "normal": "GAU", "tight": "GAU_TIGHT"}
    output_xyz = f"{Path(input_xyz).stem}_optim.xyz"
    params = OptParams(
        maxiter=qc_params.get("opt_cycles") if qc_params.get("opt_cycles") is not None else 300,
        convergence_set=threshold.get(qc_params.get("opt_threshold"), "GAU"),
        xyzout=output_xyz,
    )
    directory_root = Path(f"{Path(input_xyz).stem}.tmp")
    directory_root.mkdir(exist_ok=True)
    max_segments = int(qc_params.get("adaptive_max_segments", 64))
    if max_segments < 1:
        raise ValueError("adaptive_max_segments must be a positive integer")

    for segment in range(max_segments):
        molecule.xyzs = [coordinates.reshape((-1, 3)).copy() * Bohr]
        molecule.build_topology()
        internal = DelocalizedInternalCoordinates(
            molecule, build=True, connect=False, addcart=False
        )
        directory = str(directory_root / f"segment_{segment:03d}")
        Path(directory).mkdir(exist_ok=True)
        optimizer = AdaptiveOptimizer(coordinates, molecule, internal, engine, directory, params)
        try:
            progress = optimizer.optimizeGeometry()
        except ReleaseCandidateReached:
            progress = optimizer.progress
            progress.write(output_xyz)
            _write_optimization_status("release_candidate")
            return optimizer
        except GeomOptNotConvergedError:
            # geomeTRIC may have appended a rejected final trial. The controller
            # checkpoint is written only at accepted geometries, so use it as
            # the authoritative recovery point and refresh its physical result.
            checkpoint = json.loads(Path(_CONTROLLER_STATE_FILE).read_text())
            accepted_angstrom = np.asarray(checkpoint["positions_angstrom"], dtype=float)
            accepted_bohr = accepted_angstrom.reshape(-1) / Bohr
            engine.update_atoms(accepted_bohr)
            engine.clearCalcs()
            result = engine.calc(accepted_bohr, optimizer.dirname)
            optimizer.X = accepted_bohr.copy()
            optimizer.E = result["energy"]
            optimizer.gradx = result["gradient"]
            optimizer.G = optimizer.IC.calcGrad(optimizer.X, optimizer.gradx).flatten()
            progress = optimizer.progress
            if progress.xyzs:
                progress.xyzs[-1] = accepted_angstrom.copy()
                progress.qm_energies[-1] = optimizer.E
                progress.qm_grads[-1] = optimizer.gradx.copy()
            progress.write(output_xyz)
            _write_optimization_status("cycle_exceeded")
            return optimizer

        progress.write(output_xyz)
        decision = calculator.bias_controller.decision
        if (decision is None or decision.alpha >= calculator.alpha_max
                or np.isclose(decision.alpha, calculator.alpha_max,
                              rtol=1e-10, atol=1e-14)):
            _write_optimization_status("converged")
            return optimizer

        # Do not checkpoint a newly proposed alpha when the segment budget is
        # exhausted: it would describe a strength that was never evaluated.
        if segment + 1 >= max_segments:
            _write_optimization_status("cycle_exceeded")
            return optimizer

        # A biased local minimum is not the endpoint of the search. Re-estimate
        # the resistance at that optimized structure, add another load
        # increment, then restart with a fresh geomeTRIC Hessian.
        coordinates = optimizer.X.copy()
        engine.update_atoms(coordinates)
        calculator.accept_geometry(engine.ase_atoms, observe_release=False)

    tracker = calculator.release_tracker
    _write_optimization_status(
        "release_candidate" if tracker is not None and tracker.evidence.state == "CANDIDATE"
        else "cycle_exceeded"
    )
    return optimizer


def _write_optimization_status(status):
    state_path = Path("pyar_geometric_state.json")
    if state_path.exists():
        state = json.loads(state_path.read_text())
        state["optimization_status"] = status
        state_path.write_text(json.dumps(state, indent=2, sort_keys=True))


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    run_adaptive_optimization(sys.argv[1], json.loads(sys.argv[2]))
