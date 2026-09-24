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
        changed = self.engine.calculator.accept_geometry(self.engine.ase_atoms)
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

    def optimizeGeometry(self):
        progress = super().optimizeGeometry()
        controller = self.engine.calculator.bias_controller
        decision = controller.decision
        if (decision is not None and controller.policy == "adaptive"
                and decision.alpha < controller.alpha_max
                and not np.isclose(decision.alpha, controller.alpha_max,
                                   rtol=1e-10, atol=1e-14)):
            raise RuntimeError(
                "Adaptive bias stalled below its force ceiling: the driving "
                "force is too small for the optimizer convergence tolerances. "
                "Increase --bias-alpha-margin or tighten --opt-threshold. "
                "Increasing --bias-max alone does not increase the applied force."
            )
        return progress


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
    internal = DelocalizedInternalCoordinates(
        molecule, build=True, connect=False, addcart=False
    )
    threshold = {"loose": "GAU_LOOSE", "normal": "GAU", "tight": "GAU_TIGHT"}
    output_xyz = f"{Path(input_xyz).stem}_optim.xyz"
    params = OptParams(
        maxiter=qc_params.get("opt_cycles") if qc_params.get("opt_cycles") is not None else 300,
        convergence_set=threshold.get(qc_params.get("opt_threshold"), "GAU"),
        xyzout=output_xyz,
    )
    directory = f"{Path(input_xyz).stem}.tmp"
    Path(directory).mkdir(exist_ok=True)
    optimizer = AdaptiveOptimizer(coordinates, molecule, internal, engine, directory, params)
    try:
        progress = optimizer.optimizeGeometry()
    except GeomOptNotConvergedError:
        # geomeTRIC has already evaluated and persisted the last accepted
        # geometry. Preserve that endpoint so the reaction workflow can run
        # its unbiased relaxation instead of discarding a useful candidate.
        progress = optimizer.progress
        progress.write(output_xyz)
        state_path = Path("pyar_geometric_state.json")
        if state_path.exists():
            state = json.loads(state_path.read_text())
            state["optimization_status"] = "cycle_exceeded"
            state_path.write_text(json.dumps(state, indent=2, sort_keys=True))
        return optimizer
    progress.write(output_xyz)
    return optimizer


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    run_adaptive_optimization(sys.argv[1], json.loads(sys.argv[2]))
