"""Version-pinned accepted-step observer in an isolated geomeTRIC process.

Uses the standard driver/parameters. It adds immutable accepted observations;
no bias, force modification, Hessian switching or orbital adoption is performed.
"""
import json
from pathlib import Path
import sys

import geometric
import geometric.optimize as driver
import numpy as np
from geometric.errors import GeomOptNotConvergedError
from geometric.optimize import OPT_STATE, Optimizer
from geometric.nifty import bohr2ang

from pyar.persistence import write_json
from pyar.recovery.models import OptimizationRequest


class AcceptedOptimizer(Optimizer):
    def _record(self):
        settings = dict(self.engine.calculator.qc_params, fragments=self.engine.calculator.fragment_indices)
        request = OptimizationRequest(
            Path(sys.argv[1]).stem.removeprefix('trial_'), tuple(self.molecule.elem),
            self.molecule.xyzs[0], 'xtb', 'geometric', 'minimum', json.dumps(settings))
        gradient = self.gradx.reshape((-1, 3))
        norms = np.linalg.norm(gradient, axis=1)
        value = dict(schema_version=1, version=geometric.__version__,
            objective_fingerprint=request.objective_fingerprint, symbols=self.molecule.elem,
            coordinates_angstrom=(self.X.reshape((-1, 3)) * bohr2ang).tolist(),
            energy_hartree=float(self.E), gradient_hartree_per_bohr=gradient.tolist(),
            gradient_rms=float(np.sqrt(np.mean(norms**2))), gradient_max=float(max(norms)),
            step=int(self.Iteration))
        write_json(f'ore_accepted_{self.Iteration:06d}.json', value, overwrite=False)
        self._ore_files = getattr(self, '_ore_files', []) + [f'ore_accepted_{self.Iteration:06d}.json']
        write_json('ore_accepted.json', {'schema_version': 1, 'files': self._ore_files})

    def prepareFirstStep(self):
        super().prepareFirstStep()
        self._record()

    def evaluateStep(self):
        trial = self.X.copy()
        super().evaluateStep()
        # A limit is checked before the last trial's acceptance decision in
        # geomeTRIC 1.1.1. It is deliberately NOT checkpointed as accepted.
        if self.state != OPT_STATE.FAILED and np.array_equal(self.X, trial):
            self._record()


def main():
    if geometric.__version__ != '1.1.1':
        raise ValueError('SAFE recovery requires geomeTRIC 1.1.1')
    arguments = json.loads(sys.argv[2])
    settings = arguments['qc_params']
    from pyar.recovery.policy import validate_settings
    validate_settings(settings)
    # This interpreter owns one optimization only; the standard driver keeps
    # all coordinate setup, logging, thresholds and convergence checks.
    driver.Optimizer = AcceptedOptimizer
    try:
        driver.run_optimizer(input=sys.argv[1], engine='ase',
            ase_class='pyar.backends.geometric.PyarGeometricCalculator',
            ase_kwargs=json.dumps(arguments), coordsys='tric', maxiter=settings['opt_cycles'],
            converge=['set', {'loose': 'GAU_LOOSE', 'normal': 'GAU', 'tight': 'GAU_TIGHT'}[settings['opt_threshold']]])
    except GeomOptNotConvergedError:
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
