"""Evidence helpers for the risk-gated advanced recovery pilot.

These helpers do not select a new Hamiltonian, alter requested convergence
criteria, or authorize an electronic-root change. Optimizer retries remain
limited to routes with a validated accepted-frame contract.
"""
from dataclasses import dataclass
import math

import numpy as np


@dataclass(frozen=True)
class GradientConsistency:
    evaluations: int
    rms_absolute_error_hartree_per_bohr: float
    max_absolute_error_hartree_per_bohr: float
    relative_rms_error: float
    consistent: bool


def finite_difference_gradient(evaluate_energy, coordinates_bohr, analytic_gradient,
                               *, step_bohr=1e-4, max_evaluations=512,
                               relative_tolerance=.02):
    """Compare a calculator gradient with bounded central energy differences.

    ``evaluate_energy`` receives an N-by-3 Cartesian array in Bohr and returns
    energy in Hartree. This is an explicit diagnostic; callers must budget the
    6N single-point evaluations and must not use it to weaken an optimizer gate.
    """
    coordinates = np.asarray(coordinates_bohr, dtype=float)
    gradient = np.asarray(analytic_gradient, dtype=float)
    if coordinates.ndim != 2 or coordinates.shape[1:] != (3,) or not coordinates.size:
        raise ValueError('Coordinates must have shape (N, 3)')
    if gradient.shape != coordinates.shape:
        raise ValueError('Analytic gradient shape must match coordinates')
    if not np.isfinite(coordinates).all() or not np.isfinite(gradient).all():
        raise ValueError('Coordinates and analytic gradient must be finite')
    if (isinstance(step_bohr, bool) or not math.isfinite(step_bohr) or step_bohr <= 0
            or type(max_evaluations) is not int or max_evaluations < 1
            or isinstance(relative_tolerance, bool) or not math.isfinite(relative_tolerance)
            or relative_tolerance <= 0):
        raise ValueError('Finite-difference controls must be finite and positive')
    required = 2 * coordinates.size
    if required > max_evaluations:
        raise ValueError(f'Finite-difference check needs {required} evaluations; cap is {max_evaluations}')
    numerical = np.empty_like(gradient)
    for atom in range(len(coordinates)):
        for axis in range(3):
            plus = coordinates.copy()
            minus = coordinates.copy()
            plus[atom, axis] += step_bohr
            minus[atom, axis] -= step_bohr
            raw_plus, raw_minus = evaluate_energy(plus), evaluate_energy(minus)
            if isinstance(raw_plus, bool) or isinstance(raw_minus, bool):
                raise ValueError('Finite-difference energies must be finite Hartree scalars')
            e_plus, e_minus = float(raw_plus), float(raw_minus)
            if not math.isfinite(e_plus) or not math.isfinite(e_minus):
                raise ValueError('Finite-difference energies must be finite Hartree scalars')
            numerical[atom, axis] = (e_plus - e_minus) / (2 * step_bohr)
    error = gradient - numerical
    rms = float(np.sqrt(np.mean(error**2)))
    maximum = float(np.max(np.abs(error)))
    scale = max(float(np.sqrt(np.mean(numerical**2))), 1e-8)
    relative = rms / scale
    return GradientConsistency(required, rms, maximum, relative,
                               relative <= relative_tolerance)
