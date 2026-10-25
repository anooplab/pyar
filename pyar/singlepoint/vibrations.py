"""Numerical Cartesian Hessians and projected harmonic vibrational analysis.

Coordinates are never optimized or linearized. Hessians use Hartree/bohr**2;
atomic masses are in unified atomic mass units. Rigid modes are removed by an
SVD of mass-weighted translations and rotations at the evaluated geometry.
"""

import numpy as np
from scipy import constants


def gradient_metrics(gradient, maximum=4.5e-4, rms=3.0e-4):
    norms = np.linalg.norm(gradient, axis=1)
    largest = float(np.max(norms))
    average = float(np.sqrt(np.mean(norms**2)))
    return {
        "maximum_gradient_hartree_per_bohr": largest,
        "rms_gradient_hartree_per_bohr": average,
        "gradient_norm_definition": "per-atom Cartesian vector norms",
        "maximum_gradient_threshold": maximum,
        "rms_gradient_threshold": rms,
        "stationary": largest < maximum and average < rms,
    }


def numerical_hessian(evaluate, coordinates_bohr, step):
    flat = np.asarray(coordinates_bohr, dtype=float).reshape(-1)
    hessian = np.empty((flat.size, flat.size))
    for index in range(flat.size):
        plus, minus = flat.copy(), flat.copy()
        plus[index] += step
        minus[index] -= step
        gp = evaluate(plus.reshape(-1, 3), f"displacement_{index:04d}_plus")[1]
        gm = evaluate(minus.reshape(-1, 3), f"displacement_{index:04d}_minus")[1]
        hessian[:, index] = (gp.reshape(-1) - gm.reshape(-1)) / (2 * step)
    if not np.isfinite(hessian).all():
        raise ValueError("Non-finite numerical Hessian")
    asymmetry = float(np.max(np.abs(hessian - hessian.T)))
    return (hessian + hessian.T) / 2, asymmetry


def vibrational_basis(coordinates, masses):
    centered = coordinates - np.average(coordinates, axis=0, weights=masses)
    weights = np.sqrt(masses)[:, None]
    translations = [np.tile(axis, (len(masses), 1)) * weights for axis in np.eye(3)]
    rotations = [np.cross(centered, axis) * weights for axis in np.eye(3)]
    rigid = np.column_stack([mode.reshape(-1) for mode in translations + rotations])
    u, singular, _ = np.linalg.svd(rigid, full_matrices=True)
    rank = int(np.sum(singular > singular[0] * 1e-8))
    return u[:, rank:], rank


def analyze_hessian(hessian, coordinates, masses, threshold):
    masses = np.asarray(masses, dtype=float)
    if np.any(masses <= 0) or not np.isfinite(masses).all():
        raise ValueError("Vibrational analysis requires positive finite atomic masses")
    basis, rank = vibrational_basis(np.asarray(coordinates), masses)
    weights = np.repeat(np.sqrt(masses), 3)
    weighted = hessian / np.outer(weights, weights)
    values, vectors = np.linalg.eigh(basis.T @ weighted @ basis)
    factor = np.sqrt(
        constants.physical_constants["Hartree energy"][0]
        / (constants.physical_constants["Bohr radius"][0] ** 2 * constants.atomic_mass)
    )
    factor /= 2 * np.pi * constants.c * 100
    frequencies = np.sign(values) * np.sqrt(np.abs(values)) * factor
    modes = (basis @ vectors) / weights[:, None]
    return {
        "frequencies_cm-1": frequencies.tolist(),
        "imaginary_frequencies_cm-1": frequencies[frequencies < -threshold].tolist(),
        "imaginary_frequency_count": int(np.sum(frequencies < -threshold)),
        "small_negative_frequencies_cm-1": frequencies[
            (frequencies < 0) & (frequencies >= -threshold)
        ].tolist(),
        "imaginary_frequency_threshold_cm-1": threshold,
        "rigid_mode_count": rank,
        "vibrational_mode_count": len(values),
        "masses_amu": masses.tolist(),
        "normal_modes_cartesian": modes.T.reshape(-1, len(masses), 3).tolist(),
    }


def qualification(stationary, imaginary):
    if not stationary:
        return "nonstationary"
    return {0: "minimum_confirmed", 1: "first_order_saddle_candidate"}.get(
        imaginary, "higher_order_saddle_candidate"
    )
