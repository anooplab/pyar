"""Soft-minimum collective-coordinate bias for intermolecular reactions."""

from itertools import product
import math

import autograd.numpy as np
from autograd import grad

import pyar.data.units
from pyar.biases.afir import get_covalent_radius

__all__ = ["softmin"]


def softmin(fragment_indices, atoms_list, coordinates, gamma, beta=1.0):
    """Return soft-min bias energy and forces for two fragments.

    ``gamma`` is the AFIR-calibrated bias strength in kJ mol-1. ``beta`` is the
    inverse-length localization parameter in bohr units. Lower values spread
    the bias over more interfragment contacts; higher values increasingly focus
    it on the shortest surface gap.
    """
    try:
        beta = float(beta)
    except (TypeError, ValueError):
        raise ValueError(f"Invalid soft-min beta: {beta!r}") from None
    if not math.isfinite(beta) or beta <= 0.0:
        raise ValueError("Soft-min beta must be a finite positive number")

    epsilon = pyar.data.units.kilojoules2atomic_units(1.0061)
    r_zero = pyar.data.units.angstrom2bohr(3.8164)
    gamma = pyar.data.units.kilojoules2atomic_units(gamma)
    if gamma == 0.0:
        alpha = 0.0
    else:
        alpha = gamma / (
            (2 ** (-1.0 / 6.0) - (1 + np.sqrt(1 + gamma / epsilon)) ** (-1.0 / 6.0))
            * r_zero
        )

    fragment_one, fragment_two = [coordinates[indices, :] for indices in fragment_indices]
    atom_symbols = np.array(atoms_list)
    symbols_one, symbols_two = [atom_symbols[indices] for indices in fragment_indices]
    radii_one = [get_covalent_radius(symbol) for symbol in symbols_one]
    radii_two = [get_covalent_radius(symbol) for symbol in symbols_two]

    def restraint_energy_for_fragments(f_one, f_two):
        distances = np.array([np.linalg.norm(a - b) for a, b in product(f_one, f_two)])
        radii = np.array([a + b for a, b in product(radii_one, radii_two)])
        scaled_gaps = -beta * (distances - radii)
        max_scaled_gap = np.max(scaled_gaps)
        soft_minimum = -(
            max_scaled_gap + np.log(np.mean(np.exp(scaled_gaps - max_scaled_gap)))
        ) / beta
        return alpha * soft_minimum

    bias_energy = restraint_energy_for_fragments(fragment_one, fragment_two)
    gradient_one = grad(restraint_energy_for_fragments, argnum=0)(fragment_one, fragment_two)
    gradient_two = grad(restraint_energy_for_fragments, argnum=1)(fragment_one, fragment_two)
    bias_gradient = -np.concatenate((gradient_one, gradient_two))
    return bias_energy, bias_gradient
