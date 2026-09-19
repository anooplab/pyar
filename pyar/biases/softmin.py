"""Soft-minimum collective-coordinate bias for intermolecular reactions."""

import math

import numpy as np

import pyar.data.units
from pyar.biases.collective_coordinates import evaluate_contact_coordinate

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

    q, dq_dR, _ = evaluate_contact_coordinate(
        fragment_indices, atoms_list, coordinates, kind="softmin", beta=beta
    )
    return alpha * q, -alpha * np.asarray(dq_dR, dtype=float)
