"""Minimal local compatibility stub for the external MLatom package."""

from __future__ import annotations

from . import constants, data, interface_MLatomF, models, plot, shell_cmd, simulations, stats, stopper, xyz


def optimize_geometry(*args, **kwargs):
    """Placeholder geometry optimization entry point."""
    return None


def irc(*args, **kwargs):
    """Placeholder IRC entry point."""
    return None


def freq(*args, **kwargs):
    """Placeholder frequency entry point."""
    return None


def thermochemistry(*args, **kwargs):
    """Placeholder thermochemistry entry point."""
    return None


def generate_initial_conditions(*args, **kwargs):
    """Placeholder initial-condition generator."""
    return None


__all__ = [
    "constants",
    "data",
    "freq",
    "generate_initial_conditions",
    "interface_MLatomF",
    "irc",
    "models",
    "optimize_geometry",
    "plot",
    "shell_cmd",
    "simulations",
    "stats",
    "stopper",
    "thermochemistry",
    "xyz",
]
