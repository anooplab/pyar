"""Compatibility entrypoints for :mod:`pyar.backends`."""

import importlib

from pyar.backends import SF, require_executable, run_command, run_output, which, write_xyz

__all__ = [
    "ANI",
    "ANICalculationFailed",
    "ANIInterface",
    "SF",
    "require_executable",
    "which",
    "write_xyz",
    "run_command",
    "run_output",
]


def __getattr__(name):
    """Lazily preserve optional ANI exports from the legacy package."""
    if name in {"ANI", "ANICalculationFailed", "ANIInterface"}:
        ani = importlib.import_module("pyar.backends.ani")
        globals().update(
            ANI=ani.ANI,
            ANICalculationFailed=ani.ANICalculationFailed,
            ANIInterface=ani.ANIInterface,
        )
        return globals()[name]
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
