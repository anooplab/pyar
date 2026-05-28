"""Compatibility wrapper for :mod:`pyar.state.solvation`.

Legacy code imported solvation restart-state classes from ``pyar.solvation_state``.
The real implementation now lives in :mod:`pyar.state.solvation`, and this
module re-exports it so the old import path remains valid.
"""

from pyar.state.solvation import *  # noqa: F401,F403
