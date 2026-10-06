"""Compatibility spelling for modern solute-centred microsolvation."""

import warnings

from pyar.scripts.modern_microsolvate import build_parser, main as _main


def main(argv=None, *, prog=None):
    warnings.warn(
        "`pyar solvate` is a compatibility alias; use `pyar microsolvate` for "
        "solute-centred first-shell sampling or `pyar grow` for generic growth.",
        FutureWarning, stacklevel=2,
    )
    return _main(argv, prog=prog or "pyar solvate")


def modern_main(argv=None, *, prog=None):
    return main(argv, prog=prog)
