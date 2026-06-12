# -*- coding: utf-8 -*-
"""
    pyar
    ~~~~

    A Python code for Aggregation and Reaction.
    This includes modules for automating the tasks
    for the following:

    * exploring the unknown chemical reactions
      between two molecules
    * predicting the geometries of molecular
      aggregates and atomic clusters.
    * in addition, there are some analysis modules

    :copyright: © 2010 by AnoopLab.
    :license: GPL v3.

"""

__docformat__ = 'restructuredtext'

from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
import sys
from types import ModuleType

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - exercised on Python 3.10
    import tomli as tomllib

from pyar.core.molecule import Molecule


class _PyarPackage(ModuleType):
    """Package module that keeps ``pyar.Molecule`` bound to the class."""

    def __getattribute__(self, name):
        if name == "Molecule":
            return ModuleType.__getattribute__(self, "_canonical_molecule")
        return ModuleType.__getattribute__(self, name)


def _local_version_fallback():
    pyproject = Path(__file__).resolve().parent.parent / "pyproject.toml"
    try:
        with pyproject.open("rb") as handle:
            return tomllib.load(handle)["project"]["version"]
    except (OSError, KeyError, tomllib.TOMLDecodeError):
        return "0+unknown"


try:
    __version__ = version("pyar-chem")
except PackageNotFoundError:
    __version__ = _local_version_fallback()

sys.modules[__name__]._canonical_molecule = Molecule
sys.modules[__name__].__class__ = _PyarPackage

__author__ = 'Anakuthil Anoop'
__credits__ = 'IIT Kharagpur'

__all__ = [
    "Molecule",
    "__version__",
]
