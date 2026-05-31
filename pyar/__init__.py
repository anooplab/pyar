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
import tomllib

from pyar.core.molecule import Molecule


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

__author__ = 'Anakuthil Anoop'
__credits__ = 'IIT Kharagpur'

__all__ = [
    "Molecule",
    "__version__",
]
