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


try:
    __version__ = version("pyar-chem")
except PackageNotFoundError:
    __version__ = "1.1.0"
__author__ = 'Anakuthil Anoop'
__credits__ = 'IIT Kharagpur'
