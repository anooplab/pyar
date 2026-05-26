"""Compatibility alias for :mod:`pyar.biases.afir`."""

import sys
from importlib import import_module

sys.modules[__name__] = import_module("pyar.biases.afir")
