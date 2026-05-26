"""Compatibility alias for :mod:`pyar.workflows.reaction`."""

import sys
from importlib import import_module

sys.modules[__name__] = import_module("pyar.workflows.reaction")
