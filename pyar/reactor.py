"""Compatibility shim for the legacy :mod:`pyar.reactor` import path.

The reaction workflow implementation lives in :mod:`pyar.workflows.reaction`.
This module re-exports that implementation so older code that imports
``pyar.reactor`` continues to work without changes.
"""

import sys
from importlib import import_module

sys.modules[__name__] = import_module("pyar.workflows.reaction")
