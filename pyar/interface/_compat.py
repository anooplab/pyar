"""Compatibility helpers for backend modules relocated to ``pyar.backends``."""

import runpy
import sys
from importlib import import_module


def expose_backend(module_name, backend_name):
    """Alias normal imports and preserve direct execution of legacy modules."""
    if module_name == "__main__":
        runpy.run_module(backend_name, run_name="__main__")
        return
    sys.modules[module_name] = import_module(backend_name)
