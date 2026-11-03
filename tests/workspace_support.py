"""Test-only adapters for legacy relative-file doubles at explicit workspaces.

Production never changes process cwd. Mocked in-process scientific engines need
this fixture to retain their patches; separate tests exercise real subprocess
isolation, including native xTB. These helpers are not production entry points.
"""

from contextlib import contextmanager
from functools import wraps
import os
from pathlib import Path


@contextmanager
def working_directory(directory):
    previous = Path.cwd()
    Path(directory).mkdir(parents=True, exist_ok=True)
    os.chdir(directory)
    try:
        yield
    finally:
        os.chdir(previous)


def in_test_workspace(function):
    @wraps(function)
    def execute(*args, directory=None, cwd=None, **kwargs):
        directory = directory if directory is not None else cwd
        if directory is None:
            return function(*args, **kwargs)
        with working_directory(directory):
            return function(*args, **kwargs)

    return execute


def run_legacy_scan(*args, **kwargs):
    """Exercise the unchanged layout-less scan contract through its new context."""
    from pyar.workflows import scan_bond
    bound = scan_bond._SCAN_SIGNATURE.bind(*args, **kwargs)
    return scan_bond._run_scan_bond(scan_bond._scan_options(bound.arguments))
