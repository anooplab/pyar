"""Compatibility alias for :mod:`pyar.backends.subprocess_utils`."""

from pyar.interface._compat import expose_backend

expose_backend(__name__, "pyar.backends.subprocess_utils")
