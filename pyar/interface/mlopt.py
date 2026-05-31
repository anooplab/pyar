"""Compatibility shim for the legacy ``pyar.interface.mlopt`` script path.

Use :mod:`pyar.backends.mlopt` instead.
"""

from pyar.backends.mlopt import main


if __name__ == "__main__":
    main()
