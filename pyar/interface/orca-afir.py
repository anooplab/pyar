"""Compatibility launcher for the moved legacy ORCA AFIR backend script."""

import runpy
from pathlib import Path

_backend_script = Path(__file__).resolve().parents[1] / "backends" / "orca-afir.py"

if __name__ == "__main__":
    runpy.run_path(str(_backend_script), run_name="__main__")
else:
    _exports = runpy.run_path(str(_backend_script))
    globals().update(
        {name: value for name, value in _exports.items() if not name.startswith("__")}
    )
