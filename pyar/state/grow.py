"""Atomic, request-validated restart state for fixed-seed growth."""
import json
import os
import tempfile
from pathlib import Path
from pyar.state.aggregate import _json_value


class GrowStateError(RuntimeError):
    """An existing growth run cannot safely be reused."""


def write_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", dir=path.parent, prefix=".state-",
                                         suffix=".json", delete=False) as fp:
            temporary = fp.name
            json.dump(_json_value(data), fp, indent=2, sort_keys=True, allow_nan=False)
            fp.write("\n")
            fp.flush()
            os.fsync(fp.fileno())
        os.replace(temporary, path)
    finally:
        if temporary and os.path.exists(temporary):
            os.unlink(temporary)
