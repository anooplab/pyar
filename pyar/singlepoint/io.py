"""Small, explicit artifact helpers for single-point calculations."""

import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.constants import physical_constants

BOHR_ANGSTROM = physical_constants["Bohr radius"][0] * 1e10


def write_json(path, value):
    Path(path).write_text(
        json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )


def write_xyz(path, symbols, coordinates):
    rows = [str(len(symbols)), "PyAR evaluated geometry; coordinates in Angstrom"]
    rows.extend(
        f"{s} {x:.15f} {y:.15f} {z:.15f}" for s, (x, y, z) in zip(symbols, coordinates)
    )
    Path(path).write_text("\n".join(rows) + "\n")


def geometry_hash(symbols, coordinates):
    data = {
        "symbols": list(symbols),
        "coordinates_angstrom": np.asarray(coordinates).tolist(),
    }
    return hashlib.sha256(json.dumps(data, sort_keys=True).encode()).hexdigest()


def validate_result(energy, gradient, natoms, required):
    energy = float(energy)
    if not np.isfinite(energy):
        raise ValueError("Backend returned a non-finite energy")
    if gradient is None:
        if required:
            raise ValueError("Backend did not return the requested gradient")
        return energy, None
    gradient = np.asarray(gradient, dtype=float)
    if gradient.shape != (natoms, 3) or not np.isfinite(gradient).all():
        raise ValueError(f"Backend gradient must be finite with shape ({natoms}, 3)")
    return energy, gradient
