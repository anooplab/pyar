"""Strict reader for concatenated XYZ geometry trajectories."""

from pathlib import Path
import math

import numpy as np
from ase.data import atomic_numbers


def parse_xyz_trajectory(path, expected_atoms=None, expected_symbols=None):
    """Return ``(symbols, coordinates, comment)`` frames from a multi-XYZ file.

    Empty lines and standalone ``>`` separators are accepted between frames.
    Every frame must be complete, use the same ordered elements, and contain
    finite coordinates. Atom mapping is positional; atoms are never permuted.
    """
    path = Path(path)
    lines = path.read_text(encoding="utf-8").splitlines()
    frames = []
    cursor = 0
    reference_symbols = list(expected_symbols) if expected_symbols is not None else None
    while cursor < len(lines):
        while cursor < len(lines) and lines[cursor].strip() in {"", ">"}:
            cursor += 1
        if cursor >= len(lines):
            break
        try:
            count = int(lines[cursor].strip())
        except ValueError as exc:
            raise ValueError(f"Invalid XYZ atom count in {path} at line {cursor + 1}") from exc
        if count < 1 or (expected_atoms is not None and count != expected_atoms):
            raise ValueError(f"XYZ frame has {count} atoms; expected {expected_atoms}")
        if cursor + count + 2 > len(lines):
            raise ValueError(f"Truncated XYZ frame in {path}")
        symbols = []
        coordinates = []
        for line in lines[cursor + 2:cursor + count + 2]:
            fields = line.split()
            if len(fields) < 4:
                raise ValueError(f"Malformed XYZ coordinate line in {path}")
            symbol = fields[0].capitalize()
            if symbol not in atomic_numbers or symbol == "X":
                raise ValueError(f"Invalid XYZ element {fields[0]!r} in {path}")
            try:
                coordinate = [float(value) for value in fields[1:4]]
            except ValueError as exc:
                raise ValueError(f"Malformed XYZ coordinates in {path}") from exc
            if not all(math.isfinite(value) for value in coordinate):
                raise ValueError(f"Non-finite XYZ coordinates in {path}")
            symbols.append(symbol)
            coordinates.append(coordinate)
        if reference_symbols is None:
            reference_symbols = symbols
        if symbols != reference_symbols:
            raise ValueError("XYZ trajectory atom symbols or order do not match")
        frames.append((symbols, np.asarray(coordinates, dtype=float), lines[cursor + 1]))
        cursor += count + 2
    if not frames:
        raise ValueError(f"XYZ trajectory contains no frames: {path}")
    return frames
