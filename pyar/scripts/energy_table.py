#!/usr/bin/env python3
"""Importable ``pyar-energy-table`` entrypoint.

This utility prints a relative-energy table for one or more XYZ files. It is
used as a lightweight inspection tool for comparing geometries without
running a full workflow.
"""

import argparse
import sys
import json
import math
from pathlib import Path

import numpy as np

from pyar.core.molecule import Molecule, parse_xyz
from pyar.selection import reports as selection_reports


def build_parser(prog=None, *, modern=False):
    """Build the shared positional-file parser."""
    parser = argparse.ArgumentParser(
        prog=prog or "pyar-energy-table",
        description="Print a relative-energy table from one or more XYZ files.",
    )
    parser.add_argument(
        "input_files",
        metavar="xyz",
        nargs="+",
        help="XYZ files to report in energy order",
    )
    if modern:
        parser.add_argument("--json", action="store_true", help="Print machine-readable JSON")
    return parser


def main(argv=None, *, prog=None, modern=False):
    """Legacy and modern interfaces share loading and human reporting."""
    parser = build_parser(prog, modern=modern)
    args = parser.parse_args(argv)

    molecules = []
    for xyz_file in args.input_files:
        try:
            _, coordinates, _, _, _ = parse_xyz(xyz_file)
            if not np.isfinite(coordinates).all():
                raise ValueError(f"Non-finite coordinates in {xyz_file}")
            molecule = Molecule.from_xyz(xyz_file)
        except (ValueError, OSError, KeyError) as exc:
            parser.error(f"Invalid XYZ input {xyz_file}: {exc}")
        try:
            molecule.energy = selection_reports.read_energy_from_xyz_file(xyz_file)
            if not math.isfinite(molecule.energy):
                raise ValueError("non-finite energy")
        except (ValueError, IndexError) as exc:
            parser.error(f"Could not read an energy from: {xyz_file}\n"
                         "Expected a numeric energy in the XYZ comment line.")
        molecule.name = Path(xyz_file).stem
        molecule.relative_path = str(Path(xyz_file))
        molecules.append(molecule)

    if getattr(args, "json", False):
        ranked = sorted(molecules, key=lambda molecule: molecule.energy)
        minimum = ranked[0].energy
        print(json.dumps({"structures": [
            {"file": molecule.relative_path, "energy_hartree": molecule.energy,
             "relative_energy_kcal_mol": (molecule.energy - minimum)
             * selection_reports.HARTREE_TO_KCAL_MOL}
            for molecule in ranked], "minimum": ranked[0].relative_path}, allow_nan=False))
        return

    selection_reports.print_energy_table(
        molecules,
        stream=sys.stdout,
        title="Relative energy table:",
    )


def modern_main(argv=None, *, prog=None):
    """Modern rendering options with the same energy loading/reporting path."""
    return main(argv, prog=prog or "pyar energies", modern=True)


if __name__ == "__main__":
    main()
