#!/usr/bin/env python3
"""Importable ``pyar-energy-table`` entrypoint.

This utility prints a relative-energy table for one or more XYZ files. It is
used as a lightweight inspection tool for comparing geometries without
running a full workflow.
"""

import argparse
import sys
from pathlib import Path

from pyar.Molecule import Molecule
from pyar.data_analysis import clustering


def main():
    """Read XYZ files, attach energies, and print a ranked energy table."""
    parser = argparse.ArgumentParser(
        prog="pyar-energy-table",
        description="Print a relative-energy table from one or more XYZ files.",
    )
    parser.add_argument(
        "input_files",
        metavar="xyz",
        nargs="+",
        help="XYZ files to report in energy order",
    )
    args = parser.parse_args()

    molecules = []
    for xyz_file in args.input_files:
        molecule = Molecule.from_xyz(xyz_file)
        molecule.energy = clustering.read_energy_from_xyz_file(xyz_file)
        molecule.name = Path(xyz_file).stem
        molecules.append(molecule)

    clustering.print_energy_table(
        molecules,
        stream=sys.stdout,
        title="Relative energy table:",
    )


if __name__ == "__main__":
    main()
