"""Selection reporting helpers."""

from __future__ import annotations

import operator
import re

import numpy as np

__all__ = [
    "plot_energy_histogram",
    "print_energy_table",
    "read_energy_from_xyz_file",
]


def print_energy_table(molecules, stream=None, title=None):
    """Report energies with relative values against the global minimum."""
    entries = []
    for molecule in molecules:
        name = molecule.name
        energy = float(molecule.energy)
        relative_path = getattr(molecule, "relative_path", None)
        if relative_path is None:
            relative_path = getattr(molecule, "result_path", None)
        entries.append((name, energy, relative_path))

    if not entries:
        return

    ref = min(energy for _, energy, _ in entries)
    lines = []
    if title:
        lines.append(title)
    show_paths = any(relative_path for _, _, relative_path in entries)
    header = f"{'Name':>35}  {'Energy':>12}  {'R. E. (kcal/mol)':>18}"
    if show_paths:
        header += "  Relative path"
    lines.append(header)
    for name, energy, relative_path in sorted(entries, key=operator.itemgetter(1)):
        row = f"{name:>35}  {energy:12.6f}  {(energy - ref) * 627.51:18.2f}"
        if show_paths:
            row += f"  {relative_path or name}"
        lines.append(row)
    min_name, _, min_path = min(entries, key=lambda entry: entry[1])
    minimum_label = f"{min_name}"
    if min_path:
        minimum_label += f" ({min_path})"
    lines.append(f"Global minimum: {minimum_label} ({ref:12.6f} Eh)")
    lines.append("")

    if stream is not None:
        for line in lines:
            print(line, file=stream)
    else:
        from pyar.selection import clustering

        for line in lines:
            clustering.cluster_logger.info(line)


def read_energy_from_xyz_file(xyz_file):
    """Return the trailing numeric energy token from an XYZ comment line."""
    try:
        with open(xyz_file, "r") as fr:
            lines = fr.readlines()
            second_line = lines[1].strip()
            energy_matches = re.findall(r"-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?", second_line)
            if energy_matches:
                energy = float(energy_matches[-1])
            else:
                raise ValueError("No numeric energy found")
    except (IndexError, ValueError):
        with open(xyz_file, "r") as fr:
            comments_line = fr.readlines()[1].rstrip()
            energy = float(re.split(r":|=|\s+", comments_line)[1])

    return energy


def plot_energy_histogram(molecules):
    """Plot a simple energy histogram for a selected pool."""
    energies = [i.energy for i in molecules]
    ref = min(energies)
    relative_energies = [(float(energy) - ref) for energy in energies]
    histogram_bin = np.linspace(0, max(relative_energies), 10)
    import matplotlib.pyplot as plt

    plt.hist(relative_energies, histogram_bin)
    plt.xlabel("Energy")
    plt.ylabel("Population")
    plt.title("Histogram of energies")
