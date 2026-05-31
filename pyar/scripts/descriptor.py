#!/usr/bin/env python3
"""Importable ``pyar-descriptor`` entrypoint.

This utility computes a compact cluster-shape descriptor for each XYZ file
and writes per-structure descriptor files alongside summary CSV output.
"""

import argparse
import glob
import os
import sys
import warnings

import numpy as np
from scipy.spatial import ConvexHull

from pyar.optional_dependencies import optional_dependency_error


def _require_descriptor_dependencies():
    try:
        import MDAnalysis as mda  # noqa: F401
        import pandas as pd
        from ase.io import read, write
    except ImportError as exc:
        module_name = getattr(exc, "name", None) or "MDAnalysis"
        raise optional_dependency_error(
            module_name,
            feature="descriptor script",
            extra="selection",
        ) from exc
    warnings.filterwarnings("ignore", category=UserWarning, module="MDAnalysis.analysis.base")
    return pd, read, write


def calculate_properties(atoms):
    """Return the basic shape descriptors used to characterize one cluster."""
    cluster_size = len(atoms)
    points = atoms.get_positions()
    hull = ConvexHull(points)
    volume = hull.volume
    surface_area = hull.area
    distances = np.linalg.norm(points[:, np.newaxis, :] - points[np.newaxis, :, :], axis=-1)
    max_length = np.max(distances)
    rgyr = np.sqrt(np.mean(np.sum((points - np.mean(points, axis=0))**2, axis=1)))
    return cluster_size, volume, surface_area, max_length, rgyr


def create_combined_descriptor(properties):
    """Combine the basic descriptors into a single normalized score."""
    normalized = np.array(properties) / np.sum(properties)
    return np.prod(normalized)


def main(args=None):
    """Compute and write descriptors for one or more XYZ files."""
    if args is None:
        parser = argparse.ArgumentParser(description="Analyze molecular cluster XYZ files.")
        parser.add_argument("input_files", metavar='files', type=str, nargs='+',
                            help='input coordinate files (supports wildcards)')
        args = parser.parse_args()
    xyz_files = []
    for pattern in args.input_files:
        xyz_files.extend(glob.glob(pattern))
    if not xyz_files:
        print("No XYZ files found.")
        sys.exit(1)
    pd, read, write = _require_descriptor_dependencies()
    data = []
    unique_descriptors = {}
    unique_atoms = []
    duplicate_atoms = []
    for filename in xyz_files:
        atoms = read(filename)
        properties = calculate_properties(atoms)
        combined_descriptor = create_combined_descriptor(properties)
        basename = os.path.splitext(os.path.basename(filename))[0]
        with open(f"{basename}.mb", 'w') as f:
            f.write(f"{combined_descriptor}\n")
        if combined_descriptor not in unique_descriptors:
            unique_descriptors[combined_descriptor] = filename
            data.append([filename] + list(properties) + [combined_descriptor])
            unique_atoms.append(atoms)
        else:
            duplicate_atoms.append(atoms)
    columns = ["Filename", "Cluster Size", "Volume (Å³)", "Surface Area (Å²)", "Maximum Length (Å)", "Radius of Gyration (Å)", "Combined Descriptor"]
    pd.DataFrame(data, columns=columns).to_csv("unique_cluster_properties.csv", index=False)
    write("unique_files.xyz", unique_atoms)
    write("duplicate_files.xyz", duplicate_atoms)
    print(f"Processed {len(xyz_files)} XYZ files.")


if __name__ == "__main__":
    main()
