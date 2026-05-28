#!/usr/bin/env python3
"""Benchmark clustering/selection algorithms on one or more XYZ pools.

This utility compares the geometry-selection algorithms exposed by
``pyar.data_analysis.clustering`` and reports selection quality, diversity,
and runtime for one or more pools of XYZ files.
"""

from __future__ import annotations

import argparse
import time
from collections import Counter
from pathlib import Path
from unittest import mock

import numpy as np

from pyar.Molecule import Molecule
from pyar import representations
from pyar.data_analysis import clustering


DEFAULT_ALGORITHMS = [
    "hdbscan",
    "optics",
    "affinity",
    "agglomerative",
    "mean_shift",
    "spectral",
    "dbscan",
    "kmeans",
    "gaussian_mixture",
    "rbf_kernel",
    "maxmin",
]


def _discover_xyz_files(path: Path):
    if path.is_file():
        return [path]
    if path.is_dir():
        return sorted(path.rglob("result_*.xyz"))
    return []


def _load_pool(path: Path):
    files = _discover_xyz_files(path)
    molecules = []
    for each_file in files:
        mol = Molecule.from_xyz(str(each_file))
        mol.energy = clustering.read_energy_from_xyz_file(str(each_file))
        molecules.append(mol)
    return molecules


def _fingerprint_vector(molecule):
    vector = representations.fingerprint(molecule.atoms_list, molecule.coordinates)
    vector = np.real_if_close(vector)
    return np.asarray(vector, dtype=float).ravel()


def _pairwise_mean_distance(molecules):
    if len(molecules) < 2:
        return 0.0
    fingerprints = [_fingerprint_vector(mol) for mol in molecules]
    distances = []
    for index, left in enumerate(fingerprints):
        for right in fingerprints[index + 1:]:
            distances.append(float(np.linalg.norm(left - right)))
    return float(np.mean(distances)) if distances else 0.0


def _coverage_distance(pool, selected):
    if not pool or not selected:
        return float("inf")
    selected_fps = [_fingerprint_vector(mol) for mol in selected]
    distances = []
    for molecule in pool:
        fp = _fingerprint_vector(molecule)
        distances.append(min(float(np.linalg.norm(fp - sel_fp)) for sel_fp in selected_fps))
    return float(np.mean(distances))


def _tetrahedral_score(molecule):
    counts = Counter(molecule.atoms_list)
    if counts.get("C", 0) != 1 or counts.get("H", 0) != 4 or len(counts) != 2:
        return None

    carbon_index = molecule.atoms_list.index("C")
    carbon = molecule.coordinates[carbon_index]
    hydrogens = [coord for atom, coord in zip(molecule.atoms_list, molecule.coordinates) if atom == "H"]
    angles = []
    for i in range(len(hydrogens)):
        for j in range(i + 1, len(hydrogens)):
            v1 = hydrogens[i] - carbon
            v2 = hydrogens[j] - carbon
            cosine = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
            cosine = float(np.clip(cosine, -1.0, 1.0))
            angles.append(np.degrees(np.arccos(cosine)))
    target = 109.4712
    return float(sum(abs(angle - target) for angle in angles))


def _benchmark_algorithm(pool, algorithm, max_seeds):
    clustering._MBTR_RUNTIME_DISABLED = False
    clustering._MBTR_DISABLE_REASON = None

    with mock.patch("pyar.data_analysis.clustering._load_basin_registry", return_value=[]):
        start = time.perf_counter()
        selected = clustering.choose_geometries(
            pool,
            maximum_number_of_seeds=max_seeds,
            persist_basin_memory=False,
            algorithm=algorithm,
        )
        runtime = time.perf_counter() - start

    best_energy = min(float(molecule.energy) for molecule in selected) if selected else float("inf")
    spread = _pairwise_mean_distance(selected)
    coverage = _coverage_distance(pool, selected)
    tetra_scores = [score for score in (_tetrahedral_score(m) for m in selected) if score is not None]
    best_tetra = min(tetra_scores) if tetra_scores else None
    return {
        "algorithm": algorithm,
        "selected": len(selected),
        "best_energy": best_energy,
        "spread": spread,
        "coverage": coverage,
        "best_tetra": best_tetra,
        "runtime": runtime,
        "selected_names": [m.name for m in selected],
    }


def main():
    parser = argparse.ArgumentParser(description="Benchmark clustering algorithms on XYZ pools.")
    parser.add_argument("paths", nargs="+", help="XYZ files or directories containing result_*.xyz files")
    parser.add_argument(
        "-a",
        "--algorithms",
        nargs="*",
        default=DEFAULT_ALGORITHMS,
        help="Algorithms to benchmark",
    )
    parser.add_argument(
        "-N",
        "--max-seeds",
        type=int,
        default=8,
        help="Maximum number of selected seeds",
    )
    parser.add_argument(
        "--show-selected",
        action="store_true",
        help="Print selected molecule names for each benchmark row",
    )
    args = parser.parse_args()

    pools = []
    for raw_path in args.paths:
        path = Path(raw_path)
        molecules = _load_pool(path)
        if len(molecules) < 2:
            print(f"{path}\tno data")
            continue
        pools.append((path, molecules))

    if not pools:
        raise SystemExit("No XYZ data found.")

    header = [
        "pool",
        "algorithm",
        "selected",
        "best_energy",
        "spread",
        "coverage",
        "best_tetra",
        "runtime_s",
    ]
    if args.show_selected:
        header.append("selected_names")
    print("\t".join(header))

    for path, pool in pools:
        for algorithm in args.algorithms:
            try:
                row = _benchmark_algorithm(pool, algorithm, args.max_seeds)
                values = [
                    str(path),
                    row["algorithm"],
                    str(row["selected"]),
                    f"{row['best_energy']:.6f}",
                    f"{row['spread']:.6f}",
                    f"{row['coverage']:.6f}",
                    "-" if row["best_tetra"] is None else f"{row['best_tetra']:.6f}",
                    f"{row['runtime']:.3f}",
                ]
                if args.show_selected:
                    values.append(",".join(row["selected_names"]))
                print("\t".join(values))
            except Exception as exc:  # pragma: no cover - benchmark safety net
                values = [
                    str(path),
                    algorithm,
                    "error",
                    "-",
                    "-",
                    "-",
                    "-",
                    "-",
                ]
                if args.show_selected:
                    values.append(str(exc))
                print("\t".join(values))


if __name__ == "__main__":
    main()
