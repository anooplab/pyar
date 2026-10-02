"""Compare PyAR clustering choices on published LJ13 and Au13 structures.

The source records are minima/entries, not independently verified basin labels.
This harness therefore reports selection counts, energies, and geometric
descriptor coverage without claiming basin recall.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import tarfile
import zipfile
import warnings
from pathlib import Path

import numpy as np

from pyar.core.molecule import Molecule
from pyar.scripts.scientific_clustering_benchmark import (
    ALGORITHMS,
    FEATURES,
    _cluster_first_select,
)


JOULES_PER_EV = 1.602176634e-19


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_lj13_archive(path: Path) -> tuple[list[Molecule], list[dict]]:
    """Read CCD LJ13 minima; coordinates and energies use reduced LJ units."""
    with tarfile.open(path, "r:bz2") as archive:
        try:
            energy_rows = archive.extractfile("min.data").read().decode().splitlines()
            coordinate_rows = archive.extractfile("extractedmin").read().decode().splitlines()
        except (AttributeError, KeyError) as exc:
            raise ValueError("CCD archive must contain min.data and extractedmin") from exc
    energies = [float(row.split()[0]) for row in energy_rows if row.strip()]
    coordinates = np.asarray(
        [[float(value) for value in row.split()] for row in coordinate_rows if row.strip()],
        dtype=float,
    )
    atom_count = 13
    if coordinates.shape != (len(energies) * atom_count, 3):
        raise ValueError(
            f"CCD coordinate shape {coordinates.shape} does not match "
            f"{len(energies)} structures x {atom_count} atoms"
        )
    molecules, metadata = [], []
    for index, energy in enumerate(energies):
        xyz = coordinates[index * atom_count:(index + 1) * atom_count]
        if not np.isfinite(xyz).all() or not np.isfinite(energy):
            raise ValueError(f"Non-finite CCD values in minimum {index}")
        name = f"lj13_min_{index:04d}"
        molecules.append(Molecule(["Ar"] * atom_count, xyz, name=name, energy=energy))
        metadata.append({"source_id": index, "energy": energy, "energy_unit": "epsilon",
                         "coordinate_unit": "sigma", "label": name})
    return molecules, metadata


def _unwrap_nomad_record(record):
    archive = record.get("archive", record)
    metadata = archive.get("metadata", {})
    try:
        atoms = archive["run"][0]["system"][-1]["atoms"]
        energy_joule = archive["run"][0]["calculation"][-1]["energy"]["total"]["value"]
    except (KeyError, IndexError, TypeError) as exc:
        raise ValueError(f"Incomplete NOMAD archive entry {record.get('entry_id')}") from exc
    positions = np.asarray(atoms["positions"], dtype=float) * 1.0e10
    symbols = list(atoms["labels"])
    energy_ev = float(energy_joule) / JOULES_PER_EV
    geometry_optimization = (
        archive.get("results", {}).get("properties", {}).get("geometry_optimization", {})
    )
    force_n = geometry_optimization.get("final_force_maximum")
    energy_difference_j = geometry_optimization.get("final_energy_difference")
    displacement_m = geometry_optimization.get("final_displacement_maximum")
    if force_n is not None:
        metadata["final_force_maximum_ev_per_angstrom"] = float(force_n) / 1.602176634e-9
    if energy_difference_j is not None:
        metadata["final_energy_difference_ev"] = float(energy_difference_j) / JOULES_PER_EV
    if displacement_m is not None:
        metadata["final_displacement_maximum_angstrom"] = float(displacement_m) * 1.0e10
    if positions.shape != (len(symbols), 3) or not np.isfinite(positions).all():
        raise ValueError(f"Invalid coordinates in NOMAD entry {metadata.get('entry_id')}")
    if len(symbols) != 13 or set(symbols) != {"Au"}:
        raise ValueError(f"Expected Au13, found composition {symbols}")
    if not np.isfinite(energy_ev):
        raise ValueError(f"Invalid total energy in NOMAD entry {metadata.get('entry_id')}")
    return metadata, symbols, positions, energy_ev


def load_qcd_au13(path: Path) -> tuple[list[Molecule], list[dict]]:
    """Read the 30 QCD Au13 NOMAD records and convert SI units."""
    molecules, metadata_rows = [], []
    with zipfile.ZipFile(path) as archive:
        json_paths = sorted(name for name in archive.namelist()
                            if name.endswith(".json") and name.rsplit("/", 1)[-1] != "manifest.json")
        for name in json_paths:
            record = json.loads(archive.read(name))
            metadata, symbols, positions, energy = _unwrap_nomad_record(record)
            entry_id = metadata.get("entry_id") or record.get("entry_id")
            label = str(entry_id)
            molecules.append(Molecule(symbols, positions, name=label, energy=energy))
            metadata_rows.append({
                "source_id": label,
                "mainfile": metadata.get("mainfile"),
                "upload_id": metadata.get("upload_id"),
                "license": metadata.get("license"),
                "final_force_maximum_ev_per_angstrom": metadata.get("final_force_maximum_ev_per_angstrom"),
                "final_energy_difference_ev": metadata.get("final_energy_difference_ev"),
                "final_displacement_maximum_angstrom": metadata.get("final_displacement_maximum_angstrom"),
                "energy": energy,
                "energy_unit": "eV",
                "coordinate_unit": "angstrom",
                "label": label,
            })
    if not molecules:
        raise ValueError("NOMAD archive contains no JSON archive records")
    if len({row["source_id"] for row in metadata_rows}) != len(metadata_rows):
        raise ValueError("NOMAD archive has duplicate entry IDs")
    return molecules, metadata_rows


def pair_distance_fingerprint(molecule: Molecule) -> np.ndarray:
    """Sorted pair distances; permutation/rotation/translation invariant proxy."""
    xyz = np.asarray(molecule.coordinates, dtype=float)
    delta = xyz[:, None, :] - xyz[None, :, :]
    distances = np.sqrt(np.sum(delta * delta, axis=-1))
    return np.sort(distances[np.triu_indices(len(xyz), k=1)])


def run_benchmark(name: str, molecules: list[Molecule], metadata: list[dict],
                  output: Path, algorithms=ALGORITHMS, features=FEATURES,
                  maximum_seeds: int = 12):
    if len(molecules) != len(metadata) or not molecules:
        raise ValueError("Molecule and metadata rows must be nonempty and aligned")
    output.mkdir(parents=True, exist_ok=True)
    records = []
    energy_order = np.argsort([float(molecule.energy) for molecule in molecules])
    global_minimum_index = int(energy_order[0])
    top_five_indices = set(int(index) for index in energy_order[:min(5, len(energy_order))])
    fingerprints = np.stack([pair_distance_fingerprint(molecule) for molecule in molecules])
    fingerprint_scale = max(float(np.sqrt(np.mean(fingerprints ** 2))), 1.0e-12)
    for algorithm in algorithms:
        for feature in features:
            with warnings.catch_warnings(record=True) as caught_warnings:
                warnings.simplefilter("always")
                selected, diagnostics = _cluster_first_select(
                    molecules, algorithm, feature, maximum_seeds,
                    system_type="atomic-cluster",
                )
            selected_indices = [next(i for i, item in enumerate(molecules)
                                     if item.name == chosen.name) for chosen in selected]
            selected_fingerprints = fingerprints[selected_indices]
            nearest = np.sqrt(np.mean(
                (fingerprints[:, None, :] - selected_fingerprints[None, :, :]) ** 2,
                axis=2,
            )).min(axis=1)
            row = {
                "algorithm_requested": algorithm,
                "feature_requested": feature,
                "algorithm_used": diagnostics.get("algorithm_used"),
                "feature_used": diagnostics.get("feature_used"),
                "input_count": len(molecules),
                "cluster_minimum_count": diagnostics["cluster_minimum_count"],
                "selected_count": len(selected_indices),
                "selected_source_ids": [metadata[index]["source_id"] for index in selected_indices],
                "selected_energy_min": min(float(molecules[i].energy) for i in selected_indices),
                "selected_energy_max": max(float(molecules[i].energy) for i in selected_indices),
                "energy_unit": metadata[0]["energy_unit"],
                "lowest_source_energy_entry_retained": global_minimum_index in selected_indices,
                "best_selected_source_energy_rank": min(
                    int(np.where(energy_order == index)[0][0]) + 1 for index in selected_indices
                ),
                "top_five_source_energy_entries_retained": len(top_five_indices.intersection(selected_indices)),
                "runtime_warnings": [
                    {"category": warning.category.__name__, "message": str(warning.message)}
                    for warning in caught_warnings
                ],
                "mean_nearest_pair_spectrum_rms": float(np.mean(nearest)),
                "max_nearest_pair_spectrum_rms": float(np.max(nearest)),
                "normalized_mean_pair_spectrum_rms": float(np.mean(nearest) / fingerprint_scale),
                "distance_proxy": "RMS of sorted all-pairs distances; invariant to rigid motion and atom permutation; not structural RMSD",
                "diagnostics": diagnostics,
            }
            records.append(row)
            print(f"{name} {algorithm}/{feature}: selected {len(selected_indices)}; "
                  f"pair-spectrum mean {row['normalized_mean_pair_spectrum_rms']:.4f}", flush=True)
    report = {
        "schema_version": 1,
        "benchmark": name,
        "input_count": len(molecules),
        "seed_budget": maximum_seeds,
        "label_status": "source structures are preserved as distinct entries; no ground-truth basin labels are asserted",
        "source_energy_note": "energy ranks are within each source dataset as supplied; no cross-dataset energy comparison is made, and QCD entries are not filtered to a single documented electronic-structure protocol",
        "selection_protocol": "cluster, retain lowest-source-energy candidate per cluster; retain noise points; max-min trim cluster minima to seed budget",
        "coverage_metric": "nearest selected sorted pair-distance fingerprint RMS; scale-normalized by input fingerprint RMS scale",
        "coverage_limitations": "pair-distance spectra are rotation/translation/permutation invariant but can collide for non-congruent structures; this is a descriptive proxy, not RMSD or verified basin recall",
        "source_quality": {
            "entry_count_with_final_force_diagnostic": sum(
                row.get("final_force_maximum_ev_per_angstrom") is not None for row in metadata
            ),
            "maximum_reported_final_force_ev_per_angstrom": max(
                (row["final_force_maximum_ev_per_angstrom"] for row in metadata
                 if row.get("final_force_maximum_ev_per_angstrom") is not None),
                default=None,
            ),
            "interpretation": "source-reported optimization diagnostics are retained per structure; no common force cutoff is imposed because records may come from different calculation settings",
        },
        "conditions": records,
    }
    (output / "comparison.json").write_text(json.dumps(report, indent=2) + "\n")
    with (output / "structures.csv").open("w", newline="") as stream:
        columns = list(metadata[0])
        writer = csv.DictWriter(stream, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        writer.writerows(metadata)
    with (output / "comparison.csv").open("w", newline="") as stream:
        columns = [key for key, value in records[0].items() if not isinstance(value, (list, dict))]
        writer = csv.DictWriter(stream, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        writer.writerows({key: value for key, value in row.items() if key in columns} for row in records)
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lj-archive", type=Path, required=True)
    parser.add_argument("--au13-archive", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--max-seeds", type=int, default=12)
    parser.add_argument("--algorithms", nargs="+", choices=ALGORITHMS, default=list(ALGORITHMS))
    parser.add_argument("--features", nargs="+", choices=FEATURES, default=list(FEATURES))
    args = parser.parse_args(argv)
    if args.max_seeds < 1:
        parser.error("--max-seeds must be positive")
    datasets = (
        ("lj13", args.lj_archive, load_lj13_archive, "atomic-cluster"),
        ("qcd_au13", args.au13_archive, load_qcd_au13, "atomic-cluster"),
    )
    manifests = []
    for name, source, loader, system_type in datasets:
        molecules, metadata = loader(source)
        run_dir = args.output / name
        report = run_benchmark(name, molecules, metadata, run_dir,
                               args.algorithms, args.features, args.max_seeds)
        # Save indexed structures in explicit XYZ for independent inspection.
        xyz = run_dir / "structures.xyz"
        with xyz.open("w") as stream:
            for molecule, row in zip(molecules, metadata):
                stream.write(f"{len(molecule.atoms_list)}\n")
                stream.write(f"{molecule.name} energy={molecule.energy} {row['energy_unit']}\n")
                for symbol, point in zip(molecule.atoms_list, molecule.coordinates):
                    stream.write(f"{symbol} {point[0]:.12f} {point[1]:.12f} {point[2]:.12f}\n")
        manifests.append({
            "dataset": name,
            "system_type": system_type,
            "source_file": source.name,
            "source_sha256": _sha256(source),
            "structure_count": len(molecules),
            "analysis": str((run_dir / "comparison.json").relative_to(args.output)),
            "license_note": (
                "No explicit reuse license identified in CCD source page"
                if name == "lj13" else
                "NOMAD source records declare CC BY 4.0; entry-level attribution retained in structures.csv"
            ),
            "conditions": len(report["conditions"]),
        })
    (args.output / "manifest.json").write_text(json.dumps({
        "schema_version": 1,
        "pyar_revision": __import__("subprocess").run(
            ["git", "rev-parse", "HEAD"], capture_output=True, text=True, check=False
        ).stdout.strip() or None,
        "manifests": manifests,
    }, indent=2) + "\n")
    print(f"Wrote atomic-cluster benchmark results to {args.output}")


if __name__ == "__main__":
    main()
