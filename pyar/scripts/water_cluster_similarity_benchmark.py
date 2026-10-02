"""Build a compact, fixed-composition water-cluster similarity benchmark.

The source is the W6, within-5-kcal/mol TTM2.1-F minima file from the
University of Washington water-cluster database. Sampling is deterministic
and stratified over its energy-sorted records. Similarity references are
pairwise fragment-matched, global-alignment RMSD upper bounds; no basin labels
are inferred from the source energy windows.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

from pyar.core.molecule import Molecule
from pyar.selection.distances import pairwise_distances
from pyar.selection.features import compute_feature_matrix, standardize_features
from pyar.structure_comparison.fragment_rmsd import FragmentRMSDComparator
from pyar.structure_comparison.graph_rmsd import infer_molecular_graph


FEATURES = ("mbtr", "soap", "distance-histogram")
SOURCE_URL = (
    "https://sites.uw.edu/wdbase/files/2019/01/"
    "W6_geoms_5.0_KCal-1hgztfv.txt"
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_water_minima(path: Path) -> list[dict]:
    """Read concatenated XYZ-like Wn_geoms records with Ord_Energy comments."""
    lines = path.read_text().splitlines()
    records = []
    cursor = 0
    while cursor < len(lines):
        while cursor < len(lines) and not lines[cursor].strip():
            cursor += 1
        if cursor >= len(lines):
            break
        try:
            atom_count = int(lines[cursor].strip())
        except ValueError as exc:
            raise ValueError(f"Expected atom count on source line {cursor + 1}") from exc
        if atom_count < 1 or cursor + atom_count + 1 >= len(lines):
            raise ValueError(f"Truncated geometry beginning on source line {cursor + 1}")
        comment = lines[cursor + 1].split()
        if len(comment) != 2 or comment[0] != "Ord_Energy":
            raise ValueError(f"Unexpected energy comment on source line {cursor + 2}")
        energy = float(comment[1])
        symbols, coordinates = [], []
        for line_number, line in enumerate(lines[cursor + 2:cursor + 2 + atom_count], cursor + 3):
            fields = line.split()
            if len(fields) != 4:
                raise ValueError(f"Malformed atom row on source line {line_number}")
            symbols.append(fields[0])
            coordinates.append([float(value) for value in fields[1:]])
        coordinates = np.asarray(coordinates, dtype=float)
        if symbols.count("O") != 6 or symbols.count("H") != 12 or len(symbols) != 18:
            raise ValueError(f"Source record {len(records)} is not an (H2O)6 cluster")
        if not np.isfinite(energy) or not np.isfinite(coordinates).all():
            raise ValueError(f"Source record {len(records)} contains non-finite values")
        record_index = len(records)
        records.append({
            "source_index": record_index,
            "energy_kcal_per_mol": energy,
            "molecule": Molecule(
                symbols, coordinates, name=f"W6_source_{record_index:04d}",
                energy=energy, charge=0, multiplicity=1,
            ),
        })
        cursor += atom_count + 2
    if not records:
        raise ValueError("No water-cluster structures found")
    energies = [record["energy_kcal_per_mol"] for record in records]
    if energies != sorted(energies):
        raise ValueError("Source structures are not sorted by ascending TTM2.1-F energy")
    return records


def stratified_sample(records: list[dict], sample_size: int) -> list[dict]:
    """Choose evenly spaced source ranks, preserving both energy-range ends."""
    if sample_size < 2 or sample_size > len(records):
        raise ValueError("sample_size must be at least 2 and no larger than the source set")
    indices = np.rint(np.linspace(0, len(records) - 1, sample_size)).astype(int)
    if len(set(indices.tolist())) != sample_size:
        raise ValueError("Sample size is too large for unique stratified source indices")
    return [records[int(index)] for index in indices]


def _lower_triangle(matrix):
    return np.asarray(matrix, dtype=float)[np.triu_indices(len(matrix), k=1)]


def evaluate_feature_distances(reference: np.ndarray, molecules: list[Molecule], features):
    """Compare feature-distance rankings with the independent RMSD ranking."""
    reference_values = _lower_triangle(reference)
    pair_count = len(reference_values)
    top_count = max(1, int(np.ceil(pair_count * 0.10)))
    reference_nearest = set(np.argsort(reference_values, kind="stable")[:top_count].tolist())
    results = []
    for feature in features:
        feature_result = compute_feature_matrix(
            molecules, feature, allow_fallbacks=True,
            system_type="molecular-aggregate", algorithm="agglomerative",
        )
        distances = pairwise_distances(
            standardize_features(feature_result.values), metric="euclidean",
        )
        feature_values = _lower_triangle(distances)
        feature_nearest = set(np.argsort(feature_values, kind="stable")[:top_count].tolist())
        rho_value = float(spearmanr(reference_values, feature_values).statistic)
        results.append({
            "feature_requested": feature,
            "feature_used": feature_result.name,
            "feature_fallbacks": list(feature_result.fallbacks),
            "spearman_rank_correlation_with_fragment_rmsd": (
                rho_value if np.isfinite(rho_value) else None
            ),
            "nearest_pair_overlap_at_10_percent": len(reference_nearest & feature_nearest),
            "nearest_pair_precision_at_10_percent": len(reference_nearest & feature_nearest) / top_count,
            "pairwise_distance_matrix": distances,
        })
    return results


def build_benchmark(source_path: Path, output_dir: Path, sample_size=24,
                    features=FEATURES, resume=True):
    records = load_water_minima(source_path)
    selected = stratified_sample(records, sample_size)
    molecules = [record["molecule"] for record in selected]
    output_dir.mkdir(parents=True, exist_ok=True)

    # Water-fragment topology is checked explicitly before invoking the
    # fragment matcher. A cluster with cross-monomer covalent bonds is rejected.
    import networkx as nx

    fragment_counts = [len(list(nx.connected_components(infer_molecular_graph(molecule))))
                       for molecule in molecules]
    if any(count != 6 for count in fragment_counts):
        raise ValueError(f"Expected six coordinate-inferred water fragments, found {fragment_counts}")

    pair_csv = output_dir / "reference_pairs.csv"
    cache_manifest = output_dir / "reference_pairs.cache.json"
    cache_identity = {
        "schema": 1,
        "source_sha256": sha256_file(source_path),
        "sample_size": sample_size,
        "selected_source_indices": [record["source_index"] for record in selected],
        "selected_geometry_sha256": hashlib.sha256(json.dumps([
            {"atoms": record["molecule"].atoms_list,
             "coordinates": np.asarray(record["molecule"].coordinates, dtype=float).tolist()}
            for record in selected
        ], sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest(),
        "implementation_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "comparator": {"name": "fragment-rmsd", "atom_mode": "all",
                       "max_mappings": 10000, "max_fragment_permutations": 720},
    }
    cached = {}
    if resume and pair_csv.is_file():
        with pair_csv.open(newline="") as stream:
            for row in csv.DictReader(stream):
                cached[(int(row["left_sample_index"]), int(row["right_sample_index"]))] = row
    expected_pairs = [(left, right) for left in range(sample_size)
                      for right in range(left + 1, sample_size)]
    expected_pair_set = set(expected_pairs)
    try:
        cached_identity = json.loads(cache_manifest.read_text())
    except (OSError, ValueError):
        cached_identity = None
    valid_cached = bool(resume and cached_identity == cache_identity and pair_csv.is_file())
    if valid_cached:
        try:
            valid_cached = set(cached).issubset(expected_pair_set) and all(
                cached[key]["left_source_index"] == str(selected[key[0]]["source_index"])
                and cached[key]["right_source_index"] == str(selected[key[1]]["source_index"])
                and cached[key]["comparison_complete"].lower() == "true"
                and np.isfinite(float(cached[key]["fragment_rmsd_upper_bound_angstrom"]))
                for key in cached
            )
        except (KeyError, TypeError, ValueError):
            valid_cached = False
    if not valid_cached:
        cached = {}
        pair_csv.write_text("left_sample_index,right_sample_index,left_source_index,right_source_index,fragment_rmsd_upper_bound_angstrom,comparison_complete,rotation_hypotheses\n")
        cache_manifest.write_text(json.dumps(cache_identity, indent=2, sort_keys=True) + "\n")
    comparator = FragmentRMSDComparator(
        atom_mode="all", max_mappings=10000, max_fragment_permutations=720,
    )
    completed_pairs = set(cached) if valid_cached else set()
    with pair_csv.open("a", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        for left, right in expected_pairs:
            if (left, right) in completed_pairs:
                continue
            comparison = comparator.compare(molecules[left], molecules[right])
            if not comparison.compatible or comparison.distance is None:
                raise RuntimeError(
                    f"Could not compare structures {left} and {right}: {comparison.metadata}"
                )
            if comparison.metadata.get("comparison_complete") is not True:
                raise RuntimeError(
                    f"Reference fragment comparison was incomplete for pair {left},{right}"
                )
            writer.writerow([
                left, right, selected[left]["source_index"], selected[right]["source_index"],
                f"{comparison.distance:.12g}", True,
                comparison.metadata.get("rotation_hypotheses", ""),
            ])
            stream.flush()
            print(f"reference pairs: {len(completed_pairs) + 1}/{len(expected_pairs)}", flush=True)
            completed_pairs.add((left, right))

    distances = np.zeros((sample_size, sample_size), dtype=float)
    with pair_csv.open(newline="") as stream:
        for row in csv.DictReader(stream):
            left, right = int(row["left_sample_index"]), int(row["right_sample_index"])
            distances[left, right] = distances[right, left] = float(
                row["fragment_rmsd_upper_bound_angstrom"]
            )
    np.savetxt(output_dir / "reference_fragment_rmsd_upper_bound_angstrom.csv",
               distances, delimiter=",")
    reference_values = _lower_triangle(distances)

    feature_results = evaluate_feature_distances(distances, molecules, features)
    feature_report = []
    for result in feature_results:
        matrix = result.pop("pairwise_distance_matrix")
        filename = f"{result['feature_requested']}_euclidean_feature_distances.csv"
        np.savetxt(output_dir / filename, matrix, delimiter=",")
        result["distance_matrix_file"] = filename
        feature_report.append(result)

    structures = []
    xyz_path = output_dir / "structures.xyz"
    with xyz_path.open("w") as stream:
        for sample_index, record in enumerate(selected):
            molecule = record["molecule"]
            structures.append({
                "sample_index": sample_index,
                "source_index": record["source_index"],
                "energy_kcal_per_mol": record["energy_kcal_per_mol"],
                "source_file": source_path.name,
            })
            stream.write("18\n")
            stream.write(
                f"sample={sample_index} source_index={record['source_index']} "
                f"ttm21f_energy_kcal_mol={record['energy_kcal_per_mol']:.10f}\n"
            )
            for symbol, coordinate in zip(molecule.atoms_list, molecule.coordinates):
                stream.write(f"{symbol} {coordinate[0]:.12f} {coordinate[1]:.12f} {coordinate[2]:.12f}\n")
    (output_dir / "structures.csv").write_text(
        "sample_index,source_index,energy_kcal_per_mol,source_file\n"
        + "".join(
            f"{row['sample_index']},{row['source_index']},{row['energy_kcal_per_mol']:.10f},{row['source_file']}\n"
            for row in structures
        )
    )
    report = {
        "schema_version": 1,
        "dataset": "W6 water-cluster structural-similarity subset",
        "source_url": SOURCE_URL,
        "source_file": source_path.name,
        "source_sha256": sha256_file(source_path),
        "source_geometry_count": len(records),
        "sampled_geometry_count": sample_size,
        "sample_strategy": "evenly spaced source ranks over the energy-sorted within-5-kcal/mol source subset",
        "composition": "(H2O)6; every geometry is kept in one fixed-composition pool",
        "coordinate_unit": "angstrom",
        "energy_unit": "kcal/mol; source TTM2.1-F potential energy",
        "reference_similarity": {
            "method": "fragment-matched global Kabsch RMSD upper bound",
            "implementation": "PyAR FragmentRMSDComparator",
            "atom_mode": "all atoms; six water fragments can permute, and water hydrogens can permute within each fragment",
            "fragment_validation": "six coordinate-inferred disconnected H2O graphs per structure",
            "mapping_limits": {"max_fragment_permutations": 720, "max_mappings": 10000},
            "threshold_labels": "none; continuous pairwise distances are retained to avoid inventing class labels",
            "comparison_pairs": len(expected_pairs),
            "all_comparisons_complete": True,
            "distance_quantiles_angstrom": {
                str(percentile): float(np.percentile(reference_values, percentile))
                for percentile in (0, 10, 25, 50, 75, 90, 100)
            },
        },
        "distance_matrix_row_order": "sample_index order in structures.csv and structures.xyz",
        "structures": structures,
        "feature_distance_comparison": feature_report,
        "limitations": [
            "Only one composition and one cluster size are represented.",
            "The reference comparator reports a mapped RMSD upper bound, not a guaranteed global minimum over all possible mappings.",
            "The source structures are distinct TTM2.1-F local minima; energy-window membership is not a basin-equivalence label.",
            "The source page does not state a separate reuse license for the per-size file; attribution and original URL are retained.",
        ],
    }
    (output_dir / "manifest.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"Wrote {sample_size} structures and {len(expected_pairs)} pair distances to {output_dir}")
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--sample-size", type=int, default=24)
    parser.add_argument("--features", nargs="+", choices=FEATURES, default=list(FEATURES))
    parser.add_argument("--no-resume", action="store_true")
    args = parser.parse_args(argv)
    build_benchmark(args.source, args.output, args.sample_size, args.features,
                    resume=not args.no_resume)


if __name__ == "__main__":
    main()
