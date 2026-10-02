"""Audit and compare clustering distances on labelled geometry-family pools."""

from __future__ import annotations

import argparse
import hashlib
from importlib.metadata import PackageNotFoundError, version
import json
import time
from collections import Counter
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from pyar.selection.clusterers import cluster_molecules
from pyar.selection.distances import pairwise_distances
from pyar.selection.features import standardize_features
from pyar.selection.policy import SELECTION_POLICY_VERSION
from pyar.structure_comparison.coordinate_graph import infer_coordinate_graph


def load_dataset(manifest_path):
    """Validate stored geometries and independent transformation witnesses."""
    manifest_path = Path(manifest_path)
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("schema_version") != 1:
        raise ValueError("Unsupported clustering-distance dataset schema")
    xyz = manifest_path.parent / manifest["structures_file"]
    lines = xyz.read_text().splitlines()
    frames, index = {}, 0
    while index < len(lines):
        natoms = int(lines[index])
        identifier = lines[index + 1].removeprefix("structure_id=")
        rows = [line.split() for line in lines[index + 2:index + 2 + natoms]]
        coordinates = np.asarray([[float(value) for value in row[1:]] for row in rows])
        if natoms < 1 or len(rows) != natoms or coordinates.shape != (natoms, 3) or not np.isfinite(coordinates).all():
            raise ValueError(f"Invalid geometry {identifier}")
        if identifier in frames:
            raise ValueError(f"Repeated geometry id {identifier}")
        frames[identifier] = SimpleNamespace(name=identifier, atoms_list=[row[0] for row in rows],
                                             coordinates=coordinates)
        index += natoms + 2
    referenced = set()
    for pool in manifest["pools"]:
        identifiers = [record["id"] for record in pool["records"]]
        if len(set(identifiers)) != len(identifiers) or any(identifier not in frames for identifier in identifiers):
            raise ValueError("Pool contains duplicate or unavailable geometries")
        compositions = {tuple(sorted(Counter(frames[identifier].atoms_list).items())) for identifier in identifiers}
        if len(compositions) != 1:
            raise ValueError("Benchmark pool must have one formula")
        referenced.update(identifiers)
        # A different sorted pair-distance spectrum independently proves that
        # cross-family structures cannot be related by a rigid atom permutation.
        spectra = {}
        for identifier in identifiers:
            coordinates = frames[identifier].coordinates
            distances = np.linalg.norm(coordinates[:, None] - coordinates[None, :], axis=2)
            spectra[identifier] = np.sort(distances[np.triu_indices(len(coordinates), 1)])
        for left, record in enumerate(pool["records"]):
            for other in pool["records"][left + 1:]:
                if record["family"] != other["family"] and np.allclose(
                    spectra[record["id"]], spectra[other["id"]], rtol=0, atol=1e-6
                ):
                    raise ValueError("Different fixture families lack an independent noncongruence witness")
        for record in pool["records"]:
            if not np.isfinite(record["energy"]):
                raise ValueError("Ranking energies must be finite")
            frames[record["id"]].energy = record["energy"]
        for witness in pool["rigid_transform_witnesses"]:
            if witness["first"] not in identifiers or witness["second"] not in identifiers:
                raise ValueError("Transformation witness is outside its pool")
            first, second = frames[witness["first"]], frames[witness["second"]]
            rotation = np.asarray(witness["rotation"])
            order = np.asarray(witness["order"])
            if not np.allclose(rotation.T @ rotation, np.eye(3), atol=1e-10) or not np.isclose(np.linalg.det(rotation), 1):
                raise ValueError("Witness must specify a proper rigid rotation")
            if sorted(order.tolist()) != list(range(len(first.atoms_list))):
                raise ValueError("Witness atom order must be a permutation")
            if np.asarray(first.atoms_list)[order].tolist() != second.atoms_list:
                raise ValueError("Witness permutation does not preserve element labels")
            if not np.allclose(first.coordinates[order] @ rotation + witness["translation"],
                               second.coordinates, atol=2e-10, rtol=0):
                raise ValueError("Stored geometry violates its independent transformation witness")
            families = {record["id"]: record["family"] for record in pool["records"]}
            if families[witness["first"]] != families[witness["second"]]:
                raise ValueError("Rigid copies must have the same family label")
        if pool["system_type"] == "molecular-aggregate":
            import networkx as nx
            for identifier in identifiers:
                molecule = frames[identifier]
                graph = infer_coordinate_graph(molecule)
                parts = list(nx.connected_components(graph))
                if len(parts) < 2 or any(len(part) < 2 for part in parts):
                    raise ValueError("Aggregate fixture lost its molecular fragments")
                for left, part in enumerate(parts):
                    for other in parts[left + 1:]:
                        distances = np.linalg.norm(molecule.coordinates[list(part)][:, None] - molecule.coordinates[list(other)], axis=2)
                        if np.min(distances) < 1.2:
                            raise ValueError("Aggregate fixture has pathological intermolecular contacts")
    if referenced != set(frames):
        raise ValueError("Dataset includes unreferenced geometries")
    audit = {"geometry_count": len(frames), "pool_count": len(manifest["pools"]),
             "witness_count": sum(len(pool["rigid_transform_witnesses"]) for pool in manifest["pools"]),
             "manifest_sha256": hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
             "structures_sha256": hashlib.sha256(xyz.read_bytes()).hexdigest(),
             "label_interpretation": manifest["label_interpretation"]}
    return manifest, frames, audit


def run_benchmark(manifest_path, *, output=None):
    """Compare baselines and structural matrices using the actual clusterer."""
    from sklearn.metrics import adjusted_rand_score

    manifest, frames, audit = load_dataset(manifest_path)
    rows = []
    configurations = [("mbtr-euclidean", "mbtr", "euclidean", {}),
                      ("soap-mean-euclidean", "soap", "euclidean", {}),
                      ("histogram-euclidean", "distance-histogram", "euclidean", {}),
                      ("graph-rmsd-heavy", "auto", "graph-rmsd", {}),
                      ("graph-rmsd-all", "auto", "graph-rmsd", {"atom_mode": "all"}),
                      ("fragment-rmsd", "auto", "fragment-rmsd", {}),
                      ("soap-rematch", "auto", "soap-rematch", {})]
    for pool in manifest["pools"]:
        molecules = [frames[record["id"]] for record in pool["records"]]
        reference_labels = [record["family"] for record in pool["records"]]
        indices = {molecule.name: index for index, molecule in enumerate(molecules)}
        for name, feature, distance, options in configurations:
            started = time.perf_counter()
            result = cluster_molecules(
                molecules, feature=feature, distance_metric=distance, distance_options=options,
                algorithm="agglomerative", maximum_number_of_clusters=len(set(reference_labels)),
                system_type=pool["system_type"],
            )
            matrix = result.distance_matrix
            if matrix is None:
                matrix = pairwise_distances(standardize_features(result.feature_values), result.distance_metric)
            errors = [matrix[indices[witness["first"]], indices[witness["second"]]]
                      for witness in pool["rigid_transform_witnesses"]]
            rows.append({"pool": pool["id"], "configuration": name,
                         "adjusted_rand_index": float(adjusted_rand_score(reference_labels, result.labels)),
                         "maximum_rigid_copy_distance": float(max(errors, default=0)),
                         "runtime_seconds": time.perf_counter() - started,
                         **result.to_dict()})
    versions = {}
    for package in ("numpy", "scipy", "dscribe", "scikit-learn", "rdkit"):
        try:
            versions[package] = version(package)
        except PackageNotFoundError:
            versions[package] = None
    code_files = [Path(__file__), Path(__file__).parents[1] / "selection" / "structural_distances.py",
                  Path(__file__).parents[1] / "structure_comparison" / "fragment_rmsd.py",
                  Path(__file__).parents[1] / "selection" / "features.py",
                  Path(__file__).parents[1] / "selection" / "clusterers.py"]
    code_files.extend(Path(__file__).parents[1] / "structure_comparison" / name
                      for name in ("graph_rmsd.py", "deduplication_policy.py", "rmsd.py", "coordinate_graph.py"))
    code_files.append(Path(__file__).parents[1] / "selection" / "policy.py")
    report = {"audit": audit, "results": rows, "package_versions": versions,
              "selection_policy_version": SELECTION_POLICY_VERSION,
              "code_sha256": {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in code_files},
              "limitations": ["Constructed geometric families are not PES basin labels.",
                              "No production threshold or universal default is calibrated here.",
                              "REMatch is a local, finite-cutoff, reflection-invariant descriptor comparison.",
                              "Fragment RMSD is a symmetric conservative upper bound, not an exact global minimum."]}
    if output:
        Path(output).write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest")
    parser.add_argument("--output", default="distance_benchmark.json")
    parser.add_argument("--audit-only", action="store_true")
    args = parser.parse_args(argv)
    if args.audit_only:
        _, _, audit = load_dataset(args.manifest)
        print(json.dumps(audit, indent=2))
    else:
        report = run_benchmark(args.manifest, output=args.output)
        print(f"Validated {report['audit']['geometry_count']} geometries; {len(report['results'])} comparisons")
        print(f"Report: {args.output}")


if __name__ == "__main__":
    main()
