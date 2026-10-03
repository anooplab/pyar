"""Reproducible, basin-labelled clustering benchmark for conformer pools.

The benchmark separates reference-label generation (independent xTB local
optimization plus graph-constrained RMSD) from seed selection. Every input,
optimization, and report is retained below the chosen benchmark run directory.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import shutil
import subprocess
import tempfile
import time
from copy import deepcopy
from pathlib import Path

import numpy as np
import networkx as nx
from ase.io import read, write

from pyar.core.molecule import Molecule
from pyar.selection import clustering
from pyar.selection.diversity import _max_min_diversity_select
from pyar.selection.features import standardize_features
from pyar.structure_comparison.graph_rmsd import GraphRMSDComparator, infer_molecular_graph
from pyar.structure_comparison.rmsd import kabsch_rmsd


ALGORITHMS = ("auto", "agglomerative", "dbscan", "optics", "maxmin")
FEATURES = ("mbtr", "soap", "distance-histogram")
ENERGY_RE = re.compile(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[Ee][-+]?\d+)?")
TOTAL_ENERGY_RE = re.compile(r"TOTAL ENERGY\s+([-+0-9.Ee]+)")
GRADIENT_NORM_RE = re.compile(r"GRADIENT NORM\s+([-+0-9.Ee]+)")
NORMAL_TERMINATION_RE = re.compile(r"(?im)^\s*normal termination of xtb\s*$")
MAX_REFERENCE_GRADIENT_NORM = 1.0e-3


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _normally_terminated(log: str) -> bool:
    """Avoid accepting xTB's distinct phrase 'abnormal termination of xtb'."""
    return NORMAL_TERMINATION_RE.search(log) is not None


def _read_energy(atoms) -> float:
    comment = atoms.info.get("comment", "")
    if not comment:
        # ASE's extxyz parser treats a bare numeric comment as a key.
        comment = next(iter(atoms.info), "")
    match = ENERGY_RE.search(str(comment))
    if not match:
        raise ValueError(f"No numeric source energy in XYZ comment: {comment!r}")
    value = float(match.group())
    if not np.isfinite(value):
        raise ValueError("Source energy is not finite")
    return value


def load_ensemble(path: Path):
    frames = read(str(path), index=":")
    if not frames:
        raise ValueError(f"No geometries found in {path}")
    symbols = frames[0].get_chemical_symbols()
    composition = sorted((symbol, symbols.count(symbol)) for symbol in set(symbols))
    result = []
    for index, atoms in enumerate(frames):
        if atoms.get_chemical_symbols() != symbols:
            raise ValueError(f"Frame {index} does not preserve the atom ordering")
        if not np.isfinite(atoms.positions).all():
            raise ValueError(f"Frame {index} contains non-finite coordinates")
        result.append({
            "index": index,
            "name": f"frame_{index:04d}",
            "atoms": atoms,
            "source_energy": _read_energy(atoms),
            "composition": composition,
        })
    return result


def _write_xyz(path: Path, atoms, title: str):
    path.parent.mkdir(parents=True, exist_ok=True)
    write(str(path), atoms, format="xyz", comment=title)


def optimize_ensemble(ensemble, run_dir: Path, xtb: str, limit: int | None = None):
    executable = shutil.which(xtb) if "/" not in xtb else xtb
    if not executable or not Path(executable).is_file():
        raise RuntimeError(f"Cannot find xTB executable {xtb!r}")
    executable = str(Path(executable).resolve())
    records = []
    subset = ensemble if limit is None else ensemble[:limit]
    executable_hash = _sha256(Path(executable))
    optimization_settings = {
        "schema": 1, "executable": str(Path(executable).resolve()),
        "executable_sha256": executable_hash,
        "arguments": ["--opt", "--gfn", "2", "--chrg", "0", "--uhf", "0"],
    }
    for record in subset:
        frame_dir = run_dir / "optimizations" / record["name"]
        output_path = frame_dir / "xtbopt.xyz"
        log_path = frame_dir / "xtb.log"
        cache_path = frame_dir / "cache.json"
        input_identity = hashlib.sha256(json.dumps({
            "symbols": record["atoms"].get_chemical_symbols(),
            "coordinates": np.asarray(record["atoms"].positions, dtype=float).tolist(),
            "settings": optimization_settings,
        }, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()
        cache_valid = False
        if output_path.is_file() and log_path.is_file() and cache_path.is_file():
            try:
                cache = json.loads(cache_path.read_text())
                cache_valid = (
                    cache.get("input_identity") == input_identity
                    and cache.get("output_sha256") == _sha256(output_path)
                    and cache.get("log_sha256") == _sha256(log_path)
                    and cache.get("settings") == optimization_settings
                    and _normally_terminated(log_path.read_text(errors="replace"))
                )
            except (OSError, ValueError, TypeError):
                cache_valid = False
        if not cache_valid:
            frame_dir.mkdir(parents=True, exist_ok=True)
            # Never let old xtbopt.xyz or xtbrestart files satisfy a new run.
            # Keep each attempt for inspection, including failures.
            attempt = Path(tempfile.mkdtemp(prefix="attempt_", dir=frame_dir))
            input_path = attempt / "input.xyz"
            _write_xyz(input_path, record["atoms"], record["name"])
            completed = subprocess.run(
                [executable, input_path.name, *optimization_settings["arguments"]],
                cwd=attempt, capture_output=True, text=True, timeout=1800, check=False,
            )
            attempt_log = attempt / "xtb.log"
            attempt_log.write_text(completed.stdout + "\n--- STDERR ---\n" + completed.stderr)
            if completed.returncode != 0:
                raise RuntimeError(
                    f"xTB failed for {record['name']} (exit {completed.returncode}); see {attempt_log}"
                )
            if not (attempt / "xtbopt.xyz").is_file() or not _normally_terminated(attempt_log.read_text(errors="replace")):
                raise RuntimeError(f"xTB optimization did not terminate normally for {record['name']}")
            for name in ("input.xyz", "xtbopt.xyz", "xtb.log"):
                shutil.copy2(attempt / name, frame_dir / name)
            cache_path.write_text(json.dumps({
                "schema": 1, "input_identity": input_identity,
                "settings": optimization_settings, "output_sha256": _sha256(output_path),
                "log_sha256": _sha256(log_path),
            }, indent=2, sort_keys=True) + "\n")
        log = log_path.read_text(errors="replace")
        if not _normally_terminated(log) or not output_path.is_file():
            raise RuntimeError(f"xTB optimization did not terminate normally for {record['name']}")
        energy_matches = TOTAL_ENERGY_RE.findall(log)
        if not energy_matches or not np.isfinite(float(energy_matches[-1])):
            raise RuntimeError(f"No final TOTAL ENERGY found in xTB log for {record['name']}")
        gradient_matches = GRADIENT_NORM_RE.findall(log)
        if not gradient_matches:
            raise RuntimeError(f"No final gradient norm found in xTB log for {record['name']}")
        final_gradient_norm = float(gradient_matches[-1])
        if not np.isfinite(final_gradient_norm) or final_gradient_norm > MAX_REFERENCE_GRADIENT_NORM:
            raise RuntimeError(
                f"xTB geometry for {record['name']} is not sufficiently stationary "
                f"(gradient norm {final_gradient_norm:.3g}); see {log_path}"
            )
        optimized = read(str(output_path))
        if not np.isfinite(optimized.positions).all():
            raise ValueError(f"Nonfinite optimized coordinates for {record['name']}")
        if optimized.get_chemical_symbols() != record["atoms"].get_chemical_symbols():
            raise ValueError(f"xTB changed atom identities/order in {record['name']}")
        records.append({
            **record,
            "optimized_atoms": optimized,
            "xtb_energy_hartree": float(energy_matches[-1]),
            "final_gradient_norm_eh_per_alpha": final_gradient_norm,
            "optimization_log": str(log_path),
        })
        print(f"optimized {len(records)}/{len(subset)}: {record['name']}", flush=True)
    return records


def assign_reference_basins(records, threshold: float = 0.35):
    """Greedily assign to first verified optimized reference within graph RMSD."""
    comparator = GraphRMSDComparator(
        threshold=threshold, max_isomorphisms=20000, atom_mode="heavy",
    )
    representatives = []
    for record in records:
        molecule = _as_molecule(record["optimized_atoms"], record["name"])
        graph = infer_molecular_graph(molecule)
        record["inferred_topology_hash"] = nx.weisfeiler_lehman_graph_hash(
            graph, node_attr="element",
        )
        assigned = None
        best_distance = None
        for basin_id, representative in enumerate(representatives):
            result = comparator.compare(molecule, representative["molecule"])
            if (result.metadata.get("comparison_complete") is True
                    and result.equivalent is True):
                assigned = basin_id
                best_distance = result.distance
                break
        if assigned is None:
            assigned = len(representatives)
            representatives.append({"molecule": molecule, "record": record})
        record["basin_id"] = assigned
        record["basin_match_rmsd_angstrom"] = best_distance
    return records, representatives


def compare_perturbed_seeds(records, threshold: float):
    """Test whether each perturbed/reoptimized seed returns to its own basin."""
    comparator = GraphRMSDComparator(
        threshold=threshold, max_isomorphisms=20000, atom_mode="heavy",
    )
    for record in records:
        molecule = _as_molecule(record["optimized_atoms"], record["name"])
        original = _as_molecule(record["reference_atoms"], record["name"])
        result = comparator.compare(molecule, original)
        record["reference_basin_match_status"] = (
            "preserved" if result.equivalent is True else
            "uncertain-incomplete-mapping" if result.metadata.get("comparison_complete") is not True else
            "changed-complete-comparison"
        )
        record["reference_basin_rmsd_angstrom"] = result.distance
    return records


def _as_molecule(atoms, name, energy=None):
    return Molecule(
        atoms.get_chemical_symbols(), np.asarray(atoms.positions, dtype=float),
        name=name, energy=energy, charge=0, multiplicity=1,
    )


def _coverage(records, selected_records, threshold):
    basin_representatives = {}
    for record in records:
        basin_representatives.setdefault(record["basin_id"], record)
    selected_atoms = [item["optimized_atoms"] for item in selected_records]
    distances = []
    for record in basin_representatives.values():
        target = record["optimized_atoms"]
        target_symbols = target.get_chemical_symbols()
        heavy = [i for i, symbol in enumerate(target_symbols) if symbol != "H"]
        candidates = []
        for selected in selected_atoms:
            if selected.get_chemical_symbols() != target_symbols:
                continue
            candidates.append(float(kabsch_rmsd(target.positions[heavy], selected.positions[heavy])))
        distances.append(min(candidates) if candidates else None)
    finite = [value for value in distances if value is not None and np.isfinite(value)]
    covered = sum(value is not None and value < threshold for value in distances)
    return {
        "distance_method": "atom-order-corresponding heavy-atom Kabsch RMSD; the upstream XYZ ensemble preserves a common atom order",
        "reference_basin_count": len(basin_representatives),
        "covered_basins_at_label_threshold": covered,
        "basin_recall": covered / len(distances),
        "mean_nearest_indexed_heavy_atom_rmsd_angstrom": float(np.mean(finite)) if finite else None,
        "max_nearest_indexed_heavy_atom_rmsd_angstrom": float(max(finite)) if finite else None,
    }


def _cluster_first_select(molecules, algorithm, feature, maximum_seeds,
                          system_type="conformers"):
    """Run the cluster/minima/budget stages without the separate deduplicator."""
    cluster_algorithm = "auto" if algorithm == "maxmin" else algorithm
    result = clustering.cluster_molecules(
        molecules,
        feature=feature,
        algorithm=cluster_algorithm,
        maximum_number_of_clusters=maximum_seeds,
        distance_metric="euclidean",
        system_type=system_type,
    )
    labels = np.asarray(result.labels, dtype=int)
    minima = []
    for label in sorted(set(labels.tolist())):
        members = [i for i, value in enumerate(labels) if value == label]
        if label == -1:
            minima.extend(members)
        elif members:
            minima.append(min(members, key=lambda index: float(molecules[index].energy)))
    candidate_molecules = [molecules[index] for index in minima]
    if len(candidate_molecules) > maximum_seeds:
        standardized = standardize_features(result.feature_values)
        candidate_features = standardized[minima]
        candidate_molecules = _max_min_diversity_select(
            candidate_features, candidate_molecules, maximum_seeds,
            distance_metric="euclidean",
        )
    diagnostics = result.to_dict()
    diagnostics.update({
        "selection_algorithm_requested": algorithm,
        "selection_stage": "cluster-minima-then-maxmin-budget-trim",
        "deduplication_stage": "excluded-to-isolate-clustering-and-feature-quality",
        "cluster_minimum_count": len(minima),
        "selected_count": len(candidate_molecules),
        "selected_names": [item.name for item in candidate_molecules],
    })
    return candidate_molecules, diagnostics


def run_conditions(records, representatives, output_dir: Path, algorithms, features,
                   maximum_seeds: int, basin_threshold: float, downstream_optimize: bool,
                   downstream_perturbation: float, xtb: str):
    out = output_dir / "results"
    out.mkdir(parents=True, exist_ok=True)
    report = {
        "schema_version": 1,
        "benchmark": output_dir.name,
        "input_count": len(records),
        "reference_basin_count": len(representatives),
        "reference_label": {
            "method": "independent GFN2-xTB local optimization followed by complete element-labelled coordinate-graph isomorphism plus heavy-atom Kabsch RMSD",
            "rmsd_threshold_angstrom": basin_threshold,
            "atom_mode": "heavy",
            "terminal_hydrogen_connectivity_and_counts_checked_but_hydrogen_coordinates_excluded": True,
            "incomplete_graph_comparisons": "never merged; retained as separate labels",
            "charge": 0,
            "multiplicity": 1,
        },
        "seed_budget": maximum_seeds,
        "conditions": [],
    }
    downstream_cache = {}
    checkpoint = out / "comparison.partial.json"
    for algorithm in algorithms:
        for feature in features:
            condition_id = f"{algorithm}__{feature}"
            molecules = [
                _as_molecule(item["atoms"], item["name"], item["xtb_energy_hartree"])
                for item in records
            ]
            started = time.perf_counter()
            selected, diagnostics = _cluster_first_select(
                molecules, algorithm, feature, maximum_seeds,
            )
            elapsed = time.perf_counter() - started
            by_name = {item["name"]: item for item in records}
            selected_records = [by_name[item.name] for item in selected]
            selected_basin_ids = sorted({item["basin_id"] for item in selected_records})
            condition = {
                "algorithm_requested": algorithm,
                "feature_requested": feature,
                "runtime_seconds": elapsed,
                "selected_count": len(selected_records),
                "selected_names": [item["name"] for item in selected_records],
                "selected_reference_basin_ids": selected_basin_ids,
                "basins_retained": len(selected_basin_ids),
                "retention_fraction": len(selected_basin_ids) / len(representatives),
                "coverage": _coverage(records, selected_records, basin_threshold),
                "diagnostics": diagnostics,
                "downstream_optimization": None,
            }
            if downstream_optimize:
                perturbed_selected = []
                for selected_record in selected_records:
                    if selected_record["name"] in downstream_cache:
                        continue
                    perturbed = deepcopy(selected_record)
                    perturbed["atoms"] = selected_record["optimized_atoms"].copy()
                    perturbed["reference_atoms"] = selected_record["optimized_atoms"].copy()
                    rng_seed = selected_record["index"] + 90210
                    rng = np.random.default_rng(rng_seed)
                    perturbed["atoms"].positions += rng.normal(
                        0.0, downstream_perturbation, size=perturbed["atoms"].positions.shape,
                    )
                    perturbed_selected.append(perturbed)
                optimized_selected = optimize_ensemble(
                    perturbed_selected,
                    output_dir / "downstream" / "shared",
                    xtb,
                )
                optimized_selected = compare_perturbed_seeds(
                    optimized_selected, basin_threshold,
                )
                downstream_cache.update({item["name"]: item for item in optimized_selected})
                optimized_selected = [downstream_cache[item["name"]] for item in selected_records]
                downstream_valid = sum(
                    item.get("xtb_energy_hartree") is not None for item in optimized_selected
                )
                matched_records = [item for item in optimized_selected
                                   if item["reference_basin_match_status"] == "preserved"]
                recovered_ids = {item["basin_id"] for item in matched_records}
                recovered = len(recovered_ids)
                changed_complete = sum(
                    item["reference_basin_match_status"] == "changed-complete-comparison"
                    for item in optimized_selected
                )
                uncertain_mapping = sum(
                    item["reference_basin_match_status"] == "uncertain-incomplete-mapping"
                    for item in optimized_selected
                )
                condition["downstream_optimization"] = {
                    "successful_terminations": downstream_valid,
                    "selected_count": len(selected_records),
                    "success_fraction": downstream_valid / len(selected_records) if selected_records else 0.0,
                    "distinct_reoptimized_basins": recovered,
                    "verified_reference_basin_seed_count": len(matched_records),
                    "verified_reference_basin_recovery_fraction": (
                        len(matched_records) / len(optimized_selected) if optimized_selected else 0.0
                    ),
                    "basin_changed_after_complete_comparison_count": changed_complete,
                    "uncertain_due_to_incomplete_mapping_count": uncertain_mapping,
                    "selected_basins_retained_after_reoptimization": recovered,
                    "coordinate_perturbation_sigma_angstrom": downstream_perturbation,
                    "perturbation_seed_policy": "input frame index + 90210",
                    "unique_geometries_optimized_and_cached": len(downstream_cache),
                }
            report["conditions"].append(condition)
            checkpoint_tmp = checkpoint.with_suffix(".json.tmp")
            checkpoint_tmp.write_text(json.dumps(report, indent=2) + "\n")
            checkpoint_tmp.replace(checkpoint)
            print(
                f"{condition_id}: selected {len(selected_records)}, retained "
                f"{len(selected_basin_ids)}/{len(representatives)} basins", flush=True,
            )
    (out / "comparison.json").write_text(json.dumps(report, indent=2) + "\n")
    checkpoint.unlink(missing_ok=True)
    with (out / "comparison.csv").open("w", newline="") as stream:
        columns = ["algorithm_requested", "feature_requested", "selected_count",
                   "basins_retained", "reference_basin_count", "retention_fraction",
                   "mean_nearest_indexed_heavy_atom_rmsd_angstrom",
                   "max_nearest_indexed_heavy_atom_rmsd_angstrom",
                   "runtime_seconds", "algorithm_used", "feature_used"]
        writer = csv.DictWriter(stream, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        for row in report["conditions"]:
            writer.writerow({
                "algorithm_requested": row["algorithm_requested"],
                "feature_requested": row["feature_requested"],
                "selected_count": row["selected_count"],
                "basins_retained": row["basins_retained"],
                "reference_basin_count": report["reference_basin_count"],
                "retention_fraction": row["retention_fraction"],
                "mean_nearest_indexed_heavy_atom_rmsd_angstrom": row["coverage"]["mean_nearest_indexed_heavy_atom_rmsd_angstrom"],
                "max_nearest_indexed_heavy_atom_rmsd_angstrom": row["coverage"]["max_nearest_indexed_heavy_atom_rmsd_angstrom"],
                "runtime_seconds": row["runtime_seconds"],
                "algorithm_used": row["diagnostics"].get("algorithm_used"),
                "feature_used": row["diagnostics"].get("feature_used"),
            })
    return report


def write_reference_labels(records, representatives, run_dir: Path, source_path: Path, threshold: float):
    target = run_dir / "reference_basins.json"
    sensitivity = {}
    for cutoff in sorted({0.25, threshold, 0.50}):
        copied = [dict(record) for record in records]
        _, cutoff_representatives = assign_reference_basins(copied, cutoff)
        sensitivity[f"{cutoff:.2f}"] = {
            "basin_count": len(cutoff_representatives),
            "largest_basin_size": max(
                (sum(row["basin_id"] == basin for row in copied)
                 for basin in range(len(cutoff_representatives))),
                default=0,
            ),
        }
    payload = {
        "schema_version": 1,
        "source_file": source_path.name,
        "source_sha256": _sha256(source_path),
        "label_method": "GFN2-xTB opt then element-labelled coordinate-graph isomorphism and heavy-atom Kabsch RMSD",
        "threshold_angstrom": threshold,
        "threshold_sensitivity": sensitivity,
        "inferred_topology_groups": {
            topology: sum(row.get("inferred_topology_hash") == topology for row in records)
            for topology in sorted({row.get("inferred_topology_hash") for row in records})
        },
        "convergence_check": {
            "normal_termination_required": True,
            "finite_final_energy_and_gradient_required": True,
            "maximum_final_gradient_norm_eh_per_alpha": MAX_REFERENCE_GRADIENT_NORM,
        },
        "basins": [
            {"basin_id": i, "representative": entry["record"]["name"],
             "members": [r["name"] for r in records if r["basin_id"] == i]}
            for i, entry in enumerate(representatives)
        ],
        "structures": [
            {"name": r["name"], "basin_id": r["basin_id"],
             "source_energy": r["source_energy"],
             "source_energy_interpretation": "untrusted upstream XYZ comment; see provenance note",
             "inferred_topology_hash": r["inferred_topology_hash"],
             "xtb_energy_hartree": r["xtb_energy_hartree"],
             "final_gradient_norm_eh_per_alpha": r["final_gradient_norm_eh_per_alpha"],
             "optimized_xyz": str(Path("optimizations") / r["name"] / "xtbopt.xyz"),
             "optimization_log": str(Path("optimizations") / r["name"] / "xtb.log")}
            for r in records
        ],
    }
    target.write_text(json.dumps(payload, indent=2) + "\n")
    return target


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("ensemble", type=Path, help="Multi-frame XYZ conformer ensemble")
    parser.add_argument("--output", type=Path, required=True, help="Persistent run directory")
    parser.add_argument("--xtb", default="xtb", help="GFN2-xTB executable")
    parser.add_argument("--limit", type=int, help="Debug a leading subset; not suitable for final results")
    parser.add_argument("--basin-threshold", type=float, default=0.35, help="Graph RMSD basin threshold in Angstrom")
    parser.add_argument("--max-seeds", type=int, default=12)
    parser.add_argument("--downstream-perturbation", type=float, default=0.10,
                        help="Gaussian per-coordinate sigma in Angstrom before downstream optimization")
    parser.add_argument("--algorithms", nargs="+", choices=ALGORITHMS, default=list(ALGORITHMS))
    parser.add_argument("--features", nargs="+", choices=FEATURES, default=list(FEATURES))
    parser.add_argument("--skip-downstream-optimization", action="store_true")
    args = parser.parse_args(argv)
    if args.limit is not None and args.limit < 1:
        parser.error("--limit must be positive")
    if args.basin_threshold <= 0 or args.max_seeds < 1 or args.downstream_perturbation < 0:
        parser.error("basin threshold/seed budget must be positive and perturbation nonnegative")
    args.output.mkdir(parents=True, exist_ok=True)
    ensemble = load_ensemble(args.ensemble)
    optimized = optimize_ensemble(ensemble, args.output, args.xtb, limit=args.limit)
    optimized, representatives = assign_reference_basins(optimized, args.basin_threshold)
    label_path = write_reference_labels(optimized, representatives, args.output,
                                        args.ensemble, args.basin_threshold)
    report = run_conditions(
        optimized, representatives, args.output, args.algorithms, args.features,
        args.max_seeds, args.basin_threshold, not args.skip_downstream_optimization,
        args.downstream_perturbation, args.xtb,
    )
    version_output = subprocess.run([args.xtb, "--version"], capture_output=True,
                                    text=True, check=False).stdout
    xtb_version = next((line.strip() for line in version_output.splitlines()
                        if "xtb version" in line.lower()), "unparsed; inspect xtb --version")
    revision = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True,
                              text=True, check=False).stdout.strip() or None
    dirty = bool(subprocess.run(["git", "status", "--porcelain"], capture_output=True,
                                text=True, check=False).stdout.strip())
    comment_energy_deltas = [r["source_energy"] - r["xtb_energy_hartree"] for r in optimized]
    manifest = {
        "ensemble": args.ensemble.name,
        "ensemble_sha256": _sha256(args.ensemble),
        "frames": len(optimized),
        "source_composition": optimized[0]["composition"],
        "source_xyz_comment_energy_minus_gfn2_final_energy": {
            "interpretation": "comparison only; source comment units/meaning are ambiguous",
            "mean_absolute_delta_numeric": float(np.mean(np.abs(comment_energy_deltas))),
            "maximum_absolute_delta_numeric": float(np.max(np.abs(comment_energy_deltas))),
        },
        "reference_basins": len(representatives),
        "reference_labels": label_path.name,
        "analysis": "results/comparison.json",
        "xTB": args.xtb,
        "xTB_version": xtb_version,
        "pyar_git_revision": revision,
        "pyar_worktree_dirty": dirty,
        "args": vars(args) | {"ensemble": args.ensemble.name, "output": args.output.name},
        "report_conditions": len(report["conditions"]),
    }
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"wrote {label_path} and {args.output / 'results/comparison.json'}")


if __name__ == "__main__":
    main()
