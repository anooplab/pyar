"""Cross-domain analysis of conformer, atomic-cluster, and water-cluster pilots.

Hard classification metrics are computed only where independent operational
basin labels exist (the conformer pilots). Water metrics use explicit RMSD
proxy cutoffs and are reported separately. Atomic-cluster pilots remain
coverage/energy analyses because they have no validated similarity labels.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from sklearn.metrics import roc_auc_score

from pyar.scripts.scientific_clustering_benchmark import (
    ALGORITHMS,
    FEATURES,
    _as_molecule,
    load_ensemble,
)
from pyar.selection.distances import pairwise_distances
from pyar.selection.features import compute_feature_matrix, standardize_features


CONFORMERS = ("CAMVES_I", "FGG55", "WG01")
WATER_REFERENCE_CUTOFFS = (0.60, 0.75, 1.00)
SPECIFICITY_TARGETS = (0.95, 0.99)


def _current_algorithm_name(name):
    """Map the historical ``hybrid`` spelling onto today's ``auto`` policy."""
    return "auto" if str(name).strip().lower() == "hybrid" else name


def _safe_divide(numerator, denominator):
    return float(numerator / denominator) if denominator else None


def _binary_metrics(tp, fp, fn, tn):
    sensitivity = _safe_divide(tp, tp + fn)
    specificity = _safe_divide(tn, tn + fp)
    precision = _safe_divide(tp, tp + fp)
    return {
        "tp": int(tp), "fp": int(fp), "fn": int(fn), "tn": int(tn),
        "sensitivity": sensitivity,
        "specificity": specificity,
        "precision": precision,
        "false_positive_rate": _safe_divide(fp, fp + tn),
        "false_discovery_rate": _safe_divide(fp, fp + tp),
        "miss_rate": _safe_divide(fn, fn + tp),
        "f1": _safe_divide(2 * tp, 2 * tp + fp + fn),
    }


def _same_basin_labels(reference_path: Path, input_count: int) -> np.ndarray:
    reference = json.loads(reference_path.read_text())
    basin_by_name = {row["name"]: row["basin_id"] for row in reference["structures"]}
    names = [f"frame_{index:04d}" for index in range(input_count)]
    if (len(reference["structures"]) != len(basin_by_name)
            or set(basin_by_name) != set(names)):
        raise ValueError(f"Reference-basin labels do not cover all input frames in {reference_path}")
    return np.asarray([basin_by_name[name] for name in names], dtype=int)


def _feature_roc_auc(molecules, basin_ids, feature, system_type="conformers"):
    result = compute_feature_matrix(
        molecules, feature, allow_fallbacks=False, system_type=system_type,
        algorithm="agglomerative",
    )
    distances = pairwise_distances(standardize_features(result.values), metric="euclidean")
    indices = np.triu_indices(len(molecules), k=1)
    truth = basin_ids[indices[0]] == basin_ids[indices[1]]
    scores = -distances[indices]
    auc = float(roc_auc_score(truth, scores)) if len(set(truth.tolist())) == 2 else None
    return auc, result.name, list(result.fallbacks), distances


def analyze_conformers(run_root: Path, data_root: Path):
    confusion_rows, auc_rows, condition_rows = [], [], []
    for dataset in CONFORMERS:
        source_path = data_root / "mpconf196gen" / f"{dataset}_crest_conformers.xyz"
        manifest = json.loads((run_root / dataset / "manifest.json").read_text())
        if manifest.get("ensemble_sha256") != hashlib.sha256(source_path.read_bytes()).hexdigest():
            raise ValueError(f"Source ensemble hash mismatch for {dataset}; labels cannot be paired")
        ensemble = load_ensemble(source_path)
        molecules = [_as_molecule(row["atoms"], row["name"]) for row in ensemble]
        report_path = run_root / dataset / "results" / "comparison.json"
        reference_path = run_root / dataset / "reference_basins.json"
        report = json.loads(report_path.read_text())
        basin_ids = _same_basin_labels(reference_path, len(molecules))
        truth_counts = np.unique(basin_ids, return_counts=True)[1]
        positive_pairs = int(sum(count * (count - 1) // 2 for count in truth_counts))
        total_pairs = len(molecules) * (len(molecules) - 1) // 2

        for feature in FEATURES:
            try:
                auc, feature_used, fallbacks, _ = _feature_roc_auc(
                    molecules, basin_ids, feature,
                )
                auc_rows.append({
                    "dataset": dataset, "feature_requested": feature,
                    "feature_used": feature_used, "feature_fallbacks": json.dumps(fallbacks),
                    "status": "available", "error": "",
                    "pairwise_same_basin_auc": auc, "positive_same_basin_pairs": positive_pairs,
                    "negative_different_basin_pairs": total_pairs - positive_pairs,
                    "pair_count": total_pairs,
                })
            except Exception as exc:
                auc_rows.append({
                    "dataset": dataset, "feature_requested": feature,
                    "feature_used": None, "feature_fallbacks": "[]",
                    "status": "unavailable", "error": f"{type(exc).__name__}: {exc}",
                    "pairwise_same_basin_auc": None, "positive_same_basin_pairs": positive_pairs,
                    "negative_different_basin_pairs": total_pairs - positive_pairs,
                    "pair_count": total_pairs,
                })

        for condition in report["conditions"]:
            labels = np.asarray(condition["diagnostics"]["labels"], dtype=int)
            if labels.shape != basin_ids.shape:
                raise ValueError(f"Cluster and reference labels are misaligned for {dataset}")
            tp = fp = fn = tn = 0
            for left in range(len(labels)):
                for right in range(left + 1, len(labels)):
                    same_basin = basin_ids[left] == basin_ids[right]
                    # Noise records are deliberately treated as separate retained candidates.
                    same_cluster = (labels[left] >= 0 and labels[right] >= 0
                                    and labels[left] == labels[right])
                    if same_basin and same_cluster:
                        tp += 1
                    elif not same_basin and same_cluster:
                        fp += 1
                    elif same_basin:
                        fn += 1
                    else:
                        tn += 1
            metrics = _binary_metrics(tp, fp, fn, tn)
            unique_basin_capacity = min(report["seed_budget"], report["reference_basin_count"])
            budget_efficiency = _safe_divide(condition["basins_retained"], unique_basin_capacity)
            row = {
                "dataset": dataset,
                "algorithm_requested": _current_algorithm_name(condition["algorithm_requested"]),
                "algorithm_recorded": condition["algorithm_requested"],
                "algorithm_used": condition["diagnostics"].get("algorithm_used"),
                "feature_requested": condition["feature_requested"],
                "feature_used": condition["diagnostics"].get("feature_used"),
                "input_count": len(labels), "reference_basin_count": len(set(basin_ids)),
                "same_basin_pairs": positive_pairs,
                "different_basin_pairs": total_pairs - positive_pairs,
                "basins_retained": condition["basins_retained"],
                "maximum_distinct_basins_with_seed_budget": unique_basin_capacity,
                "seed_budget_basin_efficiency": budget_efficiency,
                **metrics,
            }
            confusion_rows.append(row)
            condition_rows.append(row)

    macro_rows = []
    for algorithm in ALGORITHMS:
        for feature in FEATURES:
            rows = [row for row in condition_rows
                    if row["algorithm_requested"] == algorithm
                    and row["feature_requested"] == feature]
            if not rows:
                continue
            macro_rows.append({
                "algorithm_requested": algorithm,
                "feature_requested": feature,
                "macro_sensitivity": _mean_available(row["sensitivity"] for row in rows),
                "macro_specificity": _mean_available(row["specificity"] for row in rows),
                "macro_precision": _mean_available(row["precision"] for row in rows),
                "macro_false_positive_rate": _mean_available(row["false_positive_rate"] for row in rows),
                "macro_false_discovery_rate": _mean_available(row["false_discovery_rate"] for row in rows),
                "macro_miss_rate": _mean_available(row["miss_rate"] for row in rows),
                "macro_f1": _mean_available(row["f1"] for row in rows),
                "macro_seed_budget_basin_efficiency": _mean_available(
                    row["seed_budget_basin_efficiency"] for row in rows
                ),
                "per_dataset": {
                    row["dataset"]: {
                        key: row[key] for key in (
                            "sensitivity", "specificity", "precision", "false_positive_rate",
                            "false_discovery_rate",
                            "miss_rate", "f1", "seed_budget_basin_efficiency",
                            "tp", "fp", "fn", "tn",
                        )
                    } for row in rows
                },
            })

    auc_macro_rows = []
    for feature in FEATURES:
        rows = [row for row in auc_rows if row["feature_requested"] == feature]
        auc_macro_rows.append({
            "feature_requested": feature,
            "macro_pairwise_same_basin_auc": _mean_available(
                row["pairwise_same_basin_auc"] for row in rows
            ),
            "per_dataset_auc": {row["dataset"]: row["pairwise_same_basin_auc"] for row in rows},
            "per_dataset_pair_counts": {row["dataset"]: row["pair_count"] for row in rows},
        })
    return confusion_rows, macro_rows, auc_rows, auc_macro_rows


def _mean_available(values):
    numbers = [float(value) for value in values if value is not None and np.isfinite(value)]
    return float(np.mean(numbers)) if numbers else None


def _conservative_operating_point(reference, candidate, target_specificity):
    reference = np.asarray(reference, dtype=float)
    candidate = np.asarray(candidate, dtype=float)
    truth = reference <= 0.75
    choices = []
    for threshold in [None, *np.unique(candidate).tolist()]:
        predicted = np.zeros_like(truth) if threshold is None else candidate <= threshold
        tp = int(np.count_nonzero(truth & predicted))
        fp = int(np.count_nonzero(~truth & predicted))
        fn = int(np.count_nonzero(truth & ~predicted))
        tn = int(np.count_nonzero(~truth & ~predicted))
        metrics = _binary_metrics(tp, fp, fn, tn)
        specificity = metrics["specificity"]
        if specificity is not None and specificity >= target_specificity:
            # Maximize sensitivity; among ties, keep the more conservative threshold.
            choices.append((metrics["sensitivity"], specificity,
                            -(threshold if threshold is not None else -np.inf),
                            threshold, metrics))
    if not choices:
        return None
    _, _, _, threshold, metrics = max(choices)
    return {"feature_threshold": threshold, **metrics}


def analyze_water(runs_root: Path):
    from sklearn.metrics import roc_auc_score

    path = runs_root / "water_similarity"
    manifest = json.loads((path / "manifest.json").read_text())
    reference = np.loadtxt(path / "reference_fragment_rmsd_upper_bound_angstrom.csv", delimiter=",")
    pair_index = np.triu_indices(len(reference), k=1)
    reference_values = reference[pair_index]
    rows = []
    for cutoff in WATER_REFERENCE_CUTOFFS:
        truth = reference_values <= cutoff
        for feature in FEATURES:
            matrix = np.loadtxt(path / f"{feature}_euclidean_feature_distances.csv", delimiter=",")
            candidate = matrix[pair_index]
            auc = (float(roc_auc_score(truth, -candidate))
                   if len(set(truth.tolist())) == 2 else None)
            row = {
                "reference_rmsd_proxy_cutoff_angstrom": cutoff,
                "feature": feature,
                "reference_near_pairs": int(truth.sum()),
                "reference_far_proxy_pairs": int((~truth).sum()),
                "proxy_pair_auc": auc,
                "specificity_target": None,
                "feature_distance_threshold": None,
                # This row reports ranking only. Zero confusion counts would
                # falsely imply that an all-negative classifier was evaluated.
                "tp": None, "fp": None, "fn": None, "tn": None,
                "sensitivity": None, "specificity": None, "precision": None,
                "false_positive_rate": None, "false_discovery_rate": None,
                "miss_rate": None, "f1": None,
                "reference_method": manifest["reference_similarity"]["method"],
            }
            rows.append(row)
            if cutoff == 0.75:
                for target in SPECIFICITY_TARGETS:
                    operating = _conservative_operating_point(
                        reference_values, candidate, target,
                    )
                    if operating is not None:
                        rows.append({
                            "reference_rmsd_proxy_cutoff_angstrom": cutoff,
                            "feature": feature,
                            "reference_near_pairs": int(truth.sum()),
                            "reference_far_proxy_pairs": int((~truth).sum()),
                            "proxy_pair_auc": auc,
                            "specificity_target": target,
                            "feature_distance_threshold": operating.pop("feature_threshold"),
                            **operating,
                            "reference_method": manifest["reference_similarity"]["method"],
                        })
    return rows


def analyze_atomic(runs_root: Path):
    rows = []
    for dataset in ("lj13", "qcd_au13"):
        report = json.loads((runs_root / "atomic_clusters" / dataset / "comparison.json").read_text())
        for condition in report["conditions"]:
            rows.append({
                "dataset": dataset,
                "input_count": report["input_count"],
                "algorithm_requested": _current_algorithm_name(condition["algorithm_requested"]),
                "algorithm_recorded": condition["algorithm_requested"],
                "algorithm_used": condition["algorithm_used"],
                "feature_requested": condition["feature_requested"],
                "feature_used": condition["feature_used"],
                "selected_count": condition["selected_count"],
                "lowest_source_energy_entry_retained": condition["lowest_source_energy_entry_retained"],
                "top_five_source_energy_entries_retained": condition["top_five_source_energy_entries_retained"],
                "normalized_mean_pair_spectrum_rms": condition["normalized_mean_pair_spectrum_rms"],
                "runtime_warning_count": len(condition.get("runtime_warnings", [])),
                "reference_label_status": "not available; source entries are not validated similarity classes",
            })
    return rows


def _write_csv(path: Path, rows):
    if not rows:
        path.write_text("")
        return
    columns = list(dict.fromkeys(key for row in rows for key in row if not isinstance(row[key], dict)))
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        writer.writerows({key: value for key, value in row.items() if key in columns} for row in rows)


def build_meta_analysis(repository_root: Path, output_dir: Path):
    benchmark_root = repository_root / "benchmarks" / "clustering_scientific"
    runs_root = benchmark_root / "runs"
    data_root = benchmark_root / "data"
    output_dir.mkdir(parents=True, exist_ok=True)

    confusion, conformer_macro, auc_rows, auc_macro = analyze_conformers(runs_root, data_root)
    water_rows = analyze_water(runs_root)
    atomic_rows = analyze_atomic(runs_root)

    _write_csv(output_dir / "conformer_pair_confusion.csv", confusion)
    _write_csv(output_dir / "conformer_condition_macro.csv", conformer_macro)
    _write_csv(output_dir / "conformer_feature_roc_auc.csv", auc_rows)
    _write_csv(output_dir / "water_proxy_roc_auc_confusion.csv", water_rows)
    _write_csv(output_dir / "atomic_cluster_proxy_metrics.csv", atomic_rows)
    summary = {
        "schema_version": 2,
        "scope": ["molecular conformers", "LJ13 atomic clusters", "QCD Au13 clusters", "W6 water clusters"],
        "conformer_pair_confusion": confusion,
        "conformer_condition_macro": conformer_macro,
        "conformer_feature_roc_auc": auc_rows,
        "conformer_feature_roc_auc_macro": auc_macro,
        "water_proxy_roc_auc_confusion": water_rows,
        "atomic_cluster_proxy_metrics": atomic_rows,
        "pooling_policy": "Metrics are not pooled across datasets with different label definitions or without validated labels. Conformer AUCs are macro-averaged across three molecular systems; water uses RMSD-proxy cutoffs and atomic clusters report coverage only.",
    }
    (output_dir / "meta_analysis.json").write_text(json.dumps(summary, indent=2) + "\n")
    report = _render_report(summary)
    (output_dir / "report.md").write_text(report)
    print(f"Wrote cross-domain clustering analysis to {output_dir}")
    return summary


def _fmt(value, digits=3):
    return "—" if value is None else f"{value:.{digits}f}"


def _render_report(summary):
    feature_rows = [row for row in summary["conformer_feature_roc_auc_macro"]
                    if row["macro_pairwise_same_basin_auc"] is not None]
    best_feature = max(feature_rows, key=lambda row: row["macro_pairwise_same_basin_auc"])["feature_requested"] if feature_rows else None
    lines = [
        "# Cross-domain clustering meta-analysis",
        "",
        "This analysis compares the molecular conformer, atomic-cluster, and water-cluster pilots. It keeps each domain's reference label and outcome separate; a pooled AUC across these unlike targets would be misleading.",
        "",
        "## Reference evidence by domain",
        "",
        "| Domain | Available reference | What can be measured | Limitation |",
        "|---|---|---|---|",
        "| Molecular conformers | Independent GFN2-xTB optimizations; graph identity plus heavy-atom Kabsch RMSD below 0.35 Å | Basin-pair ROC/AUC and confusion matrices for each clusterer-feature setting | Three macrocycles; pair counts are dependent and class balance differs |",
        "| LJ13 and Au13 atomic clusters | Distinct source minima/entries, energies, and pair-spectrum coverage | Source-energy retention and normalized geometric coverage | No independently validated same/different similarity labels; no sensitivity/specificity/AUC |",
        "| W6 water clusters | All-atom fragment-matched Kabsch RMSD upper-bound distances | Feature ROC/AUC and confusion matrices at declared distance cutoffs | One composition; proxy labels, not verified basin labels; reference distances are upper bounds |",
        "",
        "## Conformer feature ROC/AUC",
        "",
        "The positive class is a pair assigned to the same operational reference basin. Lower feature distances predict a positive pair. Values are ROC AUC by source pool and the unweighted macro-average.",
        "",
        "| Feature | CAMVES_I | FGG55 | WG01 | Macro AUC |",
        "|---|---:|---:|---:|---:|",
    ]
    for row in summary["conformer_feature_roc_auc_macro"]:
        per = row["per_dataset_auc"]
        lines.append(
            f"| {row['feature_requested']} | {_fmt(per['CAMVES_I'])} | {_fmt(per['FGG55'])} | {_fmt(per['WG01'])} | {_fmt(row['macro_pairwise_same_basin_auc'])} |"
        )
    lines += [
        "",
        (f"{best_feature} has the highest observed same-basin pair AUC in this three-pool comparison. "
         "AUC ranks pair distances; it does not by itself establish a safe merge threshold."
         if best_feature else "No feature was available across the conformer pools for an AUC comparison."),
        "",
        "## Cluster-assignment confusion summary",
        "",
        "For each clusterer-feature condition, a pair is predicted positive if it has the same non-noise cluster label. Noise points are treated as separate retained candidates. The table below macro-averages metrics across the three conformer pools; undefined metrics in a pool are omitted from that metric's mean.",
        "",
        "| Algorithm | Feature | Sensitivity | Specificity | Precision | False-positive rate | False-discovery rate | Miss rate | Basin budget efficiency |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in sorted(
        summary["conformer_condition_macro"],
        key=lambda item: (-(item["macro_specificity"] or 0), -(item["macro_sensitivity"] or 0)),
    ):
        lines.append(
            f"| {row['algorithm_requested']} | {row['feature_requested']} | {_fmt(row['macro_sensitivity'])} | {_fmt(row['macro_specificity'])} | {_fmt(row['macro_precision'])} | {_fmt(row['macro_false_positive_rate'])} | {_fmt(row['macro_false_discovery_rate'])} | {_fmt(row['macro_miss_rate'])} | {_fmt(row['macro_seed_budget_basin_efficiency'])} |"
        )
    lines += [
        "",
        "The class imbalance matters: FGG55 and WG01 contain many more different-basin pairs than same-basin pairs, so specificity can look high while precision remains low. Use the per-pool 2×2 counts in `conformer_pair_confusion.csv` alongside these rates. The seed-budget basin efficiency divides recovered basins by the maximum distinct basins the fixed 12-seed budget could retain; it avoids calling 12/98 a 12% failure when the budget itself is 12.",
        "",
        "## Water-cluster proxy ROC/AUC and conservative operating points",
        "",
        "Reference-positive pairs are those with fragment-RMSD upper bound at or below the listed cutoff. For the 0.75 Å proxy cutoff, feature thresholds are selected on this same 24-structure sample to maximize sensitivity subject to the listed minimum specificity. These are apparent, in-sample operating points, not validated deployment thresholds.",
        "",
        "| Feature | Reference cutoff (Å) | Positive pairs | AUC | Specificity target | Actual specificity | Sensitivity | TP | FP | FN | TN |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in summary["water_proxy_roc_auc_confusion"]:
        if row["specificity_target"] is None and row["reference_rmsd_proxy_cutoff_angstrom"] != 0.75:
            continue
        lines.append(
            f"| {row['feature']} | {row['reference_rmsd_proxy_cutoff_angstrom']:.2f} | {row['reference_near_pairs']} | {_fmt(row['proxy_pair_auc'])} | {_fmt(row['specificity_target'], 2)} | {_fmt(row['specificity'])} | {_fmt(row['sensitivity'])} | {_fmt(row['tp'], 0)} | {_fmt(row['fp'], 0)} | {_fmt(row['fn'], 0)} | {_fmt(row['tn'], 0)} |"
        )
    lines += [
        "",
        "At the 0.75 Å proxy cutoff there are 15 positive pairs and 261 proxy-negative pairs. The high-specificity points deliberately accept false negatives: with a keep-in-doubt rule, an uncertain pair is retained rather than merged. Because the reference distance is an upper bound, proxy-negative pairs are not proven dissimilar; the resulting specificity is against this operational proxy only.",
        "",
        "## Atomic-cluster outcomes",
        "",
        "All 15 tested conditions on each atomic dataset retained its lowest-source-energy entry. Best normalized pair-spectrum coverage differs by dataset:",
        "",
        "| Dataset | Best setting by pair-spectrum coverage | Selected | Normalized mean nearest-spectrum RMS |",
        "|---|---|---:|---:|",
    ]
    for dataset in ("lj13", "qcd_au13"):
        rows = [row for row in summary["atomic_cluster_proxy_metrics"] if row["dataset"] == dataset]
        best = min(rows, key=lambda row: row["normalized_mean_pair_spectrum_rms"])
        lines.append(
            f"| {dataset} | {best['algorithm_requested']} / {best['feature_requested']} | {best['selected_count']} | {best['normalized_mean_pair_spectrum_rms']:.5f} |"
        )
    lines += [
        "",
        "Sorted pair-distance spectra are not unique structural identifiers, so these coverage results cannot be converted into confusion matrices without new independent labels.",
        "",
        "## Provisional policy",
        "",
        "1. **Deduplication:** retain the conservative `in doubt, keep` rule. Use feature distances only to find candidate neighbors; merge only after a complete, chemistry-aware structural comparison confirms equivalence. Incomplete graph mappings or threshold-borderline matches remain separate.",
        "2. **Conformer feature:** SOAP is the strongest tested basin-pair ranking feature across the three macrocycle pools. Use it as the leading candidate for further held-out testing, not as a universal default yet.",
        "3. **Clustering algorithm:** no algorithm-feature pair dominates the cluster confusion, specificity, and seed-budget metrics. The HDBSCAN-first `auto` policy is not independently validated against DBSCAN or OPTICS by this analysis; `maxmin` shares the automatic/HDBSCAN cluster labels and is a cluster-then-trim selection mode, not a separate partitioner.",
        "4. **Atomic and aggregate systems:** keep feature/algorithm choices system-aware. Two atomic datasets disagree on their best proxy setting, and the W6 benchmark shows descriptor ranking depends on whether the objective is global rank correlation or retrieval of the nearest pairs.",
        "5. **Operating point:** the high-specificity water settings miss most proxy-positive pairs in this sample. Do not adopt these in-sample thresholds; evaluate retrieval recall and false-discovery rate on held-out cluster sizes and a chemically diverse aggregate set before choosing an operating point.",
        "6. **Cutoff interpretation:** the 0.35 Å conformer reference label is not the production duplicate-removal cutoff. Current generation and backend-final reduction floors are 0.50 and 0.75 Å, respectively; this benchmark does not establish those broader removals as safe. Keep the distinction explicit until direct deletion-safety labels are available.",
        "",
        "## Reproducibility and caution",
        "",
        "The pairwise entries are not independent observations, and this is a small exploratory set. Conformer reference basins are operational labels produced by the saved graph/Kabsch protocol. Water thresholds are RMSD-proxy categories. Atomic-cluster similarity labels remain unavailable. No single pooled accuracy number is reported across these different targets.",
        "",
        "Rebuild this analysis with `python -m pyar.scripts.meta_analyze_clustering_benchmarks --output benchmarks/clustering_scientific/meta_analysis`.",
    ]
    return "\n".join(lines) + "\n"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository-root", type=Path, default=Path.cwd())
    parser.add_argument("--output", type=Path,
                        default=Path("benchmarks/clustering_scientific/meta_analysis"))
    args = parser.parse_args(argv)
    output_dir = args.output if args.output.is_absolute() else args.repository_root / args.output
    build_meta_analysis(args.repository_root, output_dir)


if __name__ == "__main__":
    main()
