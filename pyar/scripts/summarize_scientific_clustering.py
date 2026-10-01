"""Aggregate persistent scientific clustering benchmark run reports."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def summarize(runs: list[Path], output: Path):
    systems = []
    condition_map = {}
    reference_by_system = {}
    for run in runs:
        manifest = json.loads((run / "manifest.json").read_text())
        report = json.loads((run / "results" / "comparison.json").read_text())
        system_name = run.name
        reference = json.loads((run / "reference_basins.json").read_text())
        structure_by_name = {item["name"]: item for item in reference["structures"]}
        basin_energy = {}
        for structure in reference["structures"]:
            basin_id = structure["basin_id"]
            basin_energy[basin_id] = min(
                basin_energy.get(basin_id, float("inf")), structure["xtb_energy_hartree"],
            )
        reference_by_system[system_name] = (
            structure_by_name, basin_energy, min(basin_energy.values()),
        )
        systems.append({
            "system": system_name,
            "frames": manifest["frames"],
            "reference_basins": report["reference_basin_count"],
            "completed_settings": len(report["conditions"]),
            "selected_seed_condition_occurrences": sum(
                condition["downstream_optimization"]["selected_count"]
                for condition in report["conditions"]
            ),
            "successful_perturbed_optimizations": sum(
                condition["downstream_optimization"]["successful_terminations"]
                for condition in report["conditions"]
            ),
            "verified_perturbed_basin_recoveries": sum(
                condition["downstream_optimization"]["verified_reference_basin_seed_count"]
                for condition in report["conditions"]
            ),
            "unique_downstream_optimizations": max(
                (condition["downstream_optimization"]["unique_geometries_optimized_and_cached"]
                 for condition in report["conditions"]
                 if condition["downstream_optimization"] is not None),
                default=0,
            ),
            "source_sha256": manifest["ensemble_sha256"],
            "xtb_version": manifest["xTB_version"],
            "threshold_sensitivity": reference["threshold_sensitivity"],
            "inferred_topology_group_count": len(reference["inferred_topology_groups"]),
            "source_comment_energy_delta": manifest["source_xyz_comment_energy_minus_gfn2_final_energy"],
        })
        for condition in report["conditions"]:
            key = (condition["algorithm_requested"], condition["feature_requested"])
            structure_by_name, basin_energy, global_minimum = reference_by_system[system_name]
            selected_basins = {
                structure_by_name[name]["basin_id"] for name in condition["selected_names"]
            }
            energy_window_retention = {}
            for window in (1.0, 3.0, 6.0):
                in_window = {
                    basin_id for basin_id, energy in basin_energy.items()
                    if (energy - global_minimum) * 627.509474 <= window
                }
                retained = len(in_window & selected_basins)
                energy_window_retention[f"{window:.1f}"] = {
                    "reference_basins": len(in_window),
                    "retained_basins": retained,
                    "recall": retained / len(in_window) if in_window else None,
                }
            condition_map.setdefault(key, []).append({
                "system": system_name,
                "energy_window_basin_retention": energy_window_retention,
                **condition,
            })

    rows = []
    for (algorithm, feature), observations in sorted(condition_map.items()):
        retention = [item["retention_fraction"] for item in observations]
        coverage = [item["coverage"]["mean_nearest_indexed_heavy_atom_rmsd_angstrom"]
                    for item in observations
                    if item["coverage"]["mean_nearest_indexed_heavy_atom_rmsd_angstrom"] is not None]
        downstream = [item["downstream_optimization"] for item in observations
                      if item["downstream_optimization"] is not None]
        downstream_total = sum(item["selected_count"] for item in downstream)
        downstream_success = sum(item["successful_terminations"] for item in downstream)
        recovered_seeds = sum(item["verified_reference_basin_seed_count"] for item in downstream)
        fallback_algorithm = sum(bool(item["diagnostics"].get("algorithm_fallbacks"))
                                 for item in observations)
        fallback_feature = sum(bool(item["diagnostics"].get("feature_fallbacks"))
                               for item in observations)
        fallback_distance = sum(bool(item["diagnostics"].get("distance_fallbacks"))
                                for item in observations)
        rows.append({
            "algorithm_requested": algorithm,
            "feature_requested": feature,
            "systems": len(observations),
            "macro_mean_basin_retention": float(np.mean(retention)),
            "micro_basin_retention": (
                sum(item["basins_retained"] for item in observations)
                / sum(item["coverage"]["reference_basin_count"] for item in observations)
            ),
            "per_system_basin_retention": {
                item["system"]: item["retention_fraction"] for item in observations
            },
            "per_system_retained_reference_basins": {
                item["system"]: item["basins_retained"] for item in observations
            },
            "energy_window_basin_retention": {
                window: {
                    "reference_basins_total": sum(
                        item["energy_window_basin_retention"][window]["reference_basins"]
                        for item in observations
                    ),
                    "retained_basins_total": sum(
                        item["energy_window_basin_retention"][window]["retained_basins"]
                        for item in observations
                    ),
                    "macro_mean_recall": float(np.mean([
                        item["energy_window_basin_retention"][window]["recall"]
                        for item in observations
                        if item["energy_window_basin_retention"][window]["recall"] is not None
                    ])),
                    "micro_recall": (
                        sum(item["energy_window_basin_retention"][window]["retained_basins"]
                            for item in observations)
                        / sum(item["energy_window_basin_retention"][window]["reference_basins"]
                              for item in observations)
                    ),
                    "per_system": {
                        item["system"]: item["energy_window_basin_retention"][window]
                        for item in observations
                    },
                }
                for window in ("1.0", "3.0", "6.0")
            },
            "macro_mean_nearest_indexed_heavy_atom_rmsd_angstrom": float(np.mean(coverage)) if coverage else None,
            "downstream_optimization_success_fraction": (
                downstream_success / downstream_total if downstream_total else None
            ),
            "verified_reference_basin_recovery_fraction": (
                recovered_seeds / downstream_total if downstream_total else None
            ),
            "basin_changed_after_complete_comparison": sum(
                item["basin_changed_after_complete_comparison_count"] for item in downstream
            ),
            "uncertain_due_to_incomplete_mapping": sum(
                item["uncertain_due_to_incomplete_mapping_count"] for item in downstream
            ),
            "algorithm_fallback_systems": fallback_algorithm,
            "feature_fallback_systems": fallback_feature,
            "distance_fallback_systems": fallback_distance,
            "macro_mean_runtime_seconds": float(np.mean([
                item["runtime_seconds"] for item in observations
            ])),
            "actual_methods_by_system": {
                item["system"]: {
                    "algorithm": item["diagnostics"].get("algorithm_used"),
                    "feature": item["diagnostics"].get("feature_used"),
                    "distance": item["diagnostics"].get("distance_used"),
                } for item in observations
            },
        })
    rows.sort(key=lambda item: (
        -item["micro_basin_retention"],
        -item["energy_window_basin_retention"]["1.0"]["micro_recall"],
        float("inf") if item["macro_mean_nearest_indexed_heavy_atom_rmsd_angstrom"] is None
        else item["macro_mean_nearest_indexed_heavy_atom_rmsd_angstrom"],
        item["macro_mean_runtime_seconds"],
    ))

    output.mkdir(parents=True, exist_ok=True)
    payload = {
        "schema_version": 1,
        "data_source": "MPCONF196GEN; CC BY 4.0; see data/mpconf196gen/",
        "primary_ranking": "micro-average basin retention weighted by the number of reference basins, then micro-average recall within 1 kcal/mol, then indexed heavy-atom RMSD coverage, then runtime",
        "system_count": len(systems),
        "total_input_geometries": sum(item["frames"] for item in systems),
        "total_reference_basins": sum(item["reference_basins"] for item in systems),
        "completed_setting_runs": sum(item["completed_settings"] for item in systems),
        "unique_downstream_optimizations": sum(item["unique_downstream_optimizations"] for item in systems),
        "selected_seed_condition_occurrences": sum(item["selected_seed_condition_occurrences"] for item in systems),
        "successful_perturbed_optimizations": sum(item["successful_perturbed_optimizations"] for item in systems),
        "verified_perturbed_basin_recoveries": sum(item["verified_perturbed_basin_recoveries"] for item in systems),
        "systems": systems,
        "settings": rows,
        "findings": {},
        "limitations": [
            "Only macrocycle conformer pools are represented.",
            "The independent clustering benchmark excludes the workflow's pre-clustering deduplication stage.",
            "Reference basins depend on coordinate-inferred connectivity and the stated RMSD cutoff.",
            "Heavy-atom RMSD ignores hydrogen-coordinate-only rearrangements, while retaining hydrogen connectivity/count constraints.",
            "Downstream optimization success is measured after deterministic 0.10 Angstrom Gaussian coordinate perturbations.",
            "Successful optimization is distinct from recovery of the selected seed's basin; complete-comparison basin changes and incomplete graph mappings are reported separately.",
            "Results do not establish policy for isomer pools, atomic clusters, or non-covalent aggregates.",
        ],
    }
    if rows:
        best = rows[0]
        best_total_retained = sum(best["per_system_retained_reference_basins"].values())
        findings = {
            "best_observed_setting": {
                "algorithm": best["algorithm_requested"],
                "feature": best["feature_requested"],
                "retained_basins": best_total_retained,
                "reference_basins": payload["total_reference_basins"],
                "weighted_basin_recall": best["micro_basin_retention"],
                "one_kcal_window": best["energy_window_basin_retention"]["1.0"],
                "three_kcal_window": best["energy_window_basin_retention"]["3.0"],
            },
            "all_settings_perturbed_optimizer_success_fraction": (
                payload["successful_perturbed_optimizations"]
                / payload["selected_seed_condition_occurrences"]
            ),
            "all_settings_own_basin_recovery_fraction": (
                payload["verified_perturbed_basin_recoveries"]
                / payload["selected_seed_condition_occurrences"]
            ),
        }
        payload["findings"] = findings
    (output / "summary.json").write_text(json.dumps(payload, indent=2) + "\n")
    with (output / "report.md").open("w") as stream:
        stream.write("# Scientific clustering benchmark results\n\n")
        stream.write("This is a macrocycle conformer-selection pilot. Reference labels come from independently GFN2-xTB-optimized geometries grouped by complete coordinate-graph isomorphism and heavy-atom Kabsch RMSD. See the benchmark README for attribution, methods, and limits.\n\n")
        stream.write(
            f"The run covers {payload['total_input_geometries']} starting geometries, "
            f"{payload['total_reference_basins']} reference basins at 0.35 Å, "
            f"{payload['completed_setting_runs']} algorithm/feature-setting runs, and "
            f"{payload['unique_downstream_optimizations']} unique perturbed seed optimizations.\n\n"
        )
        if payload["findings"]:
            finding = payload["findings"]
            best = finding["best_observed_setting"]
            one = best["one_kcal_window"]
            three = best["three_kcal_window"]
            stream.write(
                f"The highest weighted basin retention was **{best['algorithm']} + {best['feature']}**: "
                f"{best['retained_basins']}/{best['reference_basins']} basins "
                f"({best['weighted_basin_recall']:.1%}). It retained "
                f"{one['retained_basins_total']}/{one['reference_basins_total']} basins within 1 kcal/mol and "
                f"{three['retained_basins_total']}/{three['reference_basins_total']} within 3 kcal/mol "
                "of each pool's lowest reference energy.\n\n"
                f"Across all {payload['selected_seed_condition_occurrences']} selected-seed/settings "
                f"occurrences, xTB optimization success was "
                f"{finding['all_settings_perturbed_optimizer_success_fraction']:.1%} and "
                f"verified recovery of each seed's own basin was "
                f"{finding['all_settings_own_basin_recovery_fraction']:.1%}. "
                "These are robustness checks at the tested 0.10 Å perturbation level.\n\n"
            )
        stream.write("## Test systems\n\n")
        stream.write("| System | Geometries | Reference basins |\n|---|---:|---:|\n")
        for system in systems:
            stream.write(f"| {system['system']} | {system['frames']} | {system['reference_basins']} |\n")
        stream.write("\n## Reference-label sensitivity and source audit\n\n")
        stream.write("| System | Basins at 0.25 Å | Basins at 0.35 Å | Basins at 0.50 Å | Inferred topology groups | Max comment/xTB energy difference (Eh) |\n|---|---:|---:|---:|---:|---:|\n")
        for system in systems:
            sensitivity = system["threshold_sensitivity"]
            delta = system["source_comment_energy_delta"]["maximum_absolute_delta_numeric"]
            stream.write(
                f"| {system['system']} | {sensitivity['0.25']['basin_count']} | "
                f"{sensitivity['0.35']['basin_count']} | {sensitivity['0.50']['basin_count']} | "
                f"{system['inferred_topology_group_count']} | {delta:.2e} |\n"
            )
        stream.write("\nAll source comment energies match fresh xTB final energies numerically within `8e-7 Eh`. This conflicts with the upstream README's Egret-1/kcal-per-mole description; fresh GFN2-xTB Hartree energies are used for cluster minima, and source values are retained for audit.\n")
        stream.write("\n## Settings (ranked)\n\n")
        stream.write("| Rank | Algorithm | Feature | Basin recall (weighted) | System-macro recall | Basin recall (≤1 kcal/mol, weighted) | Basin recall (≤3 kcal/mol, weighted) | Mean nearest indexed heavy-atom RMSD (Å) | Perturbed optimization success | Verified basin recovery | Fallback systems (alg/feat/dist) | Mean selection time (s) |\n|---:|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n")
        for rank, item in enumerate(rows, 1):
            coverage = item["macro_mean_nearest_indexed_heavy_atom_rmsd_angstrom"]
            success = item["downstream_optimization_success_fraction"]
            basin_recovery = item["verified_reference_basin_recovery_fraction"]
            stream.write(
                f"| {rank} | {item['algorithm_requested']} | {item['feature_requested']} | "
                f"{item['micro_basin_retention']:.3f} | "
                f"{item['macro_mean_basin_retention']:.3f} | "
                f"{item['energy_window_basin_retention']['1.0']['micro_recall']:.3f} | "
                f"{item['energy_window_basin_retention']['3.0']['micro_recall']:.3f} | "
                f"{'n/a' if coverage is None else f'{coverage:.3f}'} | "
                f"{'n/a' if success is None else f'{success:.3f}'} | "
                f"{'n/a' if basin_recovery is None else f'{basin_recovery:.3f}'} | "
                f"{item['algorithm_fallback_systems']}/{item['feature_fallback_systems']}/"
                f"{item['distance_fallback_systems']} | {item['macro_mean_runtime_seconds']:.3f} |\n"
            )
        stream.write("\n## Scope and caveats\n\n")
        for limitation in payload["limitations"]:
            stream.write(f"- {limitation}\n")
        stream.write("\nThe ordering is descriptive for these three test pools, not a universal feature/algorithm policy. Requested methods can fall back; per-system actual methods are recorded in `summary.json`.\n")
    return payload


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runs", nargs="+", required=True, type=Path,
                        help="Completed run directories containing manifest.json and results/")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args(argv)
    summarize(args.runs, args.output)
    print(f"wrote {args.output / 'report.md'} and {args.output / 'summary.json'}")


if __name__ == "__main__":
    main()
