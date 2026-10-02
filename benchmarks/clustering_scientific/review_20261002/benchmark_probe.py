"""Read benchmark evidence and demonstrate cache-key defects in isolated fixtures.

Run from repository root with .venv/bin/python <this file>.
No production or original benchmark artifacts are modified.
"""
import hashlib
import json
import sys
from pathlib import Path
from unittest.mock import patch

import numpy as np
from ase import Atoms
from ase.io import read, write

from pyar.scripts import scientific_clustering_benchmark as conformer
from pyar.scripts import water_cluster_similarity_benchmark as water
from pyar.scripts import meta_analyze_clustering_benchmarks as meta


ROOT = Path(__file__).resolve().parent
BENCH = ROOT.parent


def main():
    result = {}
    run = ROOT / "cache_fixture"
    frame = run / "optimizations" / "frame_0000"
    frame.mkdir(parents=True, exist_ok=True)
    old = Atoms("HH", positions=[[0, 0, 0], [0, 0, 0.74]])
    new = Atoms("HH", positions=[[0, 0, 0], [0, 0, 4.0]])
    write(frame / "xtbopt.xyz", old)
    (frame / "xtb.log").write_text(
        "TOTAL ENERGY -1.0\nGRADIENT NORM 0.0001\nnormal termination of xtb\n"
    )
    record = {"name": "frame_0000", "atoms": new}
    with patch.object(conformer.subprocess, "run", side_effect=AssertionError("reran optimizer")):
        reused = conformer.optimize_ensemble([record], run, sys.executable)[0]
    result["conformer_stale_cache_demonstration"] = {
        "new_input_bond_length": 4.0,
        "returned_cached_bond_length": float(reused["optimized_atoms"].get_distance(0, 1)),
        "optimizer_was_not_called": True,
    }

    source = BENCH / "data/water_clusters/W6_geoms_5.0_KCal-1hgztfv.txt"
    if not source.exists():
        source = next((BENCH / "data/water_clusters").rglob("W6_geoms_5.0_KCal-1hgztfv.txt"))
    fixture_source = ROOT / "water_source_fixture.txt"
    fixture_source.write_text(source.read_text())
    fixture_run = ROOT / "water_cache_fixture"
    fixture_run.mkdir(exist_ok=True)
    calls = []

    class FakeComparator:
        def __init__(self, **kwargs):
            pass

        def compare(self, first, second):
            calls.append(True)
            class Result:
                compatible = True
                distance = 1.234
                metadata = {"comparison_complete": True}
            return Result()

    with patch.object(water, "FragmentRMSDComparator", FakeComparator):
        first = water.build_benchmark(fixture_source, fixture_run, sample_size=2,
                                      features=(), resume=False)
        lines = fixture_source.read_text().splitlines()
        atoms = lines[2].split()
        atoms[1] = str(float(atoms[1]) + 0.01)
        lines[2] = " ".join(atoms)
        fixture_source.write_text("\n".join(lines) + "\n")
        second = water.build_benchmark(fixture_source, fixture_run, sample_size=2,
                                       features=(), resume=True)
    result["water_stale_cache_demonstration"] = {
        "source_hash_changed": first["source_sha256"] != second["source_sha256"],
        "total_comparator_calls_for_two_runs": len(calls),
        "expected_calls_if_inputs_checked": 2,
        "manifest_replaced_with_new_source_hash": True,
    }
    result["conformer_label_diagnostics"] = {}
    for name in ("CAMVES_I", "FGG55", "WG01"):
        reference = json.loads((BENCH / "runs" / name / "reference_basins.json").read_text())
        groups = {}
        for row in reference["structures"]:
            groups.setdefault(row["basin_id"], []).append(row["xtb_energy_hartree"])
        result["conformer_label_diagnostics"][name] = {
            "threshold_sensitivity": reference["threshold_sensitivity"],
            "maximum_same_label_energy_spread_kcal_mol": max(max(v)-min(v) for v in groups.values()) * 627.509474,
        }
        source = BENCH / "data/mpconf196gen" / f"{name}_crest_conformers.xyz"
        frames = read(source, index=":")
        deviation = max(float(np.abs(
            read(BENCH / "runs" / name / "optimizations" / f"frame_{i:04d}" / "input.xyz").positions
            - atoms.positions
        ).max()) for i, atoms in enumerate(frames))
        result["conformer_label_diagnostics"][name].update({
            "source_hash_matches_reference": hashlib.sha256(source.read_bytes()).hexdigest() == reference["source_sha256"],
            "maximum_current_source_vs_saved_optimizer_input_deviation": deviation,
        })
    saved = json.loads((BENCH / "meta_analysis" / "meta_analysis.json").read_text())
    water_rows = meta.analyze_water(BENCH / "runs")
    conformer_rows = meta.analyze_conformers(BENCH / "runs", BENCH / "data")
    result["independent_recomputation"] = {
        "water_metrics_match_saved": water_rows == saved["water_proxy_roc_auc_confusion"],
        "conformer_confusion_matches_saved": conformer_rows[0] == saved["conformer_pair_confusion"],
        "conformer_auc": conformer_rows[-1],
    }
    (ROOT / "benchmark_probe_results.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
