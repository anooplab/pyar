from pathlib import Path
import json
from types import SimpleNamespace
from unittest import mock

import numpy as np
from ase import Atoms

from pyar.scripts.scientific_clustering_benchmark import (
    _normally_terminated,
    assign_reference_basins,
    compare_perturbed_seeds,
    load_ensemble,
    optimize_ensemble,
)
from pyar.scripts.summarize_scientific_clustering import summarize


DATA = Path(__file__).parents[1] / "benchmarks" / "clustering_scientific" / "data" / "mpconf196gen"


def test_xTB_termination_parser_rejects_abnormal_phrase():
    assert _normally_terminated("normal termination of xtb\n")
    assert not _normally_terminated("abnormal termination of xtb\n")
    assert not _normally_terminated("xTB terminated without a final status\n")


def test_optimization_cache_is_invalidated_by_changed_geometry(tmp_path):
    xtb = tmp_path / "xtb"
    xtb.write_text("test executable identity")
    record = {"name": "frame_0000", "atoms": Atoms(["H", "H"], positions=[[0, 0, 0], [0, 0, .74]])}
    calls = []

    def fake_run(command, *, cwd, **kwargs):
        calls.append(command)
        from ase.io import write
        write(str(Path(cwd) / "xtbopt.xyz"), record["atoms"], format="xyz")
        return SimpleNamespace(returncode=0, stdout=(
            "TOTAL ENERGY -0.5\nTOTAL ENERGY -1.0\nGRADIENT NORM 0.0001\nnormal termination of xtb\n"
        ), stderr="")

    with mock.patch("pyar.scripts.scientific_clustering_benchmark.subprocess.run", side_effect=fake_run):
        result = optimize_ensemble([record], tmp_path, str(xtb))
        assert result[0]["xtb_energy_hartree"] == -1.0
        optimize_ensemble([record], tmp_path, str(xtb))
        record["atoms"].positions[1, 2] = .80
        optimize_ensemble([record], tmp_path, str(xtb))
        log = tmp_path / 'optimizations/frame_0000/xtb.log'
        log.write_text(log.read_text().replace('-1.0', '-2.0'))
        result = optimize_ensemble([record], tmp_path, str(xtb))
        assert result[0]['xtb_energy_hartree'] == -1.0

    assert len(calls) == 3


def test_curated_source_ensembles_have_stable_atom_order_and_energies():
    expected = {"CAMVES_I": 10, "FGG55": 127, "WG01": 72}
    for system, count in expected.items():
        frames = load_ensemble(DATA / f"{system}_crest_conformers.xyz")
        assert len(frames) == count
        assert len({tuple(frame["atoms"].get_chemical_symbols()) for frame in frames}) == 1
        assert np.isfinite([frame["source_energy"] for frame in frames]).all()


def test_basin_labels_merge_rigid_permutations_only_after_complete_graph_match():
    coords = np.array([
        [0.0, 0.0, 0.0],
        [0.629, 0.629, 0.629],
        [-0.629, -0.629, 0.629],
        [-0.629, 0.629, -0.629],
        [0.629, -0.629, -0.629],
    ])
    symbols = ["C", "H", "H", "H", "H"]
    rotated = coords @ np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]]) + [3.0, -2.0, 0.5]
    records = [
        {"name": "first", "optimized_atoms": Atoms(symbols, positions=coords)},
        {"name": "rigid-copy", "optimized_atoms": Atoms(symbols, positions=rotated)},
    ]
    assigned, representatives = assign_reference_basins(records, threshold=0.35)
    assert len(representatives) == 1
    assert [row["basin_id"] for row in assigned] == [0, 0]


def test_perturbed_seed_recovery_compares_with_its_own_unperturbed_geometry():
    coords = np.array([
        [0.0, 0.0, 0.0],
        [0.629, 0.629, 0.629],
        [-0.629, -0.629, 0.629],
        [-0.629, 0.629, -0.629],
        [0.629, -0.629, -0.629],
    ])
    records = [{
        "name": "seed",
        "optimized_atoms": Atoms(["C", "H", "H", "H", "H"], positions=coords),
        "reference_atoms": Atoms(["C", "H", "H", "H", "H"], positions=coords.copy()),
        "basin_id": 0,
    }]
    compare_perturbed_seeds(records, threshold=0.35)
    assert records[0]["reference_basin_match_status"] == "preserved"


def test_summary_calculates_micro_and_energy_window_basin_recall(tmp_path):
    run = tmp_path / "one_system"
    (run / "results").mkdir(parents=True)
    (run / "manifest.json").write_text('''{
      "frames": 2,
      "ensemble_sha256": "abc",
      "xTB_version": "GFN2-xTB test",
      "source_xyz_comment_energy_minus_gfn2_final_energy": {
        "maximum_absolute_delta_numeric": 0.0
      }
    }''')
    (run / "reference_basins.json").write_text('''{
      "threshold_sensitivity": {"0.25": {"basin_count": 2}, "0.35": {"basin_count": 2}, "0.50": {"basin_count": 1}},
      "inferred_topology_groups": {"same-graph": 2},
      "structures": [
        {"name": "frame_0000", "basin_id": 0, "xtb_energy_hartree": -1.0},
        {"name": "frame_0001", "basin_id": 1, "xtb_energy_hartree": -0.997}
      ]
    }''')
    condition = {
        "algorithm_requested": "agglomerative",
        "feature_requested": "soap",
        "selected_names": ["frame_0000"],
        "retention_fraction": 0.5,
        "basins_retained": 1,
        "coverage": {
            "reference_basin_count": 2,
            "mean_nearest_indexed_heavy_atom_rmsd_angstrom": 0.2,
            "max_nearest_indexed_heavy_atom_rmsd_angstrom": 0.3,
        },
        "downstream_optimization": {
            "selected_count": 1, "successful_terminations": 1,
            "unique_geometries_optimized_and_cached": 1,
            "verified_reference_basin_seed_count": 1,
            "basin_changed_after_complete_comparison_count": 0,
            "uncertain_due_to_incomplete_mapping_count": 0,
        },
        "diagnostics": {
            "algorithm_fallbacks": [], "feature_fallbacks": [], "distance_fallbacks": [],
            "algorithm_used": "agglomerative", "feature_used": "soap", "distance_used": "euclidean",
        },
        "runtime_seconds": 0.1,
    }
    (run / "results" / "comparison.json").write_text(
        '{"reference_basin_count": 2, "conditions": ['
        + json.dumps(condition) + ']}'
    )
    result = summarize([run], tmp_path / "analysis")
    setting = result["settings"][0]
    assert setting["micro_basin_retention"] == 0.5
    assert setting["energy_window_basin_retention"]["1.0"]["micro_recall"] == 1.0
    assert setting["energy_window_basin_retention"]["3.0"]["micro_recall"] == 0.5
    assert setting["verified_reference_basin_recovery_fraction"] == 1.0


def test_changed_input_cannot_reuse_old_optimizer_output_or_restart(tmp_path):
    import pytest
    xtb = tmp_path / 'xtb'
    xtb.write_text('executable')
    record = {'name': 'frame_0000', 'atoms': Atoms('H2', positions=[[0, 0, 0], [0, 0, .74]])}
    frame = tmp_path / 'optimizations' / record['name']
    frame.mkdir(parents=True)
    (frame / 'xtbopt.xyz').write_text('stale optimized geometry')
    (frame / 'xtbrestart').write_text('stale restart')
    def run(command, *, cwd, **kwargs):
        assert not (Path(cwd) / 'xtbrestart').exists()
        assert not (Path(cwd) / 'xtbopt.xyz').exists()
        return SimpleNamespace(returncode=0, stdout='normal termination of xtb', stderr='')
    with mock.patch('pyar.scripts.scientific_clustering_benchmark.subprocess.run', side_effect=run):
        with pytest.raises(RuntimeError, match='did not terminate normally'):
            optimize_ensemble([record], tmp_path, str(xtb))
    assert not (frame / 'cache.json').exists()
