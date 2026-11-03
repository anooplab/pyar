"""Analytical coverage and bounded report/cache orchestration, without QC."""

from copy import deepcopy
import hashlib
import json
from unittest.mock import Mock
from types import SimpleNamespace

from ase import Atoms
import numpy as np
import pytest

from pyar.scripts import scientific_clustering_benchmark as benchmark


def records():
    result = []
    for index, distance in enumerate([1.0, 3.0]):
        atoms = Atoms("CC", positions=[[-distance / 2, 0, 0], [distance / 2, 0, 0]])
        result.append(
            dict(
                name=f"frame-{index}",
                index=index,
                atoms=atoms,
                optimized_atoms=atoms.copy(),
                basin_id=index,
                xtb_energy_hartree=-1.0 - index,
                source_energy=10.0 + index,
                final_gradient_norm_eh_per_alpha=0.0,
            )
        )
    return result


@pytest.mark.parametrize(
    "selection,covered,mean,maximum",
    [
        ([0], 1, 0.5, 1.0),
        ([0, 1], 2, 0.0, 0.0),
        ([], 0, None, None),
    ],
)
def test_coverage_counts_basins_once_and_uses_heavy_atom_distance(
    selection, covered, mean, maximum
):
    ensemble = records()
    ensemble.append(dict(ensemble[0], name="same-basin"))
    result = benchmark._coverage(ensemble, [ensemble[i] for i in selection], 0.5)
    assert result["reference_basin_count"] == 2
    assert result["covered_basins_at_label_threshold"] == covered
    assert result["basin_recall"] == covered / 2
    expected_mean = pytest.approx(mean) if mean is not None else None
    expected_maximum = pytest.approx(maximum) if maximum is not None else None
    assert result["mean_nearest_indexed_heavy_atom_rmsd_angstrom"] == expected_mean
    assert result["max_nearest_indexed_heavy_atom_rmsd_angstrom"] == expected_maximum


def test_coverage_never_compares_incompatible_ordered_elements():
    ensemble = records()
    selected = {"optimized_atoms": Atoms("NN", positions=[[0, 0, 0], [1, 0, 0]])}
    result = benchmark._coverage(ensemble, [selected], 0.5)
    assert result["covered_basins_at_label_threshold"] == 0
    assert result["mean_nearest_indexed_heavy_atom_rmsd_angstrom"] is None


@pytest.mark.parametrize("downstream", [False, True])
def test_conditions_preserve_budget_labels_random_stream_and_shared_cache(
    tmp_path, monkeypatch, downstream
):
    ensemble = records()
    original = [item["atoms"].positions.copy() for item in ensemble]
    selected = Mock(
        side_effect=lambda molecules, *args: (
            [molecules[0]],
            {"algorithm_used": "fixture", "feature_used": "fixture"},
        )
    )
    monkeypatch.setattr(benchmark, "_cluster_first_select", selected)
    optimized_inputs = []

    def optimize(batch, output, executable):
        optimized_inputs.append(deepcopy(batch))
        assert output == tmp_path / "downstream/shared"
        assert executable == "unused-boundary"
        return [
            dict(item, optimized_atoms=item["reference_atoms"].copy()) for item in batch
        ]

    optimize_mock = Mock(side_effect=optimize)
    monkeypatch.setattr(benchmark, "optimize_ensemble", optimize_mock)
    report = benchmark.run_conditions(
        ensemble,
        [object(), object()],
        tmp_path,
        ["first", "second"],
        ["feature"],
        1,
        0.5,
        downstream,
        0.01,
        "unused-boundary",
    )
    assert report["seed_budget"] == 1
    assert len(report["conditions"]) == 2
    assert selected.call_count == 2
    for call in selected.call_args_list:
        assert call.args[3] == 1
    for condition in report["conditions"]:
        assert condition["selected_names"] == ["frame-0"]
        assert condition["selected_reference_basin_ids"] == [0]
        assert condition["retention_fraction"] == 0.5
        if downstream:
            details = condition["downstream_optimization"]
            assert details["successful_terminations"] == 1
            assert details["verified_reference_basin_recovery_fraction"] == 1.0
            assert details["unique_geometries_optimized_and_cached"] == 1
        else:
            assert condition["downstream_optimization"] is None
    if downstream:
        assert [len(batch) for batch in optimized_inputs] == [1, 0]
        expected = original[0] + np.random.default_rng(90210).normal(
            0.0, 0.01, size=(2, 3)
        )
        np.testing.assert_array_equal(
            optimized_inputs[0][0]["atoms"].positions, expected
        )
    else:
        optimize_mock.assert_not_called()
    for item, coordinates in zip(ensemble, original):
        np.testing.assert_array_equal(item["atoms"].positions, coordinates)
    assert json.loads((tmp_path / "results/comparison.json").read_text()) == report
    assert not (tmp_path / "results/comparison.partial.json").exists()


def test_interrupted_condition_retains_only_completed_checkpoint_then_clean_rerun(
    tmp_path, monkeypatch
):
    ensemble = records()

    def select(molecules, algorithm, *args):
        if algorithm == "interrupt":
            raise RuntimeError("injected interruption")
        return [molecules[0]], {}

    monkeypatch.setattr(benchmark, "_cluster_first_select", select)
    arguments = (ensemble, [object(), object()], tmp_path)
    with pytest.raises(RuntimeError, match="injected interruption"):
        benchmark.run_conditions(
            *arguments,
            ["first", "interrupt"],
            ["feature"],
            1,
            0.5,
            False,
            0.01,
            "unused",
        )
    partial = tmp_path / "results/comparison.partial.json"
    checkpoint = json.loads(partial.read_text())
    assert len(checkpoint["conditions"]) == 1
    assert checkpoint["conditions"][0]["algorithm_requested"] == "first"
    assert not (tmp_path / "results/comparison.json").exists()
    # This API reruns conditions; it does not claim checkpoint-based resume.
    benchmark.run_conditions(
        *arguments, ["first"], ["feature"], 1, 0.5, False, 0.01, "unused"
    )
    assert not partial.exists()
    assert (tmp_path / "results/comparison.csv").is_file()


def test_reference_labels_bind_source_bytes_and_preserve_upstream_energy_interpretation(
    tmp_path,
):
    ensemble, representatives = benchmark.assign_reference_basins(records(), 0.35)
    assert len(representatives) == 2
    source = tmp_path / "source.xyz"
    source.write_bytes(b"synthetic source provenance witness\n")
    before = [(item["basin_id"], item["atoms"].positions.copy()) for item in ensemble]
    target = benchmark.write_reference_labels(
        ensemble, representatives, tmp_path, source, 0.35
    )
    report = json.loads(target.read_text())
    assert report["source_sha256"] == hashlib.sha256(source.read_bytes()).hexdigest()
    assert report["source_file"] == "source.xyz"
    assert report["threshold_angstrom"] == 0.35
    assert set(report["threshold_sensitivity"]) == {"0.25", "0.35", "0.50"}
    assert [basin["members"] for basin in report["basins"]] == [
        ["frame-0"],
        ["frame-1"],
    ]
    assert [row["source_energy"] for row in report["structures"]] == [10.0, 11.0]
    assert all(
        "untrusted upstream" in row["source_energy_interpretation"]
        for row in report["structures"]
    )
    for item, (basin, coordinates) in zip(ensemble, before):
        assert item["basin_id"] == basin
        np.testing.assert_array_equal(item["atoms"].positions, coordinates)


@pytest.mark.parametrize(
    "mutation", ["executable", "missing-output", "corrupt-cache", "output-tamper"]
)
def test_invalid_cache_uses_a_fresh_attempt_and_cannot_adopt_stale_output(
    tmp_path, monkeypatch, mutation
):
    from ase.io import write

    executable = tmp_path / "xtb"
    executable.write_text("synthetic executable identity")
    ensemble = records()[:1]
    attempts = []

    def execute(command, *, cwd, **kwargs):
        attempts.append(cwd)
        assert command[2:] == ["--opt", "--gfn", "2", "--chrg", "0", "--uhf", "0"]
        if len(attempts) == 1:
            write(cwd / "xtbopt.xyz", ensemble[0]["atoms"], format="xyz")
        return SimpleNamespace(
            returncode=0,
            stdout="TOTAL ENERGY -1.0\nGRADIENT NORM 0.0\nnormal termination of xtb\n",
            stderr="",
        )

    monkeypatch.setattr(benchmark.subprocess, "run", execute)
    benchmark.optimize_ensemble(ensemble, tmp_path, str(executable))
    benchmark.optimize_ensemble(ensemble, tmp_path, str(executable))
    assert len(attempts) == 1
    frame = tmp_path / "optimizations/frame-0"
    if mutation == "executable":
        executable.write_text("changed executable identity")
    elif mutation == "missing-output":
        (frame / "xtbopt.xyz").unlink()
    elif mutation == "corrupt-cache":
        (frame / "cache.json").write_text("invalid JSON")
    else:
        (frame / "xtbopt.xyz").write_text("truncated geometry")
    with pytest.raises(RuntimeError, match="did not terminate normally"):
        benchmark.optimize_ensemble(ensemble, tmp_path, str(executable))
    assert len(attempts) == 2
    assert attempts[0] != attempts[1]
    assert (attempts[0] / "xtbopt.xyz").is_file()
    assert not (attempts[1] / "xtbopt.xyz").exists()
