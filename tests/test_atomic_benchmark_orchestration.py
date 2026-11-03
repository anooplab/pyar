"""Small analytical report witnesses with only clustering replaced."""

import csv
import json
import warnings

import pytest

from pyar.backend_errors import BackendExecutionError
from pyar.core.molecule import Molecule
from pyar.scripts import atomic_cluster_benchmark as benchmark


def ensemble():
    molecules = [
        Molecule(
            ["He", "He"],
            [[0.0, 0.0, 0.0], [distance, 0.0, 0.0]],
            name=f"entry-{index}",
            energy=energy,
        )
        for index, (distance, energy) in enumerate(
            [(1.0, -3.0), (2.0, -1.0), (3.0, -2.0)]
        )
    ]
    metadata = [
        {
            "source_id": f"source-{index}",
            "energy_unit": "synthetic-unit",
            "final_force_maximum_ev_per_angstrom": 0.1 if index == 0 else None,
        }
        for index in range(3)
    ]
    return molecules, metadata


def test_report_preserves_source_ranks_budget_and_analytical_distance_proxy(
    tmp_path, monkeypatch
):
    molecules, metadata = ensemble()
    calls = []

    def select(inputs, algorithm, feature, maximum_seeds, **kwargs):
        assert inputs is molecules
        assert maximum_seeds == 2
        assert kwargs == {"system_type": "atomic-cluster"}
        calls.append((algorithm, feature))
        warnings.warn("synthetic selection warning", UserWarning)
        return [inputs[0], inputs[2]], {
            "algorithm_used": "test-selector",
            "feature_used": "test-feature",
            "cluster_minimum_count": 3,
        }

    monkeypatch.setattr(benchmark, "_cluster_first_select", select)
    output = tmp_path / "report"
    report = benchmark.run_benchmark(
        "synthetic",
        molecules,
        metadata,
        output,
        algorithms=["first", "second"],
        features=["a", "b"],
        maximum_seeds=2,
    )
    assert calls == [("first", "a"), ("first", "b"), ("second", "a"), ("second", "b")]
    assert report["seed_budget"] == 2
    assert report["source_quality"]["entry_count_with_final_force_diagnostic"] == 1
    assert (
        report["source_quality"]["maximum_reported_final_force_ev_per_angstrom"] == 0.1
    )
    for condition in report["conditions"]:
        assert condition["selected_source_ids"] == ["source-0", "source-2"]
        assert condition["selected_count"] == 2
        assert condition["lowest_source_energy_entry_retained"] is True
        assert condition["best_selected_source_energy_rank"] == 1
        assert condition["top_five_source_energy_entries_retained"] == 2
        assert condition["selected_energy_min"] == -3.0
        assert condition["selected_energy_max"] == -2.0
        assert condition["energy_unit"] == "synthetic-unit"
        assert condition["mean_nearest_pair_spectrum_rms"] == pytest.approx(1 / 3)
        assert condition["max_nearest_pair_spectrum_rms"] == 1.0
        assert condition["normalized_mean_pair_spectrum_rms"] == pytest.approx(
            1 / 42**0.5
        )
        assert condition["runtime_warnings"] == [
            {"category": "UserWarning", "message": "synthetic selection warning"}
        ]
    assert json.loads((output / "comparison.json").read_text()) == report
    with (output / "structures.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    assert [row["source_id"] for row in rows] == ["source-0", "source-1", "source-2"]
    with (output / "comparison.csv").open() as stream:
        assert len(list(csv.DictReader(stream))) == 4


@pytest.mark.parametrize("empty", [False, True])
def test_misaligned_or_empty_ensembles_fail_before_creating_output(tmp_path, empty):
    molecules, metadata = ensemble()
    with pytest.raises(ValueError, match="nonempty and aligned"):
        benchmark.run_benchmark(
            "invalid",
            [] if empty else molecules,
            [] if empty else metadata[:-1],
            tmp_path / "report",
        )
    assert list(tmp_path.iterdir()) == []


def test_failed_selection_never_publishes_a_complete_report(tmp_path, monkeypatch):
    molecules, metadata = ensemble()

    def fail(*args, **kwargs):
        raise BackendExecutionError("synthetic selection failure")

    monkeypatch.setattr(benchmark, "_cluster_first_select", fail)
    with pytest.raises(BackendExecutionError, match="selection failure"):
        benchmark.run_benchmark(
            "failed",
            molecules,
            metadata,
            tmp_path / "report",
            algorithms=["a"],
            features=["b"],
        )
    assert not (tmp_path / "report/comparison.json").exists()
    assert [molecule.name for molecule in molecules] == [
        "entry-0",
        "entry-1",
        "entry-2",
    ]
