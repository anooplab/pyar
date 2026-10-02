import pytest
from types import SimpleNamespace

from pyar.scripts.water_cluster_similarity_benchmark import (
    load_water_minima,
    stratified_sample,
)


def _geometry_frame(energy, shift=0.0):
    atoms = [
        ("O", 0.0, 0.0, 0.0), ("H", 0.96, 0.0, 0.0), ("H", -0.24, 0.93, 0.0),
        ("O", 5.0, 0.0, 0.0), ("H", 5.96, 0.0, 0.0), ("H", 4.76, 0.93, 0.0),
        ("O", 10.0, 0.0, 0.0), ("H", 10.96, 0.0, 0.0), ("H", 9.76, 0.93, 0.0),
        ("O", 15.0, 0.0, 0.0), ("H", 15.96, 0.0, 0.0), ("H", 14.76, 0.93, 0.0),
        ("O", 20.0, 0.0, 0.0), ("H", 20.96, 0.0, 0.0), ("H", 19.76, 0.93, 0.0),
        ("O", 25.0, 0.0, 0.0), ("H", 25.96, 0.0, 0.0), ("H", 24.76, 0.93, 0.0),
    ]
    rows = ["18", f"Ord_Energy {energy}"]
    rows.extend(f"{symbol} {x + shift} {y} {z}" for symbol, x, y, z in atoms)
    return "\n".join(rows)


def test_load_water_minima_parses_energy_and_fixed_composition(tmp_path):
    source = tmp_path / "w6.txt"
    source.write_text(_geometry_frame(-10.0) + "\n" + _geometry_frame(-9.0, 1.0) + "\n")
    records = load_water_minima(source)
    assert len(records) == 2
    assert records[0]["energy_kcal_per_mol"] == -10.0
    assert records[1]["molecule"].atoms_list.count("O") == 6
    assert records[1]["source_index"] == 1


def test_stratified_sample_is_deterministic_and_spans_energy_range():
    records = [{"source_index": i} for i in range(10)]
    selected = stratified_sample(records, 4)
    assert [row["source_index"] for row in selected] == [0, 3, 6, 9]
    assert selected == stratified_sample(records, 4)


@pytest.mark.parametrize("sample_size", [0, 1, 11])
def test_stratified_sample_rejects_invalid_size(sample_size):
    with pytest.raises(ValueError):
        stratified_sample([{"source_index": i} for i in range(10)], sample_size)


def test_water_pair_cache_checks_source_and_comparator_identity(tmp_path, monkeypatch):
    import numpy as np
    import pyar.scripts.water_cluster_similarity_benchmark as benchmark

    source = tmp_path / "water.xyz"
    source.write_text(_geometry_frame(-10.0) + "\n" + _geometry_frame(-9.0, 1.0) + "\n")
    comparisons = []

    class Comparator:
        def __init__(self, **kwargs):
            self.options = kwargs

        def compare(self, *_args):
            comparisons.append(self.options)
            return SimpleNamespace(
                compatible=True, distance=float(len(comparisons)),
                metadata={"comparison_complete": True, "rotation_hypotheses": 1},
            )

    monkeypatch.setattr(benchmark, "FragmentRMSDComparator", Comparator)
    monkeypatch.setattr(benchmark, "evaluate_feature_distances", lambda distances, _molecules, features: [
        {"feature_requested": feature, "feature_used": feature, "feature_fallbacks": [],
         "spearman_rank_correlation_with_fragment_rmsd": None,
         "nearest_pair_overlap_at_10_percent": 0, "nearest_pair_precision_at_10_percent": 0,
         "pairwise_distance_matrix": np.zeros_like(distances)}
        for feature in features
    ])
    output = tmp_path / "water-cache"

    benchmark.build_benchmark(source, output, sample_size=2, features=("mbtr",))
    benchmark.build_benchmark(source, output, sample_size=2, features=("mbtr",))
    assert len(comparisons) == 1

    source.write_text(_geometry_frame(-10.0) + "\n" + _geometry_frame(-9.0, 2.0) + "\n")
    benchmark.build_benchmark(source, output, sample_size=2, features=("mbtr",))
    assert len(comparisons) == 2
