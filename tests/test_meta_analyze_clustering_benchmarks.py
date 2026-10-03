import numpy as np

from pyar.scripts.meta_analyze_clustering_benchmarks import (
    _binary_metrics,
    _conservative_operating_point,
)


def test_binary_metrics_returns_none_for_undefined_rates():
    metrics = _binary_metrics(0, 0, 0, 0)

    assert metrics["tp"] == 0
    assert metrics["sensitivity"] is None
    assert metrics["specificity"] is None
    assert metrics["precision"] is None
    assert metrics["false_positive_rate"] is None
    assert metrics["false_discovery_rate"] is None


def test_false_positive_rate_and_false_discovery_rate_use_distinct_denominators():
    metrics = _binary_metrics(tp=3, fp=2, fn=12, tn=259)

    assert metrics["false_positive_rate"] == 2 / 261
    assert metrics["false_discovery_rate"] == 2 / 5
    assert metrics["precision"] == 3 / 5


def test_conservative_water_threshold_meets_specificity_and_maximizes_sensitivity():
    reference = np.array([0.2, 0.7, 0.8, 1.1])
    feature_distance = np.array([0.1, 0.4, 0.2, 0.3])

    result = _conservative_operating_point(reference, feature_distance, 1.0)

    assert result["specificity"] == 1.0
    assert result["sensitivity"] == 0.5
    assert (result["tp"], result["fp"], result["fn"], result["tn"]) == (1, 0, 1, 2)


def test_conservative_water_threshold_returns_none_when_target_is_unattainable():
    reference = np.array([0.2, 0.7, 0.8, 1.1])
    feature_distance = np.array([0.1, 0.4, 0.2, 0.3])

    assert _conservative_operating_point(reference, feature_distance, 1.1) is None


def test_reference_labels_reject_duplicate_frame_names(tmp_path):
    import json
    import pytest
    from pyar.scripts.meta_analyze_clustering_benchmarks import _same_basin_labels
    path = tmp_path / 'reference.json'
    path.write_text(json.dumps({'structures': [{'name': 'frame_0000', 'basin_id': 0},
                                               {'name': 'frame_0000', 'basin_id': 1}]}))
    with pytest.raises(ValueError, match='do not cover'):
        _same_basin_labels(path, 1)


def test_meta_analysis_rejects_changed_source_before_pairing_labels(tmp_path, monkeypatch):
    import json
    import pytest
    from pyar.scripts import meta_analyze_clustering_benchmarks as analysis
    monkeypatch.setattr(analysis, 'CONFORMERS', ('fixture',))
    (tmp_path / 'mpconf196gen').mkdir()
    (tmp_path / 'mpconf196gen/fixture_crest_conformers.xyz').write_text('changed')
    (tmp_path / 'fixture').mkdir()
    (tmp_path / 'fixture/manifest.json').write_text(json.dumps({'ensemble_sha256': 'old'}))
    with pytest.raises(ValueError, match='Source ensemble hash mismatch'):
        analysis.analyze_conformers(tmp_path, tmp_path)
