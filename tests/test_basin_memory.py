import json
from types import SimpleNamespace

import numpy as np
import pytest

from pyar.selection.basin_memory import (
    BasinMemoryError,
    REGISTRY_FORMAT,
    REGISTRY_VERSION,
    _apply_basin_memory,
    _load_basin_registry,
    _persist_basin_registry,
    migrate_basin_registry,
)


def molecule(name, atoms, coordinates, energy=0.0, charge=0, multiplicity=1):
    return SimpleNamespace(
        name=name,
        atoms_list=atoms,
        coordinates=np.asarray(coordinates, dtype=float),
        energy=energy,
        charge=charge,
        multiplicity=multiplicity,
    )


def test_persist_writes_geometry_and_versioned_descriptor(tmp_path):
    path = tmp_path / "basin_registry.json"
    geometry = molecule("water", ["O", "H", "H"], [[0, 0, 0], [0.95, 0, 0], [-0.2, 0.9, 0]], -76.0)

    _persist_basin_registry(str(path), [geometry])
    payload = json.loads(path.read_text())

    assert payload["format"] == REGISTRY_FORMAT
    assert payload["schema_version"] == REGISTRY_VERSION
    entry = payload["entries"][0]
    assert entry["geometry"]["atoms"] == ["O", "H", "H"]
    assert entry["geometry"]["charge"] == 0
    descriptor = entry["descriptors"]["geometry"]
    assert descriptor["name"] == "element-pair-distance-histogram"
    assert descriptor["version"] == 1
    assert np.isfinite(descriptor["values"]).all()
    assert "fingerprint" not in entry


def test_memory_compacts_verified_rigid_and_permuted_copies_before_capacity(tmp_path):
    path = tmp_path / "basin_registry.json"
    first = molecule("first", ["H", "H"], [[0, 0, 0], [0.74, 0, 0]])
    rotated_copy = molecule("copy", ["H", "H"], [[3, 2, 1], [3, 2.74, 1]])
    distinct = molecule("distinct", ["H", "H"], [[0, 0, 0], [1.2, 0, 0]])

    _persist_basin_registry(str(path), [first, rotated_copy, distinct], max_entries=2)
    entries = _load_basin_registry(str(path))

    assert [entry["name"] for entry in entries] == ["first", "distinct"]


def test_legacy_migration_is_explicit_and_preserves_opaque_fingerprint(tmp_path):
    path = tmp_path / "basin_registry.json"
    original = {"stoichiometry": "H2", "entries": [{"name": "old", "energy": -1.0,
                                                      "fingerprint": [0.25, 0.75]}]}
    path.write_text(json.dumps(original))

    report = migrate_basin_registry(str(path), dry_run=True)
    assert report["changed"]
    assert json.loads(path.read_text()) == original
    report = migrate_basin_registry(str(path))
    payload = json.loads(path.read_text())

    assert report["entries_with_geometry"] == 0
    assert report["entries_with_opaque_legacy_descriptor"] == 1
    migrated = payload["entries"][0]["legacy_descriptors"]["pyar-fingerprint-v1"]
    assert migrated["values"] == [0.25, 0.75]
    assert migrated["source_uncertain"] is True
    assert "geometry" not in payload["entries"][0]
    assert _load_basin_registry(str(path))[0]["legacy_descriptors"]
    second_report = migrate_basin_registry(str(path))
    assert second_report["changed"] is False


def test_legacy_fingerprint_without_geometry_never_prunes_candidates():
    candidates = [molecule(f"m{i}", ["H"], [[i, 0, 0]], energy=i) for i in range(8)]
    legacy = [{"legacy_descriptors": {"pyar-fingerprint-v1": {"values": [float(i)]}}} for i in range(4)]

    selected = _apply_basin_memory(candidates, 2, legacy)

    assert selected is candidates


def test_geometry_backed_memory_uses_current_distance_backend(tmp_path):
    candidates = [
        molecule(f"m{i}", ["H", "H"], [[0, 0, 0], [0.6 + i * 0.4, 0, 0]], energy=float(i))
        for i in range(8)
    ]
    registry = tmp_path / "registry.json"
    _persist_basin_registry(str(registry), [candidates[0], candidates[-1]])
    entries = _load_basin_registry(str(registry))

    selected = _apply_basin_memory(
        candidates, 2, entries, feature="distance-histogram",
        distance_metric="graph-rmsd", system_type="molecular",
    )

    assert len(selected) == 4
    assert candidates[0] in selected
    assert candidates[1] in selected


def test_auto_memory_feature_is_resolved_from_current_candidates(monkeypatch):
    from types import SimpleNamespace

    import pyar.selection.features as features
    import pyar.selection.structural_distances as distances
    from pyar.selection.basin_memory import _basin_novelty_scores

    candidates = [molecule("c0", ["H", "H"], [[0, 0, 0], [0.7, 0, 0]], 0.0),
                  molecule("c1", ["H", "H"], [[0, 0, 0], [1.0, 0, 0]], 1.0)]
    entry = {"name": "old", "energy": -1.0,
             "geometry": {"atoms": ["H", "H"], "coordinates": [[0, 0, 0], [0.9, 0, 0]]}}
    feature_calls = []

    def compute(pool, requested, **kwargs):
        feature_calls.append((list(pool), requested, kwargs["system_type"]))
        return SimpleNamespace(name="mbtr", system_type="molecular",
                               values=np.arange(len(pool), dtype=float)[:, None])

    monkeypatch.setattr(features, "compute_feature_matrix", compute)
    monkeypatch.setattr(features, "standardize_features", lambda values: values)
    monkeypatch.setattr(
        distances, "compute_distance_matrix",
        lambda pool, *_args, **_kwargs: SimpleNamespace(
            values=np.array([[0, 0.2, 0.4], [0.2, 0, 0.3], [0.4, 0.3, 0]], dtype=float)
        ),
    )

    scores = _basin_novelty_scores(candidates, [entry], feature="auto")

    assert feature_calls[0][0] == candidates
    assert feature_calls[0][1] == "auto"
    assert feature_calls[1][1] == "mbtr"
    assert len(feature_calls[1][0]) == 3
    assert [score[0] for score in scores] == [0.2, 0.4]


def test_failed_memory_distance_keeps_every_candidate(monkeypatch):
    candidates = [molecule(f"m{i}", ["H"], [[i, 0, 0]], energy=i) for i in range(8)]
    entries = [{"geometry": {"atoms": ["H"], "coordinates": [[-1, 0, 0]]}}]
    monkeypatch.setattr(
        "pyar.selection.basin_memory._geometry_novelty_scores",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(RuntimeError("distance failure")),
    )

    assert _apply_basin_memory(candidates, 2, entries) is candidates


@pytest.mark.parametrize("payload", [["bad"], {"schema_version": 3, "entries": []}, {"schema_version": 2, "entries": "bad"}])
def test_migration_rejects_invalid_or_future_registry_without_modifying_file(tmp_path, payload):
    path = tmp_path / "basin_registry.json"
    path.write_text(json.dumps(payload))
    original = path.read_bytes()

    with pytest.raises(BasinMemoryError):
        migrate_basin_registry(str(path))

    assert path.read_bytes() == original


def test_version_two_descriptor_is_not_misidentified_as_legacy(tmp_path):
    path = tmp_path / "basin_registry.json"
    payload = {
        "format": REGISTRY_FORMAT,
        "schema_version": REGISTRY_VERSION,
        "entries": [{"geometry": {"atoms": ["H"], "coordinates": [[0, 0, 0]]},
                     "descriptors": {"geometry": {"name": "new", "version": 1, "values": [0.0]}}}],
    }
    path.write_text(json.dumps(payload))

    entries = _load_basin_registry(str(path))

    assert entries[0]["descriptors"]["geometry"]["name"] == "new"
    assert "legacy_descriptors" not in entries[0]


def test_corrupt_or_future_registry_selection_does_not_overwrite(tmp_path, monkeypatch):
    from pyar.selection import clustering

    selected = tmp_path / "selected" / "stoichiometry_H2"
    selected.mkdir(parents=True)
    path = selected / "basin_registry.json"
    path.write_text('{"schema_version": 99, "entries": []}')
    original = path.read_bytes()
    monkeypatch.chdir(tmp_path)
    candidates = [
        molecule("a", ["H", "H"], [[0, 0, 0], [0.7, 0, 0]], energy=0.0),
        molecule("b", ["H", "H"], [[0, 0, 0], [1.0, 0, 0]], energy=1.0),
    ]
    diagnostics = {}

    result = clustering.choose_geometries(candidates, maximum_number_of_seeds=1, diagnostics=diagnostics)

    assert result == [candidates[0]]
    assert path.read_bytes() == original
    assert diagnostics["basin_memory"]["status"] == "registry-unavailable-keep-all"


def test_atomic_registry_write_keeps_previous_file_when_replace_fails(tmp_path, monkeypatch):
    import pyar.selection.basin_memory as memory

    path = tmp_path / "basin_registry.json"
    path.write_text('{"schema_version": 1, "entries": []}')
    original = path.read_bytes()
    monkeypatch.setattr(memory.os, "replace", lambda *_args: (_ for _ in ()).throw(OSError("disk error")))
    geometry = molecule("H", ["H"], [[0, 0, 0]])

    with pytest.raises(OSError, match="disk error"):
        _persist_basin_registry(str(path), [geometry])

    assert path.read_bytes() == original
    assert sorted(item.name for item in tmp_path.iterdir()) == ["basin_registry.json"]


def test_cli_dry_run_reports_legacy_geometry_limit(tmp_path, capsys):
    from pyar.selection.basin_memory import main

    path = tmp_path / "basin_registry.json"
    path.write_text(json.dumps({"entries": [{"fingerprint": [1.0]}]}))

    main(["--dry-run", str(path)])

    report = json.loads(capsys.readouterr().out)
    assert report["entries_with_geometry"] == 0
    assert report["entries_with_opaque_legacy_descriptor"] == 1
    assert "fingerprint" in json.loads(path.read_text())["entries"][0]


def test_workflow_selection_diagnostics_are_written_atomically(tmp_path, monkeypatch):
    from pyar.workflows._growth import _write_growth_selection_diagnostics

    monkeypatch.chdir(tmp_path)
    report = {"basin_memory": {"schema_version": 2, "distance_used": "graph-rmsd"}}

    _write_growth_selection_diagnostics(report, "selected/selection_diagnostics.json")

    path = tmp_path / "selected" / "selection_diagnostics.json"
    assert json.loads(path.read_text()) == report
    assert list(path.parent.iterdir()) == [path]
