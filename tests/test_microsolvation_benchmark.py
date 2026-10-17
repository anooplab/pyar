"""Scientific fixture and scope checks for the MS0 benchmark."""

import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from benchmarks.microsolvation import analyze


ROOT = Path(__file__).parents[1] / "benchmarks" / "microsolvation"


@pytest.mark.parametrize(("name", "atoms"), [
    ("methane", 5), ("water", 3), ("ammonia", 4), ("formamide", 6), ("acetate", 7),
])
def test_project_generated_fixtures_have_valid_atom_counts(name, atoms):
    symbols = analyze.read_xyz(ROOT / "fixtures" / f"{name}.xyz")
    assert len(symbols) == atoms


def test_primary_comparison_uses_clusters_large_enough_for_comparison():
    protocol = json.loads((ROOT / "protocol.json").read_text())
    assert protocol["primary_comparison_count"] == 4
    assert protocol["minimum_primary_cluster_atoms"] == {
        "methane": 17, "ammonia": 16, "formamide": 18, "acetate": 19,
    }
    assert all(case["water_counts"][-1] >= 4 for case in protocol["cases"])
    assert protocol["cases"][-1]["charge"] == -1


def test_xyz_reader_rejects_truncated_file(tmp_path):
    path = tmp_path / "short.xyz"
    path.write_text("2\ntruncated\nH 0 0 0\n")
    with pytest.raises(ValueError, match="Truncated"):
        analyze.read_xyz(path)


def test_water_count_checks_fragment_order():
    assert analyze.count_water_fragments(["C", "O", "H", "H"], 1) == 1
    with pytest.raises(ValueError, match="fragment order"):
        analyze.count_water_fragments(["C", "H", "O", "H"], 1)


def test_solute_prefix_maps_supported_fixtures():
    for name in ("methane", "ammonia", "formamide", "acetate"):
        symbols = analyze.read_xyz(ROOT / "fixtures" / f"{name}.xyz")
        assert analyze.match_solute_prefix(symbols) == name


def test_multiframe_xyz_reader_preserves_every_candidate(tmp_path):
    path = tmp_path / "ensemble.xyz"
    path.write_text("1\nfirst\nH 0 0 0\n1\nsecond\nH 1 0 0\n")
    frames = analyze.read_xyz_frames(path)
    assert len(frames) == 2
    assert frames[1][0] == ["H"]
    assert "second" in frames[1][1]


def test_common_method_rescore_records_energy_from_calculator(monkeypatch, tmp_path):
    xyz = ROOT / "fixtures" / "methane.xyz"

    def fake_run(command, **kwargs):
        assert command[1:] == ["candidate.xyz", "--gfn", "2", "--sp", "--chrg", "0", "--uhf", "0"]
        return SimpleNamespace(returncode=0, stdout="TOTAL ENERGY -12.345678 Eh\n")

    monkeypatch.setattr(analyze.subprocess, "run", fake_run)
    assert analyze.energy_for(xyz, "xtb", tmp_path / "rescore") == pytest.approx(-12.345678)


def test_common_method_rescore_preserves_anion_charge(monkeypatch, tmp_path):
    xyz = ROOT / "fixtures" / "acetate.xyz"

    def fake_run(command, **kwargs):
        assert command[-4:] == ["--chrg", "-1", "--uhf", "0"]
        return SimpleNamespace(returncode=0, stdout="TOTAL ENERGY -20.123456 Eh\n")

    monkeypatch.setattr(analyze.subprocess, "run", fake_run)
    assert analyze.energy_for(xyz, "xtb", tmp_path / "anion", charge=-1) == pytest.approx(-20.123456)
