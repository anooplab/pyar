"""Scientific selection failures reproduced during the reaction review."""

import csv
import json

import numpy as np
import pytest

from pyar.reaction_analysis import (
    _candidate_record_indices,
    _candidate_exclusion_reasons,
    _first_persistent_transition_index,
    analyse_reaction_trace,
)
from pyar.reaction_trace import ReactionTraceRecorder


def test_legacy_trace_force_screen_uses_uncancelled_vectors():
    records = [
        {"max_force": 0.001, "backend_forces_hartree_per_bohr": [[1., 0., 0.]],
         "backend_energy_hartree": -1.},
        {"max_force": 0.001, "backend_forces_hartree_per_bohr": [[0.01, 0., 0.]],
         "backend_energy_hartree": -1.},
    ]
    assert _candidate_record_indices(records, max_force=0.05) == [1]
    assert _candidate_exclusion_reasons(records, max_force=0.05) == {0: ["max_force"]}


def test_total_force_alone_cannot_establish_physical_force_eligibility():
    records = [{"max_force": 0., "backend_energy_hartree": -1.}]
    with pytest.raises(ValueError, match="excluded every"):
        _candidate_record_indices(records, max_force=0.05)
    assert _candidate_exclusion_reasons(records, max_force=0.05) == {
        0: ["physical_force_unavailable"]
    }


def test_zero_mad_outlier_filter_uses_only_electronic_energy():
    records = [{"backend_energy_hartree": energy, "total_energy_hartree": -10.}
               for energy in (-10., -10., -10., -1.)]
    assert _candidate_record_indices(records, energy_outlier_z=3.5) == [0, 1, 2]


def test_filtered_event_does_not_select_a_different_later_topology():
    records = [{"current_bonds": bonds} for bonds in
               ([], [[0, 1]], [[0, 1]], [[0, 2]], [[0, 1]])]
    assert _first_persistent_transition_index(records, [0, 3, 4]) is None
    assert _first_persistent_transition_index(records, [0, 2]) == 2


@pytest.mark.parametrize("distances", [(2.5, 1.4, 1.38), (2.5, 2.4, 2.3)])
def test_unavailable_transition_has_no_mislabelled_or_stale_xyz(tmp_path, distances):
    recorder = ReactionTraceRecorder(tmp_path)
    for index, distance in enumerate(distances):
        force = 0.01 if index == 0 else 1.0
        physical = np.array([[force, 0., 0.], [-force, 0., 0.]])
        recorder.record(
            symbols=["C", "C"], coordinates_angstrom=[[0., 0., 0.], [distance, 0., 0.]],
            backend_energy_hartree=float(index), bias_energy_hartree=-float(index),
            total_energy_hartree=0., backend_forces_hartree_per_bohr=physical,
            bias_forces_hartree_per_bohr=-physical, total_forces_hartree_per_bohr=np.zeros((2, 3)),
            backend_force_norm=float(np.linalg.norm(physical)),
            bias_force_norm=float(np.linalg.norm(physical)), total_force_norm=0., max_force=0.,
        )
    first = analyse_reaction_trace(tmp_path)
    candidate = tmp_path / "candidate_ts" / "first_topology_change.xyz"
    assert candidate.exists() == (distances[1] < 2.)
    if candidate.exists():
        assert first["first_topology_change_index"] == 1
    filtered = analyse_reaction_trace(tmp_path, max_force=0.05)
    assert filtered["first_topology_change_index"] is None
    assert not candidate.exists()
    metadata = json.loads((tmp_path / "candidate_ts" / "metadata.json").read_text())
    assert not metadata["candidate_files"]["first_topology_change"]["available"]
    assert metadata["candidate_files"]["first_topology_change"]["xyz_file"] is None
    with (tmp_path / "path_summary.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    assert float(rows[1]["backend_max_force"]) == 1.
    assert rows[1]["candidate_exclusion_reason"] == "max_force"
