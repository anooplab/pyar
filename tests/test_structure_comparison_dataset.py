"""Integrity and behavior checks for the local pair benchmark fixtures."""

import csv
from collections import Counter
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from pyar.structure_comparison import GraphRMSDComparator


DATASET = Path(__file__).parents[1] / "benchmarks" / "structure_comparison"


def _load_xyz_frames():
    lines = (DATASET / "structures.xyz").read_text().splitlines()
    frames = {}
    index = 0
    while index < len(lines):
        natoms = int(lines[index])
        metadata = dict(
            item.split("=", 1) for item in lines[index + 1].split() if "=" in item
        )
        rows = [line.split() for line in lines[index + 2:index + 2 + natoms]]
        frames[metadata["structure_id"]] = SimpleNamespace(
            atoms_list=[row[0] for row in rows],
            coordinates=np.asarray(
                [[float(value) for value in row[1:4]] for row in rows], dtype=float,
            ),
            metadata=metadata,
        )
        index += natoms + 2
    return frames


def _load_pairs():
    with (DATASET / "pairs.csv").open(newline="") as stream:
        return list(csv.DictReader(stream))


def test_benchmark_pair_manifest_resolves_every_structure_and_class():
    frames = _load_xyz_frames()
    pairs = _load_pairs()
    labels = Counter(pair["label"] for pair in pairs)

    assert len(frames) == 126
    assert len(pairs) == 189
    assert {"atomic_cluster", "molecular_non_covalent_cluster", "distorted_molecule"} <= {
        frame.metadata.get("system_class") for frame in frames.values()
    }
    assert labels["coordinate_perturbation"] == 32
    assert labels["different_cluster_arrangement"] == 2
    assert all(
        pair["structure_a"] in frames and pair["structure_b"] in frames
        for pair in pairs
    )


def test_irmsd_handles_atomic_cluster_permutation_and_cluster_motifs():
    pytest.importorskip("irmsd")
    from pyar.structure_comparison import IRMSDComparator

    frames = _load_xyz_frames()
    comparator = IRMSDComparator(threshold=0.1, inversion=2)
    permuted = comparator.compare(frames["au13_ico"], frames["au13_ico_rotperm"])
    motifs = comparator.compare(frames["au13_ico"], frames["au13_cubocta"])

    assert permuted.equivalent
    assert permuted.distance < 1e-12
    assert motifs.compatible
    assert not motifs.equivalent
    assert motifs.distance > 0.1


def test_irmsd_handles_whole_fragment_permutation_and_water_arrangements():
    pytest.importorskip("irmsd")
    from pyar.structure_comparison import IRMSDComparator

    frames = _load_xyz_frames()
    comparator = IRMSDComparator(threshold=0.1, inversion=2)
    fragment_permutation = comparator.compare(
        frames["water_dimer_hbond"], frames["water_dimer_hbond_fragperm"],
    )
    dimer_arrangement = comparator.compare(
        frames["water_dimer_hbond"], frames["water_dimer_side"],
    )
    trimer_arrangement = comparator.compare(
        frames["water_trimer_ring"], frames["water_trimer_chain"],
    )

    assert fragment_permutation.equivalent
    assert fragment_permutation.distance < 1e-12
    assert dimer_arrangement.distance > 0.1
    assert trimer_arrangement.distance > 0.1


def test_graph_rmsd_handles_permuted_disconnected_water_fragments():
    frames = _load_xyz_frames()
    result = GraphRMSDComparator(threshold=1e-8).compare(
        frames["water_dimer_hbond"], frames["water_dimer_hbond_fragperm"],
    )
    assert result.compatible
    assert result.equivalent
    assert result.distance < 1e-12
    assert result.metadata["connectivity_match"]
    assert result.metadata["comparison_complete"]


def test_graph_rmsd_returns_incomplete_result_for_excessive_cluster_symmetry():
    frames = _load_xyz_frames()
    result = GraphRMSDComparator(max_isomorphisms=10).compare(
        frames["au13_ico"], frames["au13_ico_rotperm"],
    )
    assert result.compatible
    assert result.distance is None
    assert result.equivalent is None
    assert not result.metadata["comparison_complete"]
    assert result.metadata["isomorphisms_evaluated"] == 10


def test_distorted_molecular_pairs_record_source_without_topology_assertion():
    frames = _load_xyz_frames()
    distorted_pairs = [
        pair for pair in _load_pairs()
        if pair["label"] == "coordinate_perturbation"
    ]
    assert len(distorted_pairs) == 32
    assert all(
        frames[pair["structure_b"]].metadata["system_class"] == "distorted_molecule"
        and pair["provenance"].startswith("source_geometry_cartesian_noise_sigma_")
        for pair in distorted_pairs
    )


def _audit_tools():
    import runpy
    return runpy.run_path(str(DATASET / 'systematic_test.py'))


def test_independent_dataset_audit_verifies_geometry_and_label_evidence():
    tools = _audit_tools()
    frames, pairs = tools['load_dataset']()
    result = tools['audit_dataset'](frames, pairs)
    assert result['valid'], result['errors']
    assert sum(p['binary_expected'] is True for p in result['pair_evidence']) == 26
    assert sum(p['binary_expected'] is False for p in result['pair_evidence']) == 73
    assert all(p['binary_expected'] is None for p in result['pair_evidence']
               if p['label'] == 'unreviewed_conformer_pair')


def test_audit_rejects_a_falsely_labelled_rigid_transform():
    tools = _audit_tools()
    frames, pairs = tools['load_dataset']()
    pair = next(p for p in pairs if p['label'] == 'same_structure')
    frames[pair['structure_b']].coordinates[0, 0] += 0.25
    result = tools['audit_dataset'](frames, pairs)
    assert not result['valid']
    assert any('Rigid transform not verified' in error for error in result['errors'])


def test_audit_rejects_same_identity_labelled_as_connectivity_isomer():
    tools = _audit_tools()
    frames, pairs = tools['load_dataset']()
    pair = next(p for p in pairs if p['label'] == 'different_connectivity')
    pair['structure_b'] = pair['structure_a']
    result = tools['audit_dataset'](frames, pairs)
    assert not result['valid']
    assert any('Invalid constitutional isomer label' in error for error in result['errors'])


def test_benchmark_does_not_count_native_backend_fallback_as_success(monkeypatch):
    import json
    import subprocess
    tools = _audit_tools()
    payload = dict(compatible=True, distance=0.01, equivalent=True,
                   metadata={}, seconds=0.1, inputs_unchanged=True)
    response = SimpleNamespace(returncode=0, stderr='', stdout=(
        "WARNING: Falling back to quaternion RMSD\nRESULT_JSON:" + json.dumps(payload) + '\n'
    ))
    monkeypatch.setattr(subprocess, 'run', lambda *args, **kwargs: response)
    case = dict(method='irmsd', case_id='fixture', direction='forward', a='a', b='b')
    assert tools['run_case'](case)['status'] == 'backend_fallback'
