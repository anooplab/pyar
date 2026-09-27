"""Characterization and contract tests for structure comparison services."""

from types import SimpleNamespace
from unittest import mock

import numpy as np

from pyar.structure_comparison import (
    ComparisonResult,
    IdentityResult,
    LegacyRMSDComparator,
    OpenBabelIdentityProvider,
)


def _molecule(name, atoms, coordinates):
    return SimpleNamespace(
        name=name, atoms_list=atoms, coordinates=np.asarray(coordinates, dtype=float),
    )


def test_identity_results_use_inchi_as_legacy_canonical_key():
    provider = OpenBabelIdentityProvider()
    first = {"inchi": "species-a", "smiles": "C#N"}
    equivalent = {"inchi": "species-a", "smiles": "N#C"}
    distinct = {"inchi": "species-b", "smiles": "C#N"}

    assert provider.same_identity(first, equivalent)
    assert not provider.same_identity(first, distinct)
    result = IdentityResult("species-a", "test", metadata={"source": "fixture"})
    assert result.canonical_key == "species-a"
    assert result.metadata["source"] == "fixture"


def test_openbabel_provider_wraps_existing_xyz_identity():
    provider = OpenBabelIdentityProvider()
    with mock.patch("pyar.reaction_identity.molecule_identity_from_xyz", return_value={
        "inchi": "fixture-inchi", "smiles": "fixture-smiles",
    }):
        result = provider.identify("fixture.xyz")
    assert result.canonical_key == "fixture-inchi"
    assert result.inchi == "fixture-inchi"
    assert result.smiles == "fixture-smiles"


def test_legacy_comparator_is_translation_and_rotation_invariant():
    first = _molecule("a", ["C", "H", "H"], [[0, 0, 0], [1, 0, 0], [0, 1, 0]])
    second = _molecule("b", ["C", "H", "H"], [[3, -2, 1], [3, -1, 1], [2, -2, 1]])
    comparator = LegacyRMSDComparator(threshold=0.01)
    with mock.patch("pyar.selection.clustering.calc_fingerprint_distance", return_value=0.0):
        result = comparator.compare(first, second)
    assert result.compatible
    assert result.equivalent
    assert result.distance < 1e-12
    assert result.method == comparator.method


def test_legacy_comparator_reports_chemically_incompatible_structures():
    first = _molecule("a", ["C", "H"], [[0, 0, 0], [1, 0, 0]])
    second = _molecule("b", ["C", "O"], [[0, 0, 0], [1, 0, 0]])
    result = LegacyRMSDComparator(threshold=0.1).compare(first, second)
    assert not result.compatible
    assert result.distance is None
    assert result.equivalent is None


def test_comparison_result_metadata_is_immutable():
    source = {"mode": "legacy"}
    result = ComparisonResult(True, 0.0, True, "fixture", metadata=source)
    source["mode"] = "changed"
    assert result.metadata["mode"] == "legacy"
