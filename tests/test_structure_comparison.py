"""Characterization and contract tests for structure comparison services."""

from types import SimpleNamespace
import ast
from pathlib import Path
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
    with mock.patch("pyar.structure_comparison.identity.babel.make_inchi_string_from_xyz", return_value="fixture-inchi"), \
        mock.patch("pyar.structure_comparison.identity.babel.make_smile_string_from_xyz", return_value="fixture-smiles"):
        result = provider.identify("fixture.xyz")
    assert result.canonical_key == "fixture-inchi"
    assert result.inchi == "fixture-inchi"
    assert result.smiles == "fixture-smiles"


def test_reaction_identity_wrapper_matches_provider_result():
    from pyar.reaction_identity import molecule_identity_from_xyz

    with mock.patch("pyar.structure_comparison.identity.babel.make_inchi_string_from_xyz", return_value="fixture-inchi"), \
        mock.patch("pyar.structure_comparison.identity.babel.make_smile_string_from_xyz", return_value="fixture-smiles"):
        provider_result = OpenBabelIdentityProvider().identify("fixture.xyz")
        legacy_result = molecule_identity_from_xyz("fixture.xyz")

    assert legacy_result == {
        "inchi": provider_result.inchi,
        "smiles": provider_result.smiles,
    }


def test_legacy_comparator_is_translation_and_rotation_invariant():
    first = _molecule("a", ["C", "H", "H"], [[0, 0, 0], [1, 0, 0], [0, 1, 0]])
    second = _molecule("b", ["C", "H", "H"], [[3, -2, 1], [3, -1, 1], [2, -2, 1]])
    comparator = LegacyRMSDComparator(threshold=0.01)
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


def test_rmsd_primitive_exact_permutation_and_legacy_wrapper_match():
    from pyar.selection import deduplication
    from pyar.structure_comparison.rmsd import rmsd_after_alignment

    first = _molecule("a", ["C", "H", "H"], [[0, 0, 0], [1, 0, 0], [0, 1, 0]])
    permuted = _molecule("b", ["H", "C", "H"], [[0, 1, 0], [0, 0, 0], [1, 0, 0]])
    assert rmsd_after_alignment(permuted, first) == deduplication._rmsd_after_alignment(permuted, first)
    assert rmsd_after_alignment(permuted, first) < 1e-12


def test_rmsd_primitive_exercises_hungarian_fallback():
    from pyar.structure_comparison.rmsd import exact_element_orders, rmsd_after_alignment

    coordinates = np.asarray([
        [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 2.0, 0.0],
        [0.0, 0.0, 3.0], [1.0, 2.0, 3.0], [2.0, 1.0, 0.5],
        [3.0, 0.5, 1.0],
    ])
    first = _molecule("a", ["H"] * 7, coordinates)
    permuted = _molecule("b", ["H"] * 7, coordinates[::-1])
    assert exact_element_orders(first.atoms_list, permuted.atoms_list) is None
    assert rmsd_after_alignment(permuted, first) < 1e-12


def test_legacy_comparator_prefilter_rejection_is_unchanged():
    first = _molecule("a", ["H", "H"], [[0, 0, 0], [1, 0, 0]])
    second = _molecule("b", ["H", "H"], [[0, 0, 0], [1, 0, 0]])
    comparator = LegacyRMSDComparator(
        threshold=0.1, fingerprint_distance=lambda _a, _b: 1.0,
    )
    result = comparator.compare(first, second)
    assert result.compatible
    assert result.distance is None
    assert result.equivalent is False
    assert result.metadata["prefiltered"] is True


def test_structure_comparison_has_no_selection_or_reaction_identity_imports():
    package = Path(__file__).parents[1] / "pyar" / "structure_comparison"
    forbidden = {"pyar.selection", "pyar.reaction_identity"}
    for path in package.glob("*.py"):
        tree = ast.parse(path.read_text(), filename=str(path))
        imported = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imported.update(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module:
                imported.add(node.module)
        assert not any(
            name == blocked or name.startswith(blocked + ".")
            for name in imported for blocked in forbidden
        ), f"forbidden dependency in {path.name}: {imported & forbidden}"
