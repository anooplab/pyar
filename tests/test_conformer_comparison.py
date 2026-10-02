"""Scientific invariants and abstention checks for conformer deduplication."""

import json
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import numpy as np
import pytest

from pyar.conformer.comparison import ConformerComparisons
from pyar.structure_comparison import ComparisonResult, GraphFirstDeduplicationComparator
from pyar.structure_comparison.graph_rmsd import GraphRMSDComparator, infer_molecular_graph
from pyar.workflows import conformer


def record(symbols, coordinates, *, energy=0, index=0, charge=0, multiplicity=1):
    result = conformer.ConformerRecord(1, index, energy, "converged", "mmff")
    result.molecule = SimpleNamespace(
        atoms_list=list(symbols), coordinates=np.asarray(coordinates, dtype=float),
        charge=charge, multiplicity=multiplicity,
    )
    return result


def rdkit_record(smiles, seed=7):
    Chem = pytest.importorskip("rdkit.Chem")
    AllChem = pytest.importorskip("rdkit.Chem.AllChem")
    molecule = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(molecule, randomSeed=seed) >= 0
    AllChem.MMFFOptimizeMolecule(molecule)
    return conformer.ConformerRecord(1, 0, 0, "converged", "mmff", rdkit_molecule=molecule)


@pytest.mark.parametrize("atom_mode", ["heavy", "all"])
@pytest.mark.parametrize("smiles", ["CCCC", "CCO", "CC(C)C", "c1ccccc1"])
def test_rotated_translated_permuted_duplicates_keep_lowest_energy(smiles, atom_mode):
    comparisons = ConformerComparisons(atom_mode)
    source = comparisons._molecule(rdkit_record(smiles))
    rng = np.random.default_rng(9)
    order = rng.permutation(len(source.atoms_list))
    rotation, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    rotation[:, -1] *= np.linalg.det(rotation)
    first = record(source.atoms_list, source.coordinates, energy=-1)
    second = record(np.asarray(source.atoms_list)[order],
                    source.coordinates[order] @ rotation + [4, -2, 3], energy=1, index=1)
    for left, right in ((first, second), (second, first)):
        result = comparisons.compare(left, right, .01)
        assert result.equivalent is True
        assert result.distance < 1e-10
    assert conformer._collapse_generation_records([second, first], .01, comparisons=comparisons) == [first]
    assert conformer._select_unique_records([first, second], 2, .01, comparisons=comparisons) == [first]


def test_butane_anti_and_gauche_are_distinct_conformers():
    source = rdkit_record("CCCC")
    Chem = pytest.importorskip("rdkit.Chem")
    transforms = pytest.importorskip("rdkit.Chem.rdMolTransforms")
    anti = Chem.Mol(source.rdkit_molecule)
    gauche = Chem.Mol(anti)
    transforms.SetDihedralDeg(anti.GetConformer(), 0, 1, 2, 3, 180)
    transforms.SetDihedralDeg(gauche.GetConformer(), 0, 1, 2, 3, 60)
    records = [conformer.ConformerRecord(1, 0, energy, "converged", "mmff", rdkit_molecule=mol)
               for energy, mol in enumerate((anti, gauche))]
    result = ConformerComparisons().compare(*records, .5)
    assert result.compatible
    assert result.equivalent is False
    assert conformer._collapse_generation_records(records, .5) == records


def test_heavy_mode_ignores_hydrogen_geometry_but_all_mode_retains_it():
    comparisons = ConformerComparisons()
    source = comparisons._molecule(rdkit_record("CCO"))
    coordinates = source.coordinates.copy()
    graph = infer_molecular_graph(source)
    oxygen = source.atoms_list.index("O")
    hydrogen = next(node for node in graph.neighbors(oxygen) if source.atoms_list[node] == "H")
    # Rotate only the O-H bond around the C-O axis; no atom correspondence is changed.
    carbon = next(node for node in graph.neighbors(oxygen) if source.atoms_list[node] == "C")
    axis = coordinates[oxygen] - coordinates[carbon]
    axis /= np.linalg.norm(axis)
    vector = coordinates[hydrogen] - coordinates[oxygen]
    coordinates[hydrogen] = coordinates[oxygen] + 2 * np.dot(vector, axis) * axis - vector
    first = record(source.atoms_list, source.coordinates)
    second = record(source.atoms_list, coordinates, index=1)
    assert comparisons.compare(first, second, .01).equivalent is True
    assert ConformerComparisons("all").compare(first, second, .01).equivalent is False
    assert conformer._deduplicate_records([first, second], .01) == [first]
    assert conformer._deduplicate_records([first, second], .01, atom_mode="all") == [first, second]


def test_heavy_mode_preserves_hydrogen_attachment_constraints():
    # The heavy C-O graph and coordinates match, but H is on a different element.
    first = record(["C", "O", "H"], [[0, 0, 0], [1.4, 0, 0], [-1, 0, 0]])
    second = record(["C", "O", "H"], [[0, 0, 0], [1.4, 0, 0], [2.4, 0, 0]], index=1)
    result = ConformerComparisons().compare(first, second, .5)
    assert not result.compatible
    assert conformer._deduplicate_records([first, second], .5) == [first, second]


def test_terminal_hydrogen_symmetry_does_not_exhaust_heavy_mapping_budget():
    source = ConformerComparisons()._molecule(rdkit_record("CCCCCCCC"))
    result = GraphRMSDComparator(threshold=.01, atom_mode="heavy", max_isomorphisms=2).compare(source, source)
    assert result.equivalent is True
    assert result.metadata["isomorphisms_evaluated"] == 2


@pytest.mark.parametrize("smiles", ["CC", "CCO", "c1ccccc1"])
def test_heavy_mapping_matches_exhaustive_full_graph_oracle(smiles):
    import networkx as nx
    from pyar.structure_comparison.graph_rmsd import selected_atom_indices
    from pyar.structure_comparison.rmsd import kabsch_rmsd

    first = ConformerComparisons()._molecule(rdkit_record(smiles))
    second = ConformerComparisons()._molecule(rdkit_record(smiles, seed=11))
    matcher = nx.algorithms.isomorphism.GraphMatcher(
        infer_molecular_graph(first), infer_molecular_graph(second),
        node_match=lambda left, right: left["element"] == right["element"],
    )
    indices = selected_atom_indices(first.atoms_list, "heavy")
    distances = [kabsch_rmsd(first.coordinates[indices],
                             second.coordinates[[mapping[index] for index in indices]])
                 for mapping in matcher.isomorphisms_iter()]
    assert distances
    result = GraphRMSDComparator(atom_mode="heavy").compare(first, second)
    assert result.distance == pytest.approx(min(distances), abs=1e-10)


@pytest.mark.parametrize("symbols,first_coords,second_coords", [
    (["C", "C", "H"], [[0, 0, 0], [1.8, 0, 0], [.9, 0, 0]],
     [[0, 0, 0], [1.8, 0, 0], [-1, 0, 0]]),
    (["C", "O", "H", "H"], [[0, 0, 0], [1.4, 0, 0], [5, 0, 0], [5.6, 0, 0]],
     [[0, 0, 0], [1.4, 0, 0], [5, 0, 0], [8, 0, 0]]),
])
def test_heavy_mapping_keeps_bridging_hydrogens_and_hydrogen_components(symbols, first_coords, second_coords):
    first = record(symbols, first_coords)
    second = record(symbols, second_coords, index=1)
    assert ConformerComparisons().compare(first, second, .1).compatible is False


@pytest.mark.parametrize("atom_mode", ["heavy", "all"])
def test_proper_rotation_does_not_merge_enantiomers(atom_mode):
    symbols = ["C", "F", "Cl", "Br", "H"]
    coordinates = np.array([[0, 0, 0], [.8, .8, .8], [-1, -1, 1],
                            [-1.1, 1.1, -1.1], [.6, -.6, -.6]])
    first = record(symbols, coordinates)
    second = record(symbols, coordinates * [-1, 1, 1], index=1)
    assert ConformerComparisons(atom_mode).compare(first, second, .1).equivalent is False
    assert conformer._deduplicate_records([first, second], .1, atom_mode=atom_mode) == [first, second]


@pytest.mark.parametrize("atom_mode", ["heavy", "all"])
def test_hydrogen_only_systems_use_all_atoms(atom_mode):
    first = record(["H", "H"], [[0, 0, 0], [.74, 0, 0]])
    # Both bond lengths remain on the same side of the coordinate graph cutoff.
    second = record(["H", "H"], [[0, 0, 0], [.83, 0, 0]], index=1)
    result = ConformerComparisons(atom_mode).compare(first, second, .01)
    assert result.equivalent is False
    assert result.distance == pytest.approx(.045)


@pytest.mark.parametrize("change", ["charge", "multiplicity", "connectivity", "composition"])
def test_incompatible_structures_are_retained(change):
    first = record(["C", "O"], [[0, 0, 0], [1.4, 0, 0]])
    second = record(["C", "O"], [[0, 0, 0], [1.4, 0, 0]], index=1)
    if change in {"charge", "multiplicity"}:
        setattr(second.molecule, change, getattr(first.molecule, change) + 1)
    elif change == "connectivity":
        second.molecule.coordinates[1, 0] = 3
    else:
        second.molecule.atoms_list[1] = "N"
    assert conformer._deduplicate_records([first, second], 10) == [first, second]


@pytest.mark.parametrize("failure", ["exception", "abstention", "missing", "nonfinite"])
def test_in_doubt_keep_with_observable_diagnostics(failure):
    first = record(["C", "O"], [[0, 0, 0], [1.4, 0, 0]])
    second = record(["C", "O"], [[0, 0, 0], [1.4, 0, 0]], index=1)
    comparisons = ConformerComparisons()
    if failure == "missing":
        second.molecule = None
    elif failure == "nonfinite":
        second.molecule.coordinates[0, 0] = float("nan")
    kwargs = {}
    if failure == "exception":
        kwargs["side_effect"] = RuntimeError("failed comparison")
    elif failure == "abstention":
        kwargs["return_value"] = ComparisonResult(True, None, None, "graph-rmsd", .1)
    if kwargs:
        with mock.patch.object(GraphFirstDeduplicationComparator, "compare", **kwargs):
            selected = conformer._deduplicate_records([first, second], .1, comparisons=comparisons)
    else:
        selected = conformer._deduplicate_records([first, second], .1, comparisons=comparisons)
    assert selected == [first, second]
    assert comparisons.summary()["stages"]["deduplication"]["uncertain"] == 1
    if failure in {"missing", "nonfinite"}:
        assert conformer._select_diverse_record([second], [first], comparisons=comparisons)[1] is None
    json.dumps(comparisons.summary(), allow_nan=False)


def test_zero_threshold_disables_deduplication():
    first = record(["C"], [[0, 0, 0]])
    assert conformer._deduplicate_records([first, first], 0) == [first, first]


@pytest.mark.parametrize("threshold", [float("nan"), float("inf"), -.1])
def test_invalid_threshold_fails_explicitly(threshold):
    with pytest.raises(ValueError, match="finite and nonnegative"):
        conformer._deduplicate_records([], threshold)


def test_heavy_irmsd_fallback_uses_heavy_score_and_records_verified_mapping(tmp_path, monkeypatch):
    pytest.importorskip("irmsd")
    source = ConformerComparisons()._molecule(rdkit_record("CCCC"))
    # The native worker must import the installed policy outside the checkout.
    monkeypatch.chdir(tmp_path)
    comparator = GraphFirstDeduplicationComparator(threshold=.01, atom_mode="heavy", max_isomorphisms=1)
    result = comparator.compare(source, source)
    assert result.equivalent is True
    assert result.metadata["fallback_mapping_verified"] is True
    assert result.metadata["atom_mode"] == "heavy"
    assert result.distance < 1e-8


@pytest.mark.parametrize("atom_mode", ["heavy", "all"])
def test_real_workflow_persists_policy_and_effective_thresholds(tmp_path, atom_mode):
    pytest.importorskip("rdkit")
    result = conformer.conformer_search(
        "CCCC", num_conformers=5, num_seeds=2, top_n=3,
        torsion_rounds=1, torsion_kicks_per_conformer=2,
        dedup_atom_mode=atom_mode, num_threads=1, root_directory=tmp_path,
    )
    state = json.loads((tmp_path / "conformers/state.json").read_text())
    assert state["version"] == 3
    assert state["request"]["dedup_atom_mode"] == atom_mode
    assert state["request"]["generation_dedup_rms"] == .5
    assert state["request"]["final_dedup_rms"] == .5
    policy = state["comparison_policy"]
    assert policy == result.metadata["comparison_policy"]
    assert policy["atom_mode"] == atom_mode
    assert policy["policy_version"] == 1
    assert policy["uncertainty_policy"] == "keep"
    assert policy["stages"]["generation"]["comparisons"] > 0
    assert policy["stages"]["generation"].get("uncertain", 0) == 0
    for path in result.selected_paths:
        assert Path(path).is_file()


def test_conformer_commands_forward_atom_mode():
    from pyar.scripts import conformer as command
    from pyar.scripts import conformer_benchmark as benchmark_command

    args = benchmark_command.argument_parse(["benchmark.json", "--dedup-atom-mode", "all"])
    assert benchmark_command._conformer_options(args)["dedup_atom_mode"] == "all"
    assert command.argument_parse(["CCO"]).dedup_atom_mode == "heavy"
    result = SimpleNamespace(status="completed", run_directory="run", selected_paths=())
    with mock.patch.object(command, "conformer_search", return_value=result) as search:
        command.main(["CCO", "--dedup-atom-mode", "all"])
    assert search.call_args.kwargs["dedup_atom_mode"] == "all"
