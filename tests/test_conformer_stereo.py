"""Conformer search must not silently mix or invert stereoisomers."""
import pytest

Chem = pytest.importorskip("rdkit.Chem")
AllChem = pytest.importorskip("rdkit.Chem.AllChem")

from pyar.conformer.stereo import conformer_coordinates, describe_input_stereo, stereo_signature  # noqa: E402
from pyar.workflows import conformer  # noqa: E402


def _embedded(smiles, count=6):
    molecule = Chem.AddHs(Chem.MolFromSmiles(smiles), addCoords=True)
    ids = list(AllChem.EmbedMultipleConfs(molecule, count, randomSeed=11))
    return molecule, ids


def test_unspecified_elements_are_reported_without_reference():
    molecule, _ = _embedded("CC(O)C=CC", 1)
    stereo = describe_input_stereo(molecule)
    assert stereo['available'] and stereo['reference'] is None
    assert [element['label'] for element in stereo['unspecified']] == ['C2 centre', 'C4=C5 double bond']


def test_conformer_search_stops_on_unspecified_stereochemistry(tmp_path):
    with pytest.raises(conformer.ConformerWorkflowError, match=r'not fully specified \(C2 centre\)'):
        conformer.conformer_search("CC(O)CC", num_conformers=4, top_n=2, num_seeds=1,
                                   torsion_kicks=False, root_directory=str(tmp_path))
    assert not list(tmp_path.rglob('*.xyz'))


def test_conformer_search_accepts_specified_stereochemistry(tmp_path):
    result = conformer.conformer_search("C[C@H](O)C[C@@H](C)Cl", num_conformers=6, top_n=2, num_seeds=1,
                                        torsion_kicks=False, root_directory=str(tmp_path))
    assert result.status == 'completed'
    assert result.metadata['stereochemistry']['reference'] == 'C[C@H](O)C[C@@H](C)Cl'


def test_specified_input_keeps_its_reference_through_embedding():
    molecule, ids = _embedded("C[C@H](O)CC")
    stereo = describe_input_stereo(molecule)
    assert stereo['unspecified'] == []
    signatures = {stereo_signature(molecule, conformer_coordinates(molecule, cid)) for cid in ids}
    assert signatures == {stereo['reference']}


def test_inverted_geometry_is_detected():
    molecule, ids = _embedded("C[C@H](O)CC", 1)
    coordinates = conformer_coordinates(molecule, ids[0])
    mirrored = coordinates * [1.0, 1.0, -1.0]
    assert stereo_signature(molecule, mirrored) != stereo_signature(molecule, coordinates)


def test_records_contradicting_a_specified_input_are_discarded():
    molecule, ids = _embedded("C[C@H](O)CC", 2)
    reference = describe_input_stereo(molecule)['reference']
    mirror = Chem.Mol(molecule)
    conf = mirror.GetConformer(ids[1])
    for index in range(mirror.GetNumAtoms()):
        position = conf.GetAtomPosition(index)
        conf.SetAtomPosition(index, (position.x, position.y, -position.z))
    records = [
        conformer.ConformerRecord(seed=1, source_conf_id=ids[0], rdkit_energy=0.0, rdkit_status='converged',
                                  force_field='mmff', rdkit_molecule=molecule, rdkit_conf_id=ids[0]),
        conformer.ConformerRecord(seed=1, source_conf_id=ids[1], rdkit_energy=0.1, rdkit_status='converged',
                                  force_field='mmff', rdkit_molecule=mirror, rdkit_conf_id=ids[1]),
    ]
    kept, rejected = conformer._assign_stereo_signatures(records, molecule, reference)
    assert rejected == 1 and kept == records[:1]
    assert kept[0].stereo_signature == reference


def test_unavailable_inputs_never_drop_records():
    assert describe_input_stereo(object())['available'] is False
    assert stereo_signature(None, [[0, 0, 0]]) is None
