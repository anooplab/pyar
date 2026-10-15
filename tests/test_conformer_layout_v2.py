"""Real seeded RDKit compatibility and explicit v2 result discovery."""
from pathlib import Path
from unittest.mock import patch

import pytest

from pyar.workflows import conformer
from pyar.conformer.layout_v2 import selected_artifacts
from pyar.persistence import read_json
from pyar.benchmarks.conformer import ConformerBenchmarkCase, diagnose_conformer_case


@pytest.fixture
def runs(tmp_path):
    pytest.importorskip('rdkit')
    options = dict(num_conformers=8, top_n=2, num_seeds=2, seed=13,
                   torsion_kicks=True, torsion_rounds=1, num_threads=1)
    old = conformer._conformer_search('CCCC', root_directory=tmp_path / 'old', **options)
    new = conformer.conformer_search('CCCC', root_directory=tmp_path / 'new', **options)
    return old, new, options


def test_real_generation_selection_bytes_and_scientific_metadata_match(runs):
    old, new, _ = runs
    assert [Path(path).read_bytes() for path in old.selected_paths] == [Path(path).read_bytes() for path in new.selected_paths]
    before, after = read_json(old.state_path), read_json(new.state_path)
    assert before['comparison_policy'] == after['comparison_policy']
    fields = ('seed', 'source_conf_id', 'force_field', 'rdkit_energy', 'generation_stage',
              'generation_round', 'parent_name', 'torsion_moves', 'refinement_reason', 'backend_energy')
    assert [[row[key] for key in fields] for row in before['records']] == [[row[key] for key in fields] for row in after['records']]
    root = Path(new.run_directory)
    generated = read_json(root / 'generated/records.json')['conformers']
    assert len({row['conformer_id'] for row in generated}) == len(generated)
    assert all(row['rdkit_energy_record']['unit'] == 'kcal/mol' for row in generated)
    assert all(row['backend_energy_record'] is None for row in generated)
    assert not (Path(old.run_directory) / 'layout.json').exists()
    assert len(selected_artifacts(root)) == len(new.selected_paths)


def test_benchmark_consumes_manifest_and_detects_tampering(runs):
    _, new, _ = runs
    root = Path(new.run_directory)
    case = ConformerBenchmarkCase('butane', 'CCCC', reference_xyz=new.selected_paths[0])
    diagnosis = diagnose_conformer_case(case, new)
    assert diagnosis['reference_found_selected']
    assert diagnosis['selected_count'] == len(new.selected_paths)
    Path(new.selected_paths[0]).write_text('tampered')
    with pytest.raises(ValueError, match='checksum'):
        diagnose_conformer_case(case, new)


def test_v1_and_v2_remain_nonresumable_without_embedding(runs):
    old, new, options = runs
    for result in (old, new):
        with patch.object(conformer, '_embed_and_rank_conformers') as embedding:
            with pytest.raises(conformer.ConformerWorkflowError, match='already'):
                conformer.conformer_search('CCCC', root_directory=Path(result.run_directory).parent, **options)
            embedding.assert_not_called()


def test_native_refinement_names_and_energy_units_are_preserved(tmp_path, monkeypatch):
    pytest.importorskip('rdkit')
    names = []
    def optimize(molecule, parameters):
        names.append(molecule.name)
        molecule.energy = -1.25
        target = Path(f'job_{molecule.name}/result_{molecule.name}.xyz')
        target.parent.mkdir()
        molecule.mol_to_xyz(str(target))
        return True
    monkeypatch.setattr(conformer, 'optimise', optimize)
    result = conformer.conformer_search('CCCC', num_conformers=6, num_seeds=1,
        top_n=2, torsion_kicks=False, seed=13, num_threads=1,
        qc_params={'software':'xtb'}, root_directory=tmp_path)
    root = Path(result.run_directory)
    records = read_json(root / 'generated/records.json')['conformers']
    assert names and all(name.startswith('seed_000013_conf_') for name in names)
    selected = selected_artifacts(root)
    assert selected and all(ref['backend_energy_record']['unit'] == 'hartree' for ref in selected)
    assert all(ref['backend_energy_record']['value'] == -1.25 for ref in selected)
    assert all(ref['rdkit_energy_record']['unit'] == 'kcal/mol' for ref in selected)
    assert all((root / row['backend_path']).is_file() for row in records if row['backend_path'])
    assert all(row['backend_path'].startswith('refinement/execution/job_')
               for row in records if row['backend_path'])


def test_unknown_manifest_schema_fails_closed(runs):
    import json
    _, result, _ = runs
    root = Path(result.run_directory)
    path = root / 'results/manifest.json'
    manifest = read_json(path)
    manifest['manifest_schema_version'] = 999
    path.write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match='identities differ'):
        selected_artifacts(root)
