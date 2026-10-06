"""Modern wrappers preserve canonical generation workflows and read-only checks."""
import importlib
import json
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace
from unittest.mock import patch

import pytest

from pyar import modern_cli
from pyar.scripts import modern_aggregate as aggregate, modern_grow as grow, modern_conformer as conformer
from pyar.state.aggregate import AggregateRunState
from pyar.workflows import conformer as conformer_workflow


@pytest.fixture
def inputs(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    for name, atoms in [('a', ['He']), ('b', ['He']), ('odd', ['H'])]:
        Path(name + '.xyz').write_text(f'{len(atoms)}\n{name}\n' + ''.join(f'{a} 0 0 0\n' for a in atoms))
    return ['a.xyz', 'b.xyz']


def result():
    return SimpleNamespace(status='completed', selected_paths=['selected.xyz'], run_directory='run',
                           metadata={'completed_additions': 4})


@pytest.mark.parametrize('command', ['conformer', 'aggregate', 'grow'])
def test_help(command, capsys):
    with pytest.raises(SystemExit) as exc:
        modern_cli.main([command, '--help'])
    assert exc.value.code == 0
    output = capsys.readouterr().out
    assert '--backend' in output and '--check' in output
    assert '--software' not in output


def test_top_level_help_remains_lazy():
    code = """import sys
from pyar.modern_cli import main
try: main(['--help'])
except SystemExit: pass
assert not any(m in sys.modules for m in ['rdkit', 'ase', 'geometric', 'torch', 'pyar.growth.service'])
"""
    subprocess.run([sys.executable, '-c', code], check=True, capture_output=True)


def test_minimal_aggregate(inputs):
    with patch.object(aggregate, 'aggregate', return_value=result()) as engine, patch.object(aggregate, 'preflight') as preflight:
        modern_cli.main(['aggregate', *inputs])
    args = engine.call_args.args
    assert args[1:7] == ([1, 1], 8, {}, 8, 0, 10)
    assert [m.charge for m in args[0]] == [0, 0]
    preflight.assert_not_called()


@pytest.mark.parametrize('arguments, sizes', [(['a.xyz', '--size', '6'], [6]),
                                             (['C', 'H', '--size', '1', '4'], [1, 4]),
                                             (['--formula', 'CH4'], [1, 4])])
def test_aggregate_compositions(inputs, arguments, sizes):
    with patch.object(aggregate, 'aggregate', return_value=result()) as engine:
        modern_cli.main(['aggregate', *arguments])
    assert engine.call_args.args[1] == sizes


@pytest.mark.parametrize('arguments', [['a.xyz'], ['a.xyz', 'b.xyz', '--size', '2'],
    ['a.xyz', '--size', '0'], ['a.xyz', '--size', '1'], ['--formula', 'CH4', 'a.xyz'],
    ['a.xyz', 'b.xyz', '--orientations', '0'], ['a.xyz', 'b.xyz', '--selection-feature', 'invalid'],
    ['a.xyz', 'b.xyz', '--site', '-1', '0']])
def test_aggregate_invalid(inputs, arguments):
    with patch.object(aggregate, 'aggregate') as engine, pytest.raises(SystemExit) as exc:
        modern_cli.main(['aggregate', *arguments])
    assert exc.value.code != 0
    engine.assert_not_called()
    assert not Path('aggregates').exists()


def test_selection_passthrough(inputs):
    options = ['--maximum-number-of-seeds', '3', '--number-of-pathways', '2', '--first-pathway', '1',
               '--connectivity-policy', 'off', '--selection-feature', 'distance-histogram',
               '--selection-algorithm', 'agglomerative', '--selection-distance', 'graph-rmsd',
               '--selection-system-type', 'molecular']
    with patch.object(aggregate, 'aggregate', return_value=result()) as engine:
        modern_cli.main(['aggregate', *inputs, *options])
    assert engine.call_args.args[4:7] == (3, 1, 2)
    assert engine.call_args.args[8:] == ('off', 'distance-histogram', 'agglomerative', 'graph-rmsd', 'molecular')


def test_minimal_grow(inputs):
    with patch.object(grow, 'grow', return_value=result()) as engine, patch.object(grow, 'preflight') as preflight:
        modern_cli.main(['grow', *inputs, '--count', '4'])
    request = engine.call_args.args[0]
    assert (request.count, request.number_of_orientations, request.maximum_number_of_seeds) == (4, 8, 12)
    assert request.backend_parameters == {}
    assert engine.call_args.kwargs == {'output': 'grow'}
    preflight.assert_not_called()


@pytest.mark.parametrize('charges, states', [(['-1'], [-1, -1]), (['0', '-1'], [0, -1])])
def test_grow_states_and_output(inputs, charges, states):
    with patch.object(grow, 'grow', return_value=result()) as engine:
        modern_cli.main(['grow', *inputs, '--count', '4', '--charge', *charges, '--output', 'custom', '--site', '0', '0'])
    request = engine.call_args.args[0]
    assert [request.seed.charge, request.monomer.charge] == states
    assert [request.seed.multiplicity, request.monomer.multiplicity] == [2 if c else 1 for c in states]
    assert request.site == (0, 0)
    assert engine.call_args.kwargs['output'] == 'custom'


@pytest.mark.parametrize('options', [['--count', '0'], ['--count', '-1'], ['--count', '1', '--site', '0', '1'],
    ['--count', '1', '--multiplicity', '2'], ['--count', '1', '--charge', '0', '0', '0']])
def test_grow_invalid(inputs, options):
    with patch.object(grow, 'grow') as engine, pytest.raises(SystemExit) as exc:
        modern_cli.main(['grow', *inputs, *options])
    assert exc.value.code != 0
    engine.assert_not_called()
    assert not Path('grow').exists()


@pytest.mark.parametrize('command', ['aggregate', 'grow', 'conformer'])
def test_backend_policy(inputs, command):
    arguments = {'aggregate': inputs, 'grow': [*inputs, '--count', '4'], 'conformer': ['CCO']}[command]
    module, target = {'aggregate': (aggregate, 'aggregate'), 'grow': (grow, 'grow'),
                      'conformer': (conformer_workflow, 'conformer_search')}[command]
    with patch('pyar.cli.shutil.which', return_value='/bin/xtb'), patch.object(module, target, return_value=result()) as engine:
        modern_cli.main([command, *arguments, '--backend', 'xtb'])
    qc = (engine.call_args.args[3] if command == 'aggregate' else engine.call_args.args[0].backend_parameters
          if command == 'grow' else engine.call_args.kwargs['qc_params'])
    assert qc['software'] == 'xtb' and qc['xtb_model'] == 'gfn2'
    assert qc['method'] is None and qc['basis'] is None


@pytest.mark.parametrize('command', ['aggregate', 'grow', 'conformer'])
def test_check_no_side_effects(inputs, command):
    arguments = {'aggregate': inputs, 'grow': [*inputs, '--count', '4'], 'conformer': ['CCO']}[command]
    module, target = {'aggregate': (aggregate, 'aggregate'), 'grow': (grow, 'grow'),
                      'conformer': (conformer_workflow, 'conformer_search')}[command]
    before = set(Path('.').iterdir())
    with patch.object(module, target) as engine:
        modern_cli.main([command, *arguments, '--check'])
    engine.assert_not_called()
    assert set(Path('.').iterdir()) == before


def test_grow_check_exits_before_workflow(inputs):
    before = set(Path('.').iterdir())
    with patch.object(grow, 'grow') as engine:
        result = modern_cli.main(['grow', *inputs, '--count', '4', '--check'])
    assert result is None
    engine.assert_not_called()
    assert set(Path('.').iterdir()) == before


@pytest.mark.parametrize('command', ['aggregate', 'grow', 'conformer'])
def test_backend_missing_fails_before_run(inputs, command, capsys):
    arguments = {'aggregate': inputs, 'grow': [*inputs, '--count', '4'], 'conformer': ['CCO']}[command]
    with patch('pyar.cli.shutil.which', return_value=None), pytest.raises(SystemExit) as exc:
        modern_cli.main([command, *arguments, '--backend', 'xtb', '--check'])
    assert exc.value.code != 0
    assert 'xtb' in capsys.readouterr().err
    assert not any(Path(d).exists() for d in ('aggregates', 'grow', 'conformers'))


def test_conformer_defaults(inputs):
    with patch.object(conformer_workflow, 'conformer_search', return_value=result()) as engine:
        modern_cli.main(['conformer', 'CCO'])
    options = engine.call_args.kwargs
    expected = dict(num_conformers=150, top_n=10, num_seeds=5, diversity_fraction=.2, compactness_fraction=.2,
                    rms_threshold=.25, use_random_coords=True, torsion_kicks=True, torsion_rounds=2,
                    torsion_kicks_per_conformer=6, torsion_max_bonds=3, torsion_dedup_rms=.5,
                    dedup_atom_mode='heavy', force_field='auto', seed=1, max_iterations=200)
    assert {key: options[key] for key in expected} == expected
    assert options['qc_params'] is None


def test_conformer_advanced_and_charged_input(inputs):
    with patch.object(conformer_workflow, 'conformer_search', return_value=result()) as engine:
        modern_cli.main(['conformer', '[NH4+]', '--num-conformers', '20', '--num-seeds', '2', '--top-n', '3',
                         '--no-use-random-coords', '--no-torsion-kicks', '--torsion-rounds', '0', '--seed', '9'])
    options = engine.call_args.kwargs
    assert options['charge'] == 1 and options['multiplicity'] == 1
    assert (options['num_conformers'], options['num_seeds'], options['top_n'], options['seed']) == (20, 2, 3, 9)
    assert not options['torsion_kicks'] and not options['use_random_coords']


def test_conformer_odd_spin_backend(inputs):
    with patch('pyar.cli.shutil.which', return_value='/bin/xtb'), patch.object(conformer_workflow, 'conformer_search', return_value=result()) as engine:
        modern_cli.main(['conformer', '[CH3]', '--backend', 'xtb'])
    options = engine.call_args.kwargs
    assert options['multiplicity'] == 2 and options['scftype'] == 'uhf'
    assert options['qc_params']['xtb_unpaired_electrons'] == 1
    with pytest.raises(SystemExit):
        modern_cli.main(['conformer', '[CH3]', '--multiplicity', '1', '--check'])


def test_missing_rdkit(inputs, capsys):
    with patch.object(conformer_workflow, '_rdkit_modules', side_effect=ImportError('RDKit required; pip install "pyar-chem[conformer]"')), pytest.raises(SystemExit):
        modern_cli.main(['conformer', 'CCO', '--check'])
    assert 'pyar-chem[conformer]' in capsys.readouterr().err
    assert not Path('conformers').exists()


def test_missing_obabel(inputs, capsys):
    with patch('pyar.scripts.modern_conformer.shutil.which', return_value=None), pytest.raises(SystemExit):
        modern_cli.main(['conformer', 'a.xyz', '--check'])
    assert 'obabel' in capsys.readouterr().err
    assert not Path('conformers').exists()


def test_aggregate_restart_read_only(inputs):
    args = aggregate.build_parser().parse_args(inputs)
    request = aggregate.resolve_request(args)
    Path('aggregates').mkdir()
    AggregateRunState.create(Path.cwd(), request.to_state_dict(), ['ab'])
    before = Path('aggregates/state.json').read_bytes()
    modern_cli.main(['aggregate', *inputs, '--check'])
    assert Path('aggregates/state.json').read_bytes() == before
    with pytest.raises(SystemExit):
        modern_cli.main(['aggregate', *inputs, '--size', '2', '1', '--check'])
    assert Path('aggregates/state.json').read_bytes() == before


def test_unsafe_outputs(inputs):
    for name, arguments in [('grow', ['grow', *inputs, '--count', '1']),
                            ('aggregates', ['aggregate', *inputs])]:
        Path(name).mkdir()
        Path(name, 'unrelated').write_text('preserve')
        with pytest.raises(SystemExit):
            modern_cli.main([*arguments, '--check'])
        assert Path(name, 'unrelated').read_text() == 'preserve'


def test_grow_restart_read_only(inputs):
    # Use the canonical workflow to create genuine snapshots and state.
    modern_cli.main(['grow', *inputs, '--count', '1', '--orientations', '1'])
    before = {p: p.read_bytes() for p in Path('grow').rglob('*') if p.is_file()}
    modern_cli.main(['grow', *inputs, '--count', '1', '--orientations', '1', '--check'])
    with pytest.raises(SystemExit):
        modern_cli.main(['grow', *inputs, '--count', '2', '--orientations', '1', '--check'])
    assert before == {p: p.read_bytes() for p in before}


def test_legacy_conformer_defaults_unchanged():
    from pyar.scripts.conformer import argument_parse
    args = argument_parse(['CCO'])
    assert args.multiplicity == 1 and args.scftype == 'rhf'
    assert args.method == 'BP86' and args.basis == 'def2-SVP'
    assert args.software is None


@pytest.mark.parametrize('command', ['grow', 'aggregate', 'conformer'])
def test_backend_check_calls_same_preflight(inputs, command):
    arguments = {'aggregate': inputs, 'grow': [*inputs, '--count', '4'], 'conformer': ['CCO']}[command]
    module = {'aggregate': aggregate, 'grow': grow, 'conformer': conformer}[command]
    with patch.object(module, 'preflight', return_value=['xtb']) as preflight:
        modern_cli.main([command, *arguments, '--backend', 'xtb', '--check'])
    assert preflight.call_args.args[0]['xtb_model'] == 'gfn2'
    assert len(preflight.call_args.args[1]) == (1 if command == 'conformer' else 2)
    assert not any(Path(d).exists() for d in ('grow', 'aggregates', 'conformers'))


def test_conformer_check_never_embeds(inputs):
    from rdkit.Chem import AllChem
    with patch.object(AllChem, 'EmbedMultipleConfs') as embed:
        modern_cli.main(['conformer', 'CCO', '--check'])
    embed.assert_not_called()
    Path('conformers').mkdir()
    Path('conformers/state.json').write_text('{"status":"completed"}')
    before = Path('conformers/state.json').read_bytes()
    with pytest.raises(SystemExit):
        modern_cli.main(['conformer', 'CCO', '--check'])
    assert Path('conformers/state.json').read_bytes() == before


@pytest.mark.parametrize('options', [['--backend', 'xtb', '--basis', 'def2-SVP'],
    ['--backend', 'unknown'], ['--method', 'BP86'], ['--backend', 'xtb', '--xtb-model', 'unknown'],
    ['--backend', 'orca', '--method', 'g-xTB'],
    ['--backend', 'orca', '--method', 'g-xTB', '--gxtb-wrapper', 'missing'],
    ['--backend', 'orca', '--gxtb-wrapper', 'unused']])
def test_backend_settings_rejected(inputs, options):
    with patch.object(grow, 'grow') as engine, pytest.raises(SystemExit):
        modern_cli.main(['grow', *inputs, '--count', '1', *options])
    engine.assert_not_called()
    assert not Path('grow').exists()


def test_missing_input_all_validated(inputs):
    with patch.object(aggregate, 'aggregate') as engine, pytest.raises(SystemExit):
        modern_cli.main(['aggregate', *inputs, 'missing.xyz'])
    engine.assert_not_called()
    assert not Path('aggregates').exists()


def test_formula_geometry_stable_for_restart_resolution():
    first = aggregate.build_parser().parse_args(['C', 'H', '--size', '1', '4'])
    second = aggregate.build_parser().parse_args(['C', 'H', '--size', '1', '4'])
    assert aggregate.resolve_request(first).to_state_dict() == aggregate.resolve_request(second).to_state_dict()


def test_xyz_conformer_check_validates_in_memory_without_writes(inputs):
    from rdkit import Chem
    molecule = Chem.AddHs(Chem.MolFromSmiles('O'))
    conformer = Chem.Conformer(molecule.GetNumAtoms())
    for index in range(molecule.GetNumAtoms()):
        conformer.SetAtomPosition(index, (float(index), 0., 0.))
    molecule.AddConformer(conformer)
    Path('water.xyz').write_text('3\nwater\nO 0 0 0\nH .8 0 .5\nH -.8 0 .5\n')
    before = set(Path('.').iterdir())
    with patch.object(conformer_workflow, '_load_xyz_with_openbabel', return_value=molecule) as convert:
        modern_cli.main(['conformer', 'water.xyz', '--check'])
    assert convert.call_args.args[1] is None
    assert set(Path('.').iterdir()) == before
