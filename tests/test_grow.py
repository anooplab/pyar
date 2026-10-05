"""Fixed-seed semantics, shared service and restart regression tests."""
import importlib
import json
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pytest

from pyar.core.molecule import Molecule
from pyar.growth.request import GrowRequest
from pyar.growth.service import select_geometry_only
from pyar.state.grow import GrowStateError
from pyar.workflows.grow import grow

workflow = importlib.import_module('pyar.workflows.grow')


def hydrogen(name='seed'):
    return Molecule(['H', 'H'], [[0, 0, 0], [0, 0, .74]], name=name)


def added(seed, monomer, step):
    return Molecule(seed.atoms_list + monomer.atoms_list,
                    np.vstack([seed.coordinates, monomer.coordinates + step * 3]),
                    name=f'stage_{step}', energy=None)


@pytest.mark.parametrize('count', [1, 3])
def test_sequential_stages(count, tmp_path):
    request = GrowRequest(hydrogen(), hydrogen('monomer'), count, maximum_number_of_seeds=2)
    calls = []
    def engine(aid, seeds, monomer, *args, **kwargs):
        calls.append(seeds)
        return [added(seeds[0], monomer, len(calls))]
    with patch.object(workflow, 'add_one', side_effect=engine):
        result = grow(request, output=tmp_path / 'grow')
    assert result.status == 'completed'
    assert len(calls) == count
    assert [len(pool[0]) for pool in calls] == [2 + 2 * i for i in range(count)]
    assert len(Molecule.from_xyz(result.selected_paths[0])) == 2 + 2 * count
    for step in range(count + 1):
        assert list((tmp_path / 'grow' / f'step_{step:03d}' / 'selected').glob('selected_*.xyz'))
    with patch.object(workflow, 'add_one') as engine:
        assert grow(request, output=tmp_path / 'grow').status == 'completed'
        engine.assert_not_called()


def test_interrupted_restart_reuses_pool(tmp_path):
    request = GrowRequest(hydrogen(), hydrogen('monomer'), 3)
    calls = []
    def interrupted(aid, seeds, monomer, *args, **kwargs):
        calls.append(aid)
        if len(calls) == 2:
            raise RuntimeError('interrupted')
        return [added(seeds[0], monomer, 1)]
    with patch.object(workflow, 'add_one', side_effect=interrupted):
        with pytest.raises(RuntimeError):
            grow(request, output=tmp_path / 'grow')
    before = (tmp_path / 'grow/step_001/selected/selected_000.xyz').read_bytes()
    with patch.object(workflow, 'add_one', side_effect=lambda aid, seeds, monomer, *a, **k: [added(seeds[0], monomer, 2)]) as engine:
        result = grow(request, output=tmp_path / 'grow')
        assert engine.call_count == 2
        assert len(engine.call_args_list[0].args[1][0]) == 4
    assert result.metadata['completed_additions'] == 3
    assert before == (tmp_path / 'grow/step_001/selected/selected_000.xyz').read_bytes()


@pytest.mark.parametrize('change', ['seed', 'monomer', 'count', 'selection_feature', 'backend_parameters', 'maximum_number_of_seeds'])
def test_incompatible_restart(change, tmp_path):
    options = dict(seed=hydrogen(), monomer=hydrogen('monomer'), count=1)
    with patch.object(workflow, 'add_one', return_value=[hydrogen()]):
        grow(GrowRequest(**options), output=tmp_path / 'grow')
    if change in {'seed', 'monomer'}:
        options[change].coordinates += 0.1
    else:
        options[change] = {'count': 2, 'selection_feature': 'soap', 'backend_parameters': {'nprocs': 2},
                           'maximum_number_of_seeds': 3}[change]
    with patch.object(workflow, 'add_one') as engine:
        with pytest.raises(GrowStateError, match='differs'):
            grow(GrowRequest(**options), output=tmp_path / 'grow')
        engine.assert_not_called()


@pytest.mark.parametrize("distance", ["euclidean", "graph-rmsd"])
def test_geometry_service_bounded_without_energy(tmp_path, monkeypatch, distance):
    monkeypatch.chdir(tmp_path)
    pool = [Molecule(['H', 'H'], [[0, 0, 0], [0, 0, 1 + i]], name=str(i)) for i in range(8)]
    result = select_geometry_only(pool, 3, feature='distance-histogram', distance=distance)
    assert len(result) <= 3
    assert all(m.energy is None for m in pool)
    assert json.loads(Path('selection_diagnostics.json').read_text())['energy_ranked'] is False


def test_geometry_end_to_end(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    request = GrowRequest(hydrogen(), Molecule(['He'], [[0, 0, 0]], name='monomer'), 3, 8,
                          maximum_number_of_seeds=2, connectivity_policy='off',
                          selection_feature='distance-histogram')
    result = grow(request)
    assert result.status == 'completed'
    assert all(stage['selected_count'] <= 2 for stage in result.metadata['stages'])
    assert all(stage['selection_mode'] == 'geometry-diversity' for stage in result.metadata['stages'])
    assert len(Molecule.from_xyz(result.selected_paths[0])) == 5
    assert Molecule.from_xyz(result.selected_paths[0]).energy is None


@pytest.mark.parametrize('kwargs', [{'count': 0}, {'number_of_orientations': 0},
                                    {'maximum_number_of_seeds': 0}, {'site': (-1, 0)},
                                    {'site': (0, 2)}, {'selection_distance': 'invalid'},
                                    {'connectivity_policy': 'invalid'},
                                    {'backend_parameters': {'software': 'unknown'}}])
def test_invalid_request(kwargs):
    options = dict(seed=hydrogen(), monomer=hydrogen(), count=1)
    options.update(kwargs)
    with pytest.raises(ValueError):
        GrowRequest(**options)


def test_invalid_spin():
    seed = hydrogen()
    seed.multiplicity = 2
    with pytest.raises(ValueError, match='multiplicity'):
        GrowRequest(seed, hydrogen(), 1)


def test_nonempty_output_rejected(tmp_path):
    (tmp_path / 'unrelated').write_text('preserve')
    with pytest.raises(GrowStateError, match='Non-empty'):
        grow(GrowRequest(hydrogen(), hydrogen(), 1), output=tmp_path)


def test_shared_engine_alias():
    from pyar.workflows import _growth
    from pyar.growth import service
    assert _growth is service
    aggregate = importlib.import_module('pyar.workflows.aggregate')
    assert aggregate.add_one is service.add_one


def test_cli(tmp_path, monkeypatch, capsys):
    from pyar.scripts.grow import main
    monkeypatch.chdir(tmp_path)
    hydrogen().mol_to_xyz('seed.xyz')
    hydrogen('monomer').mol_to_xyz('monomer.xyz')
    result = main(['seed.xyz', 'monomer.xyz', '--count', '1', '-N', '4',
                   '--maximum-number-of-seeds', '2', '--connectivity-policy', 'off'])
    assert result.status == 'completed'
    assert 'Growth completed' in capsys.readouterr().out


def test_cli_dispatch(monkeypatch):
    import pyar.cli
    monkeypatch.setattr('sys.argv', ['pyar-cli', 'grow', '--help'])
    with patch('pyar.scripts.grow.main') as main:
        pyar.cli.main()
    main.assert_called_once_with(['--help'])


def test_all_survivors_feed_next_stage(tmp_path):
    request = GrowRequest(hydrogen(), hydrogen('monomer'), 2, maximum_number_of_seeds=2)
    first_pool = [added(request.seed, request.monomer, 1), added(request.seed, request.monomer, 2)]
    with patch.object(workflow, 'add_one', side_effect=[first_pool, [added(first_pool[0], request.monomer, 3)]]) as engine:
        result = grow(request, output=tmp_path / 'grow')
    assert engine.call_args_list[1].args[1] is first_pool
    assert result.metadata['stages'][0]['selected_count'] == 2


def test_cumulative_electronic_states_reach_shared_engine(tmp_path):
    seed = Molecule(['H'], [[0, 0, 0]], name='seed', multiplicity=2)
    monomer = Molecule(['H'], [[0, 0, 0]], name='addend', charge=-1)
    request = GrowRequest(seed, monomer, 2, backend_parameters={'software': 'xtb', 'xtb_model': 'gfn2'})
    def engine(aid, seeds, monomer, orientations, qc, *args, **kwargs):
        result = added(seeds[0], monomer, 1)
        result.charge, result.multiplicity = qc['charge'], qc['multiplicity']
        return [result]
    with patch.object(workflow, 'add_one', side_effect=engine) as call:
        grow(request, output=tmp_path / 'grow')
    assert [item.args[4]['charge'] for item in call.call_args_list] == [-1, -2]
    assert all(item.args[4]['multiplicity'] == 2 for item in call.call_args_list)
    assert 'method' not in request.backend_parameters and 'basis' not in request.backend_parameters


def test_failed_comparison_conservatively_retains_geometry(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    from pyar.structure_comparison import GraphFirstDeduplicationComparator
    with patch.object(GraphFirstDeduplicationComparator, 'compare', side_effect=RuntimeError('mapping failed')):
        assert len(select_geometry_only([hydrogen('a'), hydrogen('b')], 3)) == 2


def test_snapshot_tampering_rejects_restart(tmp_path):
    request = GrowRequest(hydrogen(), hydrogen('monomer'), 1)
    with patch.object(workflow, 'add_one', return_value=[added(request.seed, request.monomer, 1)]):
        grow(request, output=tmp_path / 'grow')
    path = tmp_path / 'grow/step_001/selected/selected_000.xyz'
    path.write_text(path.read_text() + '\n')
    with pytest.raises(GrowStateError, match='modified'):
        grow(request, output=tmp_path / 'grow')


def test_geometry_only_aggregate_has_selected_paths(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    from pyar.workflows.aggregate import aggregate
    result = aggregate([hydrogen('A'), Molecule(['He'], [[0, 0, 0]], name='B')],
                       [1, 2], 4, {}, 2, 0, 1, None,
                       connectivity_policy='off', selection_feature='distance-histogram')
    assert result.workflow == 'aggregate'
    assert result.selected_paths
    structures = [Molecule.from_xyz(str(tmp_path / 'aggregates' / path)) for path in result.selected_paths]
    assert any(m.atoms_list.count('He') == 2 for m in structures)
    assert all(m.energy is None for m in structures)


def test_cli_invalid_site_no_output(tmp_path, monkeypatch):
    from pyar.scripts.grow import main
    monkeypatch.chdir(tmp_path)
    hydrogen().mol_to_xyz('a.xyz')
    hydrogen().mol_to_xyz('b.xyz')
    with pytest.raises(SystemExit):
        main(['a.xyz', 'b.xyz', '--count', '1', '--site', '-1', '0'])
    assert not Path('grow').exists()


def test_legacy_solvation_can_preserve_unbounded_geometry_generation(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    from pyar.growth import service
    pool = [hydrogen(str(i)) for i in range(4)]
    with patch.object(service, 'generate_orientations', return_value=iter(pool)):
        result = service.add_one('legacy', [hydrogen()], hydrogen(), 4, {}, 1, None,
                                 connectivity_policy='off', bound_geometry_only=False)
    assert result == pool


def test_solvation_explicitly_keeps_legacy_geometry_policy(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    from pyar.workflows import solvation
    with patch.object(solvation, 'add_one', return_value=[hydrogen()]) as engine:
        solvation.solvate([hydrogen()], hydrogen('monomer'), 1, 4, {}, 1)
    assert engine.call_args.kwargs['bound_geometry_only'] is False
