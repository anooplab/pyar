"""Physical energies must survive optimization, caching, and reaction restart."""

import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from pyar import optimiser
from pyar.backends import write_xyz
from pyar.core.molecule import Molecule
from pyar.selection.deduplication import remove_similar
from pyar.state.reaction import ReactionRunState, ReactionStateError
from pyar.workflows.reaction import build_reaction_request


def test_selection_and_cached_results_use_physical_energy(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr('pyar.backends.geometric._find_geometric_executable',
                        lambda: 'geometric-optimize')
    calls = []

    def fake_run(command, **kwargs):
        calls.append(command)
        name = Path.cwd().name.removeprefix('job_')
        physical, total = (-2., 100.) if name == 'a' else (-1., -100.)
        coordinates = [[0., 0., 0.], [2., 0., 0.]]
        write_xyz(['C', 'C'], coordinates, f'trial_{name}_optim.xyz')
        Path('pyar_geometric_state.json').write_text(json.dumps({
            'energy_hartree': total, 'total_energy_hartree': total,
            'backend_energy_hartree': physical, 'bias_energy_hartree': total-physical,
            'positions_angstrom': coordinates,
        }))
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr('pyar.backends.geometric.subp.run', fake_run)
    parameters = dict(software='xtb', geometry_optimizer='geometric', gamma=100.,
                      bias_controller='adaptive')
    molecules = [Molecule(['C', 'C'], np.array([[0., 0., 0.], [2., 0., 0.]]),
                          fragments=[[0], [1]], name=name) for name in ('a', 'b')]
    for molecule in molecules:
        assert optimiser.optimise(molecule, parameters) is True
        saved = Molecule.from_xyz(f'job_{molecule.name}/result_{molecule.name}.xyz')
        assert saved.energy == molecule.energy
    assert [m.name for m in remove_similar(molecules)] == ['a']
    # A legacy XYZ label must not override physical energy in calculator state.
    write_xyz(['C', 'C'], molecules[0].coordinates, 'job_a/result_a.xyz', energy=100.)
    assert optimiser.optimise(molecules[0], parameters) is True
    assert molecules[0].energy == -2.
    assert len(calls) == 2
    state = json.loads(Path('job_a/pyar_geometric_state.json').read_text())
    del state['backend_energy_hartree']
    Path('job_a/pyar_geometric_state.json').write_text(json.dumps(state))
    assert optimiser.optimise(molecules[0], parameters) is True
    assert molecules[0].energy == -2.
    assert len(calls) == 3
    state = json.loads(Path('job_a/pyar_geometric_state.json').read_text())
    state['backend_energy_hartree'] = -99.
    state['positions_angstrom'][1][0] = 3.
    Path('job_a/pyar_geometric_state.json').write_text(json.dumps(state))
    assert optimiser.optimise(molecules[0], parameters) is True
    assert molecules[0].energy == -2.
    assert len(calls) == 4


def test_reaction_restart_preserves_physical_energy_and_rejects_old_convention(tmp_path):
    a = Molecule(['C'], np.zeros((1, 3)), name='a')
    b = Molecule(['C'], np.ones((1, 3)), name='b')
    candidate = a + b
    candidate.energy = -2.
    parameters = dict(software='xtb', geometry_optimizer='geometric')
    request = build_reaction_request(a, b, [100., 200.], 1, parameters, None, 2.3)
    state = ReactionRunState.create(tmp_path, request, [candidate], (a, b))
    state.record_job('candidate', 100., True, [], [candidate])
    loaded = ReactionRunState.load(tmp_path, request)
    assert loaded.current_survivor_molecules()[0].energy == -2.
    del state.data['request']['backend_parameters']['reaction_energy_convention']
    state.save()
    with pytest.raises(ReactionStateError, match='does not match'):
        ReactionRunState.load(tmp_path, request)
    with pytest.raises(ReactionStateError, match='do not identify physical energies'):
        ReactionRunState.migrate_legacy(tmp_path, {'0100': [candidate]}, request)
