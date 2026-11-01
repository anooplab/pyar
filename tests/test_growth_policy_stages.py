"""Ownership and retry-budget contracts at the decomposed growth stages."""
from pathlib import Path

import pytest

from pyar.core.molecule import Molecule
from pyar.growth.policy import GrowthPolicy, GrowthRun, SeedRun
from pyar.growth import service


def test_policy_owns_settings_without_mutating_callers():
    original = {'software':'xtb', 'opt_threshold':'tight'}
    policy = GrowthPolicy(qc_params=original)
    original['opt_threshold'] = 'loose'
    assert policy.qc_params['opt_threshold'] == 'tight'
    with pytest.raises(TypeError):
        policy.qc_params['opt_threshold'] = 'normal'
    stages = service._optimization_stages(policy.qc_params)
    assert stages.loose['opt_threshold'] == 'loose'
    assert stages.refine['opt_threshold'] == 'normal'
    assert policy.qc_params['opt_threshold'] == 'tight'


def test_cycle_limited_candidates_keep_the_original_ten_round_budget(tmp_path, monkeypatch):
    molecule = Molecule(['H','H'], [[0,0,0],[0,0,.74]], name='candidate')
    policy = GrowthPolicy(qc_params={'software':'xtb'}, connectivity_policy='off')
    run = GrowthRun('addition', [molecule], molecule, policy, tmp_path, 'off',
                    service._optimization_stages(policy.qc_params))
    seed = SeedRun(run, 0, molecule)
    calls, diagnostics = [], []
    def optimize(pool, parameters, role, failure, *, directory):
        calls.append((tuple(pool), dict(parameters), role, Path(directory)))
        return ['CycleExceeded']
    monkeypatch.setattr(service, '_optimize_growth_batch', optimize)
    monkeypatch.setattr(service.clustering, 'remove_similar', lambda pool: pool)
    monkeypatch.setattr(service, '_write_growth_structure_analysis', lambda parents,pool,path: diagnostics.append(Path(path).name))
    assert service._retry_cycle_limited(seed, [molecule]) == []
    assert len(calls) == 10
    assert all(pool == (molecule,) and params['opt_threshold'] == 'loose'
               and role == 'growth_preselection' and directory == seed.directory
               for pool,params,role,directory in calls)
    assert diagnostics == [f'structural_analysis_round_{i:02d}.json' for i in range(1,11)]
