"""Adversarial evidence and bounded action tests for R3 recovery."""
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from pyar.persistence import read_json, write_json
from pyar.recovery import controller
from pyar.recovery.advanced import finite_difference_gradient
from pyar.recovery.models import (ExecutionObservation, Outcome, ProcessRecord,
                                  StepObservation)
from pyar.recovery.observations import snapshot_artifacts
from pyar.recovery.policy import RecoveryPolicy, flat_pes_evidence
from pyar.optimization_service import inspect_optimization


def test_policy_exposes_guided_separately_from_safe():
    assert RecoveryPolicy('off').mode == 'off'
    assert RecoveryPolicy('safe').mode == 'safe'
    assert RecoveryPolicy('guided').mode == 'guided'
    with pytest.raises(ValueError, match='off, safe or guided'):
        RecoveryPolicy('automatic')


def test_guided_stabilization_is_capability_gated():
    from pyar.backend_capabilities import optimization_recovery_capabilities
    assert 'guided_trust_hessian_restart' in optimization_recovery_capabilities('xtb', 'geometric')
    assert not optimization_recovery_capabilities('orca', 'native')


def test_guided_mode_fails_closed_for_unverified_orca_root(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    molecule = SimpleNamespace(name='h2', atoms_list=['H', 'H'],
        coordinates=[[0., 0., 0.], [0., 0., 1.]], charge=0, multiplicity=1,
        scftype='rhf', fragments=[])
    settings = dict(software='orca', geometry_optimizer='native', opt_target='minimum',
        opt_cycles=5, opt_threshold='normal', scf_cycles=100, method='BP86',
        basis='def2-SVP', nprocs=1)
    with pytest.raises(ValueError, match='recovery_unsupported'):
        controller.preflight_recovery([molecule], settings, RecoveryPolicy('guided'))
    assert not Path('job_h2').exists()


def test_safe_geometric_passes_only_bounded_trust_controls(monkeypatch):
    import json
    import sys
    from pyar.backends import safe_geometric
    settings = dict(software='xtb', xtb_model='gfn2', geometry_optimizer='geometric',
        opt_target='minimum', opt_cycles=6, opt_threshold='normal', nprocs=1,
        charge=0, multiplicity=1, scftype='rhf', xtb_unpaired_electrons=0,
        fragments=[], _ore_safe=True, _ore_evaluation_limit=20, _ore_initial_trust=.05)
    monkeypatch.setattr(sys, 'argv', ['safe_geometric.py', 'trial_h2.xyz', json.dumps(
        {'qc_params': settings, 'fragment_indices': [], 'opt_target': 'minimum'})])
    captured = {}
    monkeypatch.setattr(safe_geometric.driver, 'run_optimizer',
        lambda **kwargs: captured.update(kwargs))
    optimizer_class = safe_geometric.driver.Optimizer
    try:
        assert safe_geometric.main() == 0
    finally:
        safe_geometric.driver.Optimizer = optimizer_class
    assert captured['trust'] == .05
    assert captured['tmax'] == pytest.approx(.15)
    assert captured['maxiter'] == 6


def test_flat_pes_requires_multiple_accepted_steps_above_requested_threshold():
    steps = tuple(StepObservation(i, -1.0 + i * 1e-8, .004 + i * 1e-5,
                                  .006 + i * 1e-5, accepted=True) for i in range(4))
    assert flat_pes_evidence(steps, 3e-4)
    assert not flat_pes_evidence(steps[:3], 3e-4)
    assert not flat_pes_evidence(steps, .005)
    assert not flat_pes_evidence(tuple(StepObservation(
        s.index, s.effective_energy_hartree + i * 1e-4,
        s.gradient_rms_hartree_per_bohr, s.gradient_max_hartree_per_bohr, accepted=True)
        for i, s in enumerate(steps)), 3e-4)


def test_finite_difference_gradient_matches_quadratic_and_is_bounded():
    coordinates = np.array([[.2, -.1, .4], [-.3, .5, .1]])
    diagonal = np.arange(1, 7, dtype=float).reshape((2, 3))
    gradient = diagonal * coordinates
    def energy(value):
        return float(.5 * np.sum(diagonal * value**2))
    result = finite_difference_gradient(energy, coordinates, gradient)
    assert result.evaluations == 12
    assert result.consistent
    assert result.max_absolute_error_hartree_per_bohr < 1e-8

    noisy = gradient.copy()
    noisy[0, 0] += .5
    assert not finite_difference_gradient(energy, coordinates, noisy).consistent
    calls = []
    with pytest.raises(ValueError, match='cap is 11'):
        finite_difference_gradient(lambda value: calls.append(value), coordinates, gradient,
                                   max_evaluations=11)
    assert calls == []


def _write_stagnant_attempt(directory):
    intent = read_json(directory / 'input.json')
    request = controller._request(intent['request'])
    root = directory / ('job_' + request.name)
    root.mkdir()
    (directory / 'worker.out').write_text('synthetic native execution\n')
    (root / ('trial_' + request.name + '.xyz')).write_text(
        f'2\nenergy: -1.0\nH {request.coordinates_angstrom[0][0]} '
        f'{request.coordinates_angstrom[0][1]} {request.coordinates_angstrom[0][2]}\n'
        f'H {request.coordinates_angstrom[1][0]} {request.coordinates_angstrom[1][1]} '
        f'{request.coordinates_angstrom[1][2]}\n')
    files, lines = [], ['geomeTRIC started. Version: 1.1.1']
    for index in range(4):
        coordinates = [list(row) for row in request.coordinates_angstrom]
        coordinates[-1][-1] -= index * .001
        gradient = .004 + index * .00001
        energy = -1.0 + index * 1e-8
        filename = f'ore_accepted_{index:06d}.json'
        write_json(root / filename, dict(schema_version=1, version='1.1.1', step=index,
            objective_fingerprint=request.objective_fingerprint, symbols=['H', 'H'],
            coordinates_angstrom=coordinates, energy_hartree=energy,
            gradient_hartree_per_bohr=[[0, 0, gradient], [0, 0, -gradient]],
            gradient_rms=gradient, gradient_max=gradient))
        files.append(filename)
        lines.append(f'Step {index} : Grad = {gradient}/{gradient} E (change) = {energy}')
    write_json(root / 'ore_accepted.json', dict(schema_version=1, files=files))
    lines.append('Maximum iterations reached')
    (root / 'geometric.out').write_text('\n'.join(lines) + '\n')
    write_json(root / 'ore_evaluations.json', {'used': 4})
    for number in range(1, 5):
        (root / f'ore_backend_{number:06d}.out').write_text(
            'xtb version 6.7.1\nHamiltonian GFN2-xTB\nnormal termination of xtb\n')
    execution = ExecutionObservation(request.request_fingerprint, snapshot_artifacts(root),
        (ProcessRecord('python', 1),), True, 'CycleExceeded', elapsed_seconds=.1)
    write_json(directory / 'result.json', inspect_optimization(
        request, root, execution=execution).to_dict())
    return 0


def test_guided_plateau_rolls_back_and_rebuilds_geometric_model(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    molecule = SimpleNamespace(name='h2', atoms_list=['H', 'H'],
        coordinates=[[0., 0., 0.], [0., 0., 1.]], charge=0, multiplicity=1,
        scftype='rhf', fragments=[])
    settings = dict(software='xtb', xtb_model='gfn2', geometry_optimizer='geometric',
        opt_target='minimum', opt_cycles=3, opt_threshold='normal', nprocs=1)
    calls = []
    def execute(directory, timeout):
        calls.append(directory)
        if len(calls) == 1:
            return _write_stagnant_attempt(directory)
        intent = read_json(directory / 'input.json')
        assert intent['request']['settings']['_ore_initial_trust'] == .05
        assert intent['request']['settings']['_ore_recovery_action'] == 'trust_hessian_restart'
        return _complete_fake_attempt(directory)
    monkeypatch.setattr(controller, '_execute_attempt', execute)
    report = controller.run_recovery(molecule, settings, RecoveryPolicy('guided'))
    assert report.result.status is Outcome.RECOVERED_CONVERGED
    assert report.reserved_cycles == 9
    assert report.backend_evaluations == 6
    assert [attempt.action for attempt in report.result.attempts] == [
        'initial', 'trust_hessian_restart']
    assert report.result.attempts[0].request.coordinates_angstrom != \
        report.result.attempts[1].request.coordinates_angstrom
    assert report.result.attempts[0].request.objective_fingerprint == \
        report.result.attempts[1].request.objective_fingerprint
    assert report.result.attempts[0].original_scientific_fingerprint == \
        report.result.attempts[1].original_scientific_fingerprint
    assert report.result.attempts[1].request.settings['opt_threshold'] == 'normal'
    assert read_json(Path('job_h2/recovery/attempts/attempt_002/input.json'))[
        'recovery_action_policy_version'] == 2


def _complete_fake_attempt(directory):
    # Reuse the established adversarial fixture builder for a converged second run.
    intent = read_json(directory / 'input.json')
    request = controller._request(intent['request'])
    root = directory / ('job_' + request.name)
    root.mkdir()
    (directory / 'worker.out').write_text('synthetic native execution\n')
    start = request.coordinates_angstrom
    (root / ('trial_' + request.name + '.xyz')).write_text(
        f'2\nenergy: -1.0\nH {start[0][0]} {start[0][1]} {start[0][2]}\n'
        f'H {start[1][0]} {start[1][1]} {start[1][2]}\n')
    end = [list(row) for row in start]
    end[-1][-1] -= .002
    gradient, energy = .00001, -1.1
    manifest = []
    for index, coordinates, e, g in [(0, start, -1., .01), (1, end, energy, gradient)]:
        filename = f'ore_accepted_{index:06d}.json'
        write_json(root / filename, dict(schema_version=1, version='1.1.1', step=index,
            objective_fingerprint=request.objective_fingerprint, symbols=['H', 'H'],
            coordinates_angstrom=coordinates, energy_hartree=e,
            gradient_hartree_per_bohr=[[0, 0, g], [0, 0, -g]], gradient_rms=g, gradient_max=g))
        manifest.append(filename)
    write_json(root / 'ore_accepted.json', dict(schema_version=1, files=manifest))
    (root / 'geometric.out').write_text(
        'geomeTRIC started. Version: 1.1.1\nStep 0 : Gradient = .01/.01 Energy = -1.0\n'
        f'Step 1 : Grad = {gradient}/{gradient} E (change) = {energy}\nConverged!\n')
    forces = [[0, 0, -gradient * 27.211386245988 / .529177210903],
              [0, 0, gradient * 27.211386245988 / .529177210903]]
    write_json(root / 'pyar_geometric_state.json', dict(schema_version=1, software='xtb',
        opt_target='minimum', positions_angstrom=end, backend_energy_hartree=energy,
        total_energy_hartree=energy, forces_ev_per_angstrom=forces))
    (root / ('trial_' + request.name + '_optim.xyz')).write_text(
        f'2\nenergy: {energy}\nH {end[0][0]} {end[0][1]} {end[0][2]}\n'
        f'H {end[1][0]} {end[1][1]} {end[1][2]}\n')
    write_json(root / 'ore_evaluations.json', {'used': 2})
    for number in (1, 2):
        (root / f'ore_backend_{number:06d}.out').write_text(
            'xtb version 6.7.1\nHamiltonian GFN2-xTB\nnormal termination of xtb\n')
    execution = ExecutionObservation(request.request_fingerprint, snapshot_artifacts(root),
        (ProcessRecord('python', 0),), True, True, elapsed_seconds=.1)
    write_json(directory / 'result.json', inspect_optimization(
        request, root, execution=execution).to_dict())
    return 0
