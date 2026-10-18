"""Diagnosis at workflow-owned constrained scan and path boundaries."""
import logging
from pathlib import Path

import numpy as np

from pyar.persistence import write_json
from .workflow import stage_diagnostic


def persist_stage_diagnostic(directory, value, *, filename='optimization_diagnostic.json'):
    """A diagnostic sidecar is never allowed to replace the owner's outcome."""
    try:
        write_json(Path(directory) / filename, value)
    except Exception as exc:
        logging.getLogger(__name__).warning('Could not persist stage diagnostics: %s', exc)


def scan_diagnostic(molecule, request, result=None, *, exception=None):
    gates = dict(optimizer=None, constraints=None, atom_mapping=None, physical_energy=None)
    if result is not None:
        gates['optimizer'] = bool(result.success)
        if result.success:
            try:
                frames, profile = result.frames, result.profile
                count = len(frames) == len(profile) == request.n_points
                grid = np.linspace(request.start_distance_angstrom, request.end_distance_angstrom,
                                   request.n_points)
                gates['atom_mapping'] = count and all(
                    list(symbols) == list(molecule.atoms_list) and
                    np.asarray(xyz).shape == (len(symbols), 3) and np.isfinite(xyz).all()
                    for symbols, xyz, _ in frames)
                gates['constraints'] = bool(count and gates['atom_mapping'] and all(
                    np.isclose(np.linalg.norm(np.asarray(frame[1])[request.atom_i] -
                                              np.asarray(frame[1])[request.atom_j]), target,
                               atol=1e-3, rtol=0) and
                    np.isclose(point['target_distance_angstrom'], target, atol=1e-3, rtol=0)
                    for frame, point, target in zip(frames, profile, grid)))
                gates['physical_energy'] = bool(count and all(
                    np.isfinite(float(point['energy_hartree'])) for point in profile))
            except (ValueError, TypeError, KeyError, IndexError):
                gates.update(constraints=False, atom_mapping=False, physical_energy=False)
    reason = None
    if gates['constraints'] is False:
        reason = 'CONSTRAINT_VIOLATION'
    elif gates['atom_mapping'] is False or gates['physical_energy'] is False:
        reason = 'INVALID_RESULT'
    return stage_diagnostic('bond_scan_constrained', gates, exception=exception, reason=reason)


def run_scan_once(runner, molecule, request, directory, settings):
    """Observe native/generic scan execution; preserve original exception types."""
    try:
        if settings.get('workflow_recovery') is not None:
            raise ValueError('recovery_unsupported: constrained scans have no certified recovery checkpoint')
        result = runner(molecule, request, directory, settings)
    except BaseException as exc:
        # Cancellation must propagate; best-effort evidence cannot suppress it.
        persist_stage_diagnostic(directory, scan_diagnostic(molecule, request, exception=exc))
        raise
    diagnostic = scan_diagnostic(molecule, request, result)
    persist_stage_diagnostic(directory, diagnostic)
    # Reject a claimed success that violates the scan's unchanged objective.
    if result.success and diagnostic['stage_gate'] != 'passed':
        result.success = False
        result.status = 'scan_validation_failed'
    return result


def preserve_path_observations(output, label, symbols, frames, energies):
    """Preserve finite observations on failure without declaring accepted steps.

    geomeTRIC progress can include an unaccepted final trial; this file is
    explicitly excluded from checkpoint/restart authority.
    """
    try:
        from pyar.neb import _write_xyz_trajectory
        valid = [(np.asarray(frame), float(energy)) for frame, energy in zip(frames, energies)
                 if np.asarray(frame).shape == (len(symbols), 3)
                 and np.isfinite(frame).all() and np.isfinite(float(energy))]
        if valid:
            _write_xyz_trajectory(Path(output) / f'{label}_failed_observed.xyz', symbols,
                                  [p[0] for p in valid], [p[1] for p in valid])
    except Exception as exc:
        logging.getLogger(__name__).warning('Could not preserve failed path observations: %s', exc)


def preserve_band_observations(output, symbols, chain):
    """Even a partially initialized band must not replace the original error."""
    try:
        preserve_path_observations(output, 'neb', symbols,
            [structure.M.xyzs[0] for structure in chain.Structures],
            [structure.energy for structure in chain.Structures])
    except Exception as exc:
        logging.getLogger(__name__).warning('Could not inspect failed band observations: %s', exc)


def preserve_optimizer_observations(output, label, symbols, optimizer):
    try:
        progress = optimizer.progress
        preserve_path_observations(output, label, symbols, progress.xyzs, progress.qm_energies)
    except Exception as exc:
        logging.getLogger(__name__).warning('Could not inspect failed optimizer observations: %s', exc)
