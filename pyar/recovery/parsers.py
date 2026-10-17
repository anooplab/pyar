"""Read-only pilot parsers: ORCA 6.1, xTB 6.7, geomeTRIC 1.1.

Only documented success/limit/error markers are classified. Unknown versions or
incomplete evidence cannot certify convergence. No remedy is selected here.
"""
from dataclasses import dataclass, replace
import math
import re

from pyar.persistence import read_json
from .models import FailureCode, FailureDiagnostic, GeometryFrame, StepObservation, matrix
from .validation import NUMBER, close_matrix, finite_number, orca_coordinates, read_xyz_frames, xyz_energy


@dataclass(frozen=True)
class NativeObservation:
    version: str
    terminated: bool
    converged: bool
    limited: bool
    stationary: bool
    frame: GeometryFrame | None
    accepted: bool
    steps: tuple[StepObservation, ...] = ()
    diagnostics: tuple[FailureDiagnostic, ...] = ()


def _last_number(pattern, text):
    matches = re.findall(pattern, text, re.I | re.M)
    if not matches:
        raise ValueError('Missing numeric output evidence: ' + pattern)
    return finite_number(matches[-1])


def _error(text, relative, patterns, code):
    for number, line in enumerate(text.splitlines(), 1):
        if any(re.search(pattern, line, re.I) for pattern in patterns):
            return FailureDiagnostic(code, 'backend', line.strip(), relative, number, line.strip()[:300])
    return None


def _version(pattern, text, prefix):
    match = re.search(pattern, text, re.I)
    if not match or match[1] != prefix + '.1':
        raise ValueError('Unsupported or missing output version; diagnostic pilot requires ' + prefix + '.1')
    return match[1]


def _native_orca(request, path):
    relative = 'trial_' + request.name + '.out'
    text = path(relative).read_text()
    version = _version(r'Program Version\s+(\d+\.\d+(?:\.\d+)?)', text, '6.1')
    failed = _error(text, relative, (r'SCF NOT CONVERGED', r'SCF CONVERGENCE FAILURE'), FailureCode.CALC_SCF_MAX_STEPS)
    limited = bool(re.search(r'(maximum number of (?:geometry )?optimization cycles reached|THE OPTIMIZATION DID NOT CONVERGE)', text, re.I))
    terminated = bool(re.search(r'^\s*\*+ORCA TERMINATED NORMALLY\*+\s*$', text, re.M))
    converged = bool(re.search(r'^\s*\**\s*THE OPTIMIZATION HAS CONVERGED\s*\**\s*$', text, re.M))
    if failed:
        return NativeObservation(version, terminated, False, limited, False, None, False, diagnostics=(failed,))
    energy = _last_number(r'FINAL SINGLE POINT ENERGY\s+(' + NUMBER + ')', text)
    frames = read_xyz_frames(path('trial_' + request.name + '.xyz'), request.symbols)
    if len(frames) != 1:
        raise ValueError('Expected a single final ORCA geometry')
    actual, log_coordinates = orca_coordinates(text, len(request.symbols))
    if actual != request.symbols or not close_matrix(log_coordinates, frames[-1][0]):
        raise ValueError('ORCA energy/log geometry and final XYZ do not match')
    # ORCA may converge by overachieved displacement. R1 additionally requires
    # original gradient tolerances for an independent stationarity certificate.
    expected = {'normal': (1e-4, 3e-4), 'loose': (5e-4, 2e-3), 'tight': (3e-5, 1e-4)}[request.settings.get('opt_threshold', 'normal')]
    residuals, stationary = [], True
    for label, threshold in zip(('RMS', 'MAX'), expected):
        matches = re.findall(label + r' gradient\s+(' + NUMBER + r')\s+(' + NUMBER + r')\s+(YES|NO)', text, re.I)
        if not matches:
            stationary = False
            continue
        value, recorded_threshold = map(finite_number, matches[-1][:2])
        if value < 0 or recorded_threshold <= 0:
            raise ValueError('Invalid gradient residual/threshold')
        residuals.append(value)
        stationary &= value <= threshold and recorded_threshold <= threshold * (1+1e-6)
    frame = GeometryFrame(request.symbols, frames[-1][0], energy, energy,
                          source_paths=(relative, 'trial_' + request.name + '.xyz'))
    steps = (StepObservation(0, energy, *residuals, accepted=True),) if len(residuals) == 2 and converged and not limited else ()
    return NativeObservation(version, terminated, converged and not limited, limited, stationary,
                             frame, converged and not limited, steps)


def _native_xtb(request, path):
    relative = 'xtb.out'
    text = path(relative).read_text()
    version = _version(r'(?:\*\s*)?xtb version\s+(\d+\.\d+(?:\.\d+)?)', text, '6.7')
    failed = _error(text, relative, (r'\bSCC (?:failed|not converged)',), FailureCode.CALC_SCF_MAX_STEPS)
    invalid_option = _error(text, relative, (r'unknown option',), FailureCode.USER_INPUT_INVALID)
    if invalid_option:
        failed = invalid_option
    for marker in ('.sccnotconverged', 'NOT_CONVERGED'):
        if path(marker).exists():
            failed = FailureDiagnostic(FailureCode.CALC_SCF_MAX_STEPS, 'backend', 'Current SCC failure marker', marker)
    terminated = bool(re.search(r'^\s*normal termination of xtb\s*$', text, re.I | re.M))
    if failed:
        return NativeObservation(version, terminated, False, False, False, None, False, diagnostics=(failed,))
    models = re.findall(r'Hamiltonian\s+([A-Za-z0-9-]+)', text)
    if not models or models[-1].lower() != {'gfn2': 'gfn2-xtb', 'gxtb': 'g-xtb'}[request.settings.get('xtb_model', 'gfn2')]:
        raise ValueError('xTB actual Hamiltonian does not match request')
    converged = bool(re.search(r'^\s*\*{3}\s*GEOMETRY OPTIMIZATION CONVERGED AFTER \d+ ITERATIONS\s*\*{3}\s*$', text, re.M))
    limited = bool(re.search(r'GEOMETRY OPTIMIZATION (?:FAILED|DID NOT CONVERGE)|maximum.*optimization.*cycles', text, re.I))
    energy = _last_number(r'TOTAL ENERGY\s+(' + NUMBER + ')', text)
    gradient = _last_number(r'GRADIENT NORM\s+(' + NUMBER + ')', text)
    frames = read_xyz_frames(path('xtbopt.xyz'), request.symbols)
    if len(frames) != 1 or abs(xyz_energy(frames[-1][1]) - energy) > 5.1e-7:
        raise ValueError('xTB final geometry energy does not match output energy')
    threshold = {'loose': 4e-3, 'normal': 1e-3, 'tight': 8e-4}[request.settings.get('opt_threshold', 'normal')]
    if gradient < 0:
        raise ValueError('Invalid negative gradient norm')
    frame = GeometryFrame(request.symbols, frames[-1][0], energy, energy,
                          source_paths=('xtb.out', 'xtbopt.xyz'))
    return NativeObservation(version, terminated, converged and not limited, limited, gradient <= threshold,
                             frame, converged and not limited,
                             (StepObservation(0, energy, accepted=converged and not limited, gradient_norm_hartree_per_bohr=gradient),))


def _geometric_steps(text):
    steps = []
    for line in text.splitlines():
        match = re.search(r'^Step\s+(\d+)\s*:', line)
        if match:
            gradient = re.search(r'(?:Gradient|Grad)\s*=\s*(' + NUMBER + r')/(' + NUMBER + ')', line)
            energy = re.search(r'(?:Energy|E \(change\))\s*=\s*(' + NUMBER + ')', line)
            if gradient and energy:
                steps.append(StepObservation(int(match[1]), finite_number(energy[1]),
                                              finite_number(gradient[1]), finite_number(gradient[2])))
        elif 'Rejecting step' in line and steps:
            steps[-1] = replace(steps[-1], accepted=False)
    return steps


def _geometric(request, path):
    text = re.sub(r'\x1b\[[0-9;]*m', '', path('geometric.out').read_text())
    version = _version(r'geomeTRIC started\. Version:\s*(\d+\.\d+(?:\.\d+)?)', text, '1.1')
    failed = _error(text, 'geometric.out', (r'\bSCC failed', r'SCF NOT CONVERGED'), FailureCode.CALC_SCF_MAX_STEPS)
    if failed:
        return NativeObservation(version, False, False, False, False, None, False, diagnostics=(failed,))
    if request.settings.get('_ore_safe') and 'Maximum iterations reached' in text:
        from .checkpoints import accepted_checkpoint
        frame, steps = accepted_checkpoint(request, path)
        return NativeObservation(version, False, False, True, False, frame, True, steps)
    state = read_json(path('pyar_geometric_state.json'))
    if state.get('schema_version', 1) != 1:
        raise ValueError('Unknown geomeTRIC state schema version')
    if state.get('software') != request.backend or state.get('opt_target') != request.target:
        raise ValueError('geomeTRIC calculator identity does not match request')
    if float(state.get('gamma', 0)) != float(request.settings.get('gamma') or 0):
        raise ValueError('geomeTRIC bias objective does not match request')
    if state.get('bias_potential', 'afir') != request.settings.get('bias_potential', 'afir'):
        raise ValueError('geomeTRIC bias functional does not match request')
    if float(state.get('softmin_beta', 1.0)) != float(request.settings.get('softmin_beta', 1.0)):
        raise ValueError('geomeTRIC bias beta does not match request')
    if state.get('wall_configuration') != request.settings.get('microsolvation_wall'):
        raise ValueError('geomeTRIC wall objective does not match request')
    state_coordinates = matrix(state['positions_angstrom'], len(request.symbols))
    energy = finite_number(state['backend_energy_hartree'])
    effective = finite_number(state['total_energy_hartree'])
    if not (request.settings.get('gamma') or request.settings.get('microsolvation_wall')) and abs(energy-effective) > 1e-10:
        raise ValueError('Unbiased physical and effective energies disagree')
    frames = read_xyz_frames(path('trial_' + request.name + '_optim.xyz'), request.symbols)
    coordinates, comment = frames[-1]
    if not close_matrix(coordinates, state_coordinates) or abs(xyz_energy(comment) - effective) > 5.1e-8:
        raise ValueError('geomeTRIC XYZ, physical energy and objective refer to different frames')
    forces = matrix(state['forces_ev_per_angstrom'], len(request.symbols))
    frame = GeometryFrame(request.symbols, coordinates, energy, effective, forces,
                          ('geometric.out', 'pyar_geometric_state.json', 'trial_' + request.name + '_optim.xyz'))
    limited = 'Maximum iterations reached' in text or state.get('optimization_status') == 'cycle_exceeded'
    converged = bool(re.search(r'^Converged!', text, re.M)) and not limited
    steps = _geometric_steps(text)
    if not steps or abs(steps[-1].effective_energy_hartree - effective) > 5.1e-8:
        raise ValueError('Missing matching optimizer energy/gradient observation')
    # CODATA unit conversion; tolerate the optimizer log's four-digit rounding.
    gradient_norms = [math.sqrt(sum(x*x for x in force)) * 0.529177210903 / 27.211386245988 for force in forces]
    actual_rms = math.sqrt(sum(g*g for g in gradient_norms) / len(gradient_norms))
    actual_max = max(gradient_norms)
    if not all(math.isclose(actual, logged, rel_tol=6e-4, abs_tol=6e-8) for actual,logged in
               ((actual_rms, steps[-1].gradient_rms_hartree_per_bohr), (actual_max, steps[-1].gradient_max_hartree_per_bohr))):
        raise ValueError('geomeTRIC stored forces and logged gradients disagree')
    rms, maximum = {'loose': (1.7e-3, 2.5e-3), 'normal': (3e-4, 4.5e-4), 'tight': (1e-5, 1.5e-5)}[request.settings.get('opt_threshold', 'normal')]
    stationary = (steps[-1].gradient_rms_hartree_per_bohr <= rms
                  and steps[-1].gradient_max_hartree_per_bohr <= maximum
                  and actual_rms <= rms and actual_max <= maximum)
    accepted = converged
    checkpoint = path('pyar_bias_controller_state.json')
    if not converged and checkpoint.exists():
        value = read_json(checkpoint)
        configuration = value['configuration']
        settings = request.settings
        if settings.get('bias_controller') != 'adaptive' or not settings.get('gamma'):
            raise ValueError('No accepted adaptive-step contract for this request')
        # Ask the existing owner to validate its numerical state, without
        # changing any files, resuming optimization or granting restart support.
        from pyar.biases.controller import BiasController
        from pyar.biases.afir import alpha_from_gamma
        controller = BiasController('adaptive', alpha_min=settings.get('bias_alpha_min', 0.),
            safety_margin=settings.get('bias_alpha_margin'), smoothing=settings.get('bias_alpha_smoothing', 1.),
            epsilon=settings.get('bias_alpha_epsilon', 1e-12), scheduled_alpha=settings.get('bias_scheduled_alpha'))
        controller.load_state_dict(value['controller'], alpha_from_gamma(settings['gamma']))
        if value['controller'] != state.get('bias_controller'):
            raise ValueError('Accepted controller and evaluated-frame objectives disagree')
        expected = request.settings
        configured = configuration['qc_params']
        # Existing checkpoint intentionally omits these execution-only flags.
        for key in ('bias_controller_restart', 'trace_mode', 'release_retry_attempt', 'release_job_name', 'adaptive_max_segments', 'adaptive_suppress_release_candidate'):
            expected.pop(key, None)
        expected.pop('fragments', None)
        if (value['symbols'] != list(request.symbols) or configured != expected
                or configuration.get('fragment_indices') != request.settings.get('fragments')
                or configuration.get('opt_target') != request.target
                or not close_matrix(matrix(value['positions_angstrom'], len(request.symbols)), coordinates)):
            raise ValueError('Incompatible accepted-controller checkpoint')
        accepted = True
    if accepted:
        steps[-1] = replace(steps[-1], accepted=True)
    if request.settings.get('_ore_safe'):
        from .checkpoints import accepted_checkpoint
        checkpoint_frame, accepted_steps = accepted_checkpoint(request, path)
        if (not close_matrix(checkpoint_frame.coordinates_angstrom, frame.coordinates_angstrom)
                or abs(checkpoint_frame.physical_energy_hartree-energy) > 1e-9):
            raise ValueError('Converged frame and accepted checkpoint disagree')
        steps = accepted_steps
    return NativeObservation(version, converged, converged, limited, stationary, frame, accepted, tuple(steps))


def parse_observation(request, path):
    if request.optimizer == 'geometric':
        return _geometric(request, path)
    if request.backend == 'orca':
        return _native_orca(request, path)
    if request.backend == 'xtb':
        return _native_xtb(request, path)
    raise ValueError('Unsupported diagnostic route')
