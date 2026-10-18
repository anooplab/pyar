"""Validate the version-pinned driver's accepted-frame observations, read-only."""
import math
import re

from pyar.persistence import read_json
from .models import GeometryFrame, StepObservation, matrix
from .validation import close_matrix, finite_number


def accepted_checkpoint(request, path):
    evaluations = read_json(path('ore_evaluations.json'))['used']
    if type(evaluations) is not int or not 1 <= evaluations <= request.settings['_ore_evaluation_limit']:
        raise ValueError('Invalid physical evaluation count')
    for number in range(1, evaluations+1):
        output = path(f'ore_backend_{number:06d}.out').read_text()
        if (not re.search(r'xtb version\s+6\.7\.1\b', output, re.I)
                or not re.search(r'Hamiltonian\s+GFN2-xTB\b', output, re.I)
                or re.search(r'SCC (?:failed|not converged)', output, re.I)
                or 'normal termination of xtb' not in output.lower()):
            raise ValueError('Unverified physical backend evaluation')
    manifest = read_json(path('ore_accepted.json'))
    files = manifest['files']
    if manifest.get('schema_version') != 1 or not isinstance(files, list) or not 1 <= len(files) <= 10001:
        raise ValueError('Invalid accepted-step manifest')
    steps, last_frame, previous = [], None, -1
    for filename in files:
        value = read_json(path(filename))
        index = value['step']
        if (type(index) is not int or index <= previous or index > request.settings['opt_cycles']
                or filename != f'ore_accepted_{index:06d}.json'
                or value.get('schema_version') != 1 or value.get('version') != '1.1.1'
                or value['objective_fingerprint'] != request.objective_fingerprint
                or tuple(value['symbols']) != request.symbols):
            raise ValueError('Incompatible accepted-step checkpoint')
        positions = matrix(value['coordinates_angstrom'], len(request.symbols))
        gradients = matrix(value['gradient_hartree_per_bohr'], len(request.symbols))
        norms = [sum(x*x for x in row) for row in gradients]
        rms, maximum = math.sqrt(sum(norms)/len(norms)), math.sqrt(max(norms))
        if (abs(rms-finite_number(value['gradient_rms'])) > 1e-12
                or abs(maximum-finite_number(value['gradient_max'])) > 1e-12):
            raise ValueError('Checkpoint forces and residuals disagree')
        if previous == -1 and (index != 0 or not close_matrix(positions, request.coordinates_angstrom)):
            raise ValueError('Checkpoint initial frame does not match request')
        energy = finite_number(value['energy_hartree'])
        forces = tuple(tuple(-x*27.211386245988/.529177210903 for x in row) for row in gradients)
        last_frame = GeometryFrame(request.symbols, positions, energy, energy, forces, (filename,))
        trust = value.get('trust_radius_bohr')
        if trust is not None:
            trust = finite_number(trust)
            if trust <= 0:
                raise ValueError('Invalid checkpoint trust radius')
        steps.append(StepObservation(index, energy, rms, maximum, accepted=True,
                                     trust_radius_bohr=trust))
        previous = index
    return last_frame, tuple(steps)
