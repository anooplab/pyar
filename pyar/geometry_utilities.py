"""Geometry-only encounter generation and coordinate-component extraction."""

from pathlib import Path

import networkx as nx
import numpy as np

from pyar.data.defualt_parameters import values
from pyar.structure_inspection import load_structure
from pyar.structure_comparison.coordinate_graph import infer_coordinate_graph
from pyar.utility_io import validate_output, write_coordinates


DEFAULT_ORIENTATIONS = values["how_many_orientations"]


def orient_structures(first, second, *, orientations=DEFAULT_ORIENTATIONS, output='orientations',
                      distance_scale=None, sequence_offset=0):
    from pyar.sampling.trial_generator import (
        generate_trial_vectors, merge_two_molecules, write_trial_vectors,
    )
    if not isinstance(orientations, int) or orientations < 1:
        raise ValueError('--orientations must be a positive integer')
    if not isinstance(sequence_offset, int) or sequence_offset < 0:
        raise ValueError('--sequence-offset must be a nonnegative integer')
    a, b = load_structure(first), load_structure(second)
    scale = (1.2 if b.number_of_atoms == 1 else 1.5) if distance_scale is None else distance_scale
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError('--distance-scale must be finite and positive')
    directory = validate_output(output)
    vectors = generate_trial_vectors(orientations, direction_method='fibonacci',
                                     rotation_method='halton', sequence_offset=sequence_offset,
                                     use_angles=b.number_of_atoms > 1)
    geometries = [merge_two_molecules(vector, a, b, distance_scaling=scale) for vector in vectors]
    directory.mkdir(parents=True, exist_ok=True)
    write_trial_vectors(vectors, directory / 'trial_vectors.dat')
    files = []
    for index, molecule in enumerate(geometries):
        filename = directory / f'orientation_{index:03d}.xyz'
        write_coordinates(filename, molecule.atoms_list, molecule.coordinates,
                          f'Encounter geometry from {first} and {second}; energy unavailable')
        files.append(str(filename))
    return {'files': files, 'trial_vectors': str(directory / 'trial_vectors.dat'),
            'distance_scale': scale, 'sequence_offset': sequence_offset}


def split_structure(filename, *, output=None, bond_scale=1.15):
    molecule = load_structure(filename)
    graph = infer_coordinate_graph(molecule, scale=bond_scale)
    components = sorted((sorted(component) for component in nx.connected_components(graph)),
                        key=lambda component: component[0])
    if len(components) == 1:
        return {'input': str(filename), 'component_count': 1, 'files': [], 'components': components}
    stem = Path(filename).stem
    directory = validate_output(output or f'{stem}_fragments')
    directory.mkdir(parents=True, exist_ok=True)
    files = []
    for number, indices in enumerate(components, 1):
        path = directory / f'{stem}_fragment_{number:03d}.xyz'
        write_coordinates(path, [molecule.atoms_list[index] for index in indices],
                          molecule.coordinates[indices],
                          f'component {number} from {filename}; original atom indices: '
                          + ', '.join(map(str, indices)) + '; energy unavailable')
        files.append(str(path))
    return {'input': str(filename), 'component_count': len(components),
            'files': files, 'components': components}
