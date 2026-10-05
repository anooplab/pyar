"""Small shared input and safe-output helpers for XYZ utilities."""

from pathlib import Path
import shutil

from pyar.structure_inspection import load_structure, optional_energy


def load_structures(files, *, require_energy=False):
    molecules = []
    for filename in files:
        molecule = load_structure(filename)
        molecule.name = str(Path(filename))
        molecule.relative_path = str(Path(filename))
        molecule.energy = optional_energy(filename)
        if require_energy and molecule.energy is None:
            raise ValueError(f'Could not read an energy from: {filename}\n'
                             'Expected a numeric energy in the XYZ comment line.')
        molecules.append(molecule)
    return molecules


def expand_charges(charges, count):
    if charges is None:
        return [None] * count
    charges = list(charges)
    if len(charges) == 1:
        return charges * count
    if len(charges) != count:
        raise ValueError('--charge requires one value or one value per input file')
    return charges


def validate_output(directory):
    path = Path(directory)
    if path.is_symlink() or (path.exists() and (not path.is_dir() or any(path.iterdir()))):
        raise ValueError(f'Output destination {path} must be an empty directory or a new path; '
                         'choose another --output directory.')
    return path


def copy_structures(molecules, directory):
    """Copy originals only after all destinations and sources have been checked."""
    sources = [Path(molecule.relative_path) for molecule in molecules]
    names = [source.name for source in sources]
    if len(set(names)) != len(names):
        raise ValueError('Retained inputs have colliding basenames; choose distinct input filenames')
    path = validate_output(directory)
    for source in sources:
        if not source.is_file():
            raise ValueError(f'Input file is unavailable: {source}')
    if not sources:
        return []
    path.mkdir(parents=True, exist_ok=True)
    for source in sources:
        shutil.copyfile(source, path / source.name)
    return [str(path / name) for name in names]


def write_coordinates(filename, atoms, coordinates, title):
    """Write coordinate-only XYZ without fabricated energy or electronic state."""
    with Path(filename).open('x') as stream:
        stream.write(f'{len(atoms)}\n{title}\n')
        for atom, (x, y, z) in zip(atoms, coordinates):
            stream.write(f'{atom} {x:.17g} {y:.17g} {z:.17g}\n')
