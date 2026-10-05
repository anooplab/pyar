"""Read-only pairwise inspection using PyAR's coordinate graph and RMSD APIs."""

import math
from pathlib import Path

import numpy as np

from pyar.core.molecule import Molecule, parse_xyz
from pyar.selection import reports
from pyar.structure_comparison.pairwise import compare_structures as compare_molecules


def load_structure(filename):
    """Load usable coordinates without assigning an electronic state."""
    atoms, coordinates, name, title, _ = parse_xyz(filename)
    if not np.isfinite(coordinates).all():
        raise ValueError(f"Non-finite coordinates in {filename}")
    try:
        return Molecule(atoms, coordinates, name=name, title=title,
                        charge=None, multiplicity=None)
    except KeyError as exc:
        raise ValueError(f"Unknown element in {filename}: {exc}") from exc


def optional_energy(filename):
    """Missing comment-line energies do not prevent structural inspection."""
    try:
        energy = reports.read_energy_from_xyz_file(filename)
        return energy if math.isfinite(energy) else None
    except (ValueError, IndexError):
        return None


def compare_structures(first, second, *, atom_mode="heavy", bond_scale=1.15,
                       maximum_mappings=10000, rmsd_threshold=None, charges=None):
    """Compare two XYZ paths; return JSON-ready data independently of rendering.

    Indexed edge changes follow the coordinate API's identical-element-order
    convention. XYZ cannot verify correspondence of repeated elements.
    """
    a, b = load_structure(first), load_structure(second)
    if charges is not None:
        charges = list(charges)
        if len(charges) == 1:
            charges *= 2
        if len(charges) != 2:
            raise ValueError('--charge requires one value or one value per structure (two values)')
    structural = compare_molecules(a, b, atom_mode=atom_mode, bond_scale=bond_scale,
                                  maximum_mappings=maximum_mappings,
                                  rmsd_threshold=rmsd_threshold, charges=charges)
    comparable = structural['same_composition'] and structural['same_atom_count']
    chemical = structural['chemical_identity']
    same_charge = chemical['first']['charge_used'] == chemical['second']['charge_used']
    ea, eb = optional_energy(first), optional_energy(second)
    delta = None if ea is None or eb is None else eb - ea
    limitations = [
        "Connectivity is inferred from XYZ coordinates using covalent-radius adjacency; "
        "edges do not establish bond orders, chemical identity, charge, or spin.",
        "Same atom order means the same element sequence; indexed edge changes assume "
        "correspondence by row, which XYZ cannot verify for repeated elements.",
    ]
    limitations.append("Canonical SMILES are perceived from coordinates and the stated/assumed charge; "
                       "bond orders and stereochemistry are model-dependent, not definitive identity.")
    if not comparable:
        limitations.append("Absolute electronic energies of different compositions are not "
                           "directly interpretable as relative isomer/conformer energies.")
    if not same_charge:
        limitations.append('Energies for different supplied/assumed charges are not relative '
                           'isomer/conformer energies, even when composition matches.')
    if structural['comparison_complete'] is False:
        limitations.append("Graph-isomorphism limit reached; RMSD comparison is incomplete.")
    return {
        'first_file': str(Path(first)), 'second_file': str(Path(second)),
        **structural,
        'first_energy_hartree': ea, 'second_energy_hartree': eb,
        'delta_energy_hartree': delta,
        'delta_energy_kcal_mol': None if delta is None else delta * reports.HARTREE_TO_KCAL_MOL,
        'energy_difference_is_isomer_comparison': comparable and same_charge,
        'limitations': limitations,
    }


def format_formula(composition):
    """Hill formula: C/H first when carbon is present, otherwise alphabetical."""
    order = sorted(composition)
    if 'C' in composition:
        order = ['C'] + (['H'] if 'H' in composition else []) + [
            element for element in order if element not in {'C', 'H'}]
    return ''.join(element + (str(composition[element]) if composition[element] != 1 else '')
                   for element in order)


def identify_structure(molecule, *, charge=None, bond_scale=1.15):
    """Enrich an already loaded geometry without making RDKit mandatory."""
    from pyar.structure_comparison.coordinate_graph import analyze_coordinate_structure
    from pyar.structure_comparison.chemical_identity import perceive_chemical_identity

    topology = analyze_coordinate_structure(molecule, scale=bond_scale)
    return {'file': str(getattr(molecule, 'relative_path', molecule.name)),
            'atom_count': molecule.number_of_atoms, 'composition': topology['composition'],
            'formula': format_formula(topology['composition']),
            'energy_hartree': molecule.energy, 'component_count': topology['component_count'],
            'coordinate_topology': topology,
            'chemical_identity': perceive_chemical_identity(molecule, charge)}
