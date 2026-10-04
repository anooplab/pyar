"""Read-only pairwise inspection using PyAR's coordinate graph and RMSD APIs."""

from collections import Counter
import math
from pathlib import Path

import numpy as np

from pyar.core.molecule import Molecule, parse_xyz
from pyar.selection import reports
from pyar.structure_comparison.coordinate_graph import compare_coordinate_structures
from pyar.structure_comparison.graph_rmsd import GraphRMSDComparator, selected_atom_indices


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
                       maximum_mappings=10000, rmsd_threshold=None):
    """Compare two XYZ paths; return JSON-ready data independently of rendering.

    Indexed edge changes follow the coordinate API's identical-element-order
    convention. XYZ cannot verify correspondence of repeated elements.
    """
    comparator = GraphRMSDComparator(atom_mode=atom_mode, bond_scale=bond_scale,
                                    max_isomorphisms=maximum_mappings,
                                    threshold=rmsd_threshold)
    a, b = load_structure(first), load_structure(second)
    connectivity = compare_coordinate_structures(a, b, scale=bond_scale)
    comparable = connectivity['same_composition'] and connectivity['same_atom_count']
    geometry = comparator.compare(a, b) if comparable else None
    ea, eb = optional_energy(first), optional_energy(second)
    delta = None if ea is None or eb is None else eb - ea
    limitations = [
        "Connectivity is inferred from XYZ coordinates using covalent-radius adjacency; "
        "edges do not establish bond orders, chemical identity, charge, or spin.",
        "Same atom order means the same element sequence; indexed edge changes assume "
        "correspondence by row, which XYZ cannot verify for repeated elements.",
    ]
    if not comparable:
        limitations.append("Absolute electronic energies of different compositions are not "
                           "directly interpretable as relative isomer/conformer energies.")
    if geometry is not None and not geometry.metadata.get('comparison_complete', True):
        limitations.append("Graph-isomorphism limit reached; RMSD comparison is incomplete.")
    effective_mode = ('all' if len(selected_atom_indices(a.atoms_list, atom_mode))
                      == len(a.atoms_list) and all(atom == 'H' for atom in a.atoms_list)
                      else atom_mode)
    return {
        'first_file': str(Path(first)), 'second_file': str(Path(second)),
        'first_atom_count': len(a.atoms_list), 'second_atom_count': len(b.atoms_list),
        'first_composition': dict(sorted(Counter(a.atoms_list).items())),
        'second_composition': dict(sorted(Counter(b.atoms_list).items())),
        'atom_labels': a.atoms_list if connectivity['same_atom_order'] else None,
        **connectivity,
        'first_energy_hartree': ea, 'second_energy_hartree': eb,
        'delta_energy_hartree': delta,
        'delta_energy_kcal_mol': None if delta is None else delta * reports.HARTREE_TO_KCAL_MOL,
        'energy_difference_is_isomer_comparison': comparable,
        'connectivity_match': None if geometry is None else geometry.metadata.get('connectivity_match'),
        'rmsd_angstrom': None if geometry is None else geometry.distance,
        'rmsd_atom_mode': effective_mode,
        'rmsd_threshold_angstrom': rmsd_threshold,
        'geometry_equivalent_under_threshold': None if geometry is None else geometry.equivalent,
        'comparison_complete': None if geometry is None else geometry.metadata.get('comparison_complete', True),
        'mappings_evaluated': 0 if geometry is None else geometry.metadata.get('isomorphisms_evaluated', 0),
        'coordinate_model': 'covalent-radii', 'bond_scale': bond_scale,
        'limitations': limitations,
    }
