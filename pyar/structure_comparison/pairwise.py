"""Reusable molecular pairwise analysis, independent of file/reporting layers."""

from collections import Counter

from pyar.structure_comparison.chemical_identity import perceive_chemical_identity
from pyar.structure_comparison.coordinate_graph import compare_coordinate_structures
from pyar.structure_comparison.graph_rmsd import GraphRMSDComparator


def compare_structures(first, second, *, atom_mode="heavy", bond_scale=1.15,
                       maximum_mappings=10000, rmsd_threshold=None, charges=None,
                       known_charges=None):
    """Compare validated Molecules with independent coordinate/chemical evidence.

    ``known_charges`` must come from authoritative input/state, never default
    Molecule metadata. XYZ inspection does not provide reliable charge metadata.
    """
    a, b = first, second
    comparator = GraphRMSDComparator(atom_mode=atom_mode, bond_scale=bond_scale,
                                    max_isomorphisms=maximum_mappings,
                                    threshold=rmsd_threshold)
    connectivity = compare_coordinate_structures(a, b, scale=bond_scale)
    comparable = connectivity['same_composition'] and connectivity['same_atom_count']
    geometry = comparator.compare(a, b) if comparable else None
    charges = (None, None) if charges is None else tuple(charges)
    known_charges = (None, None) if known_charges is None else tuple(known_charges)
    if len(charges) != 2 or len(known_charges) != 2:
        raise ValueError('Expected one resolved charge per structure')
    identities = [perceive_chemical_identity(molecule, charge, known_charge=known)
                  for molecule, charge, known in zip((a, b), charges, known_charges)]
    successful = all(identity['perception_success'] for identity in identities)
    effective_mode = ('all' if all(atom == 'H' for atom in a.atoms_list) else atom_mode)
    return {
        'first_atom_count': len(a.atoms_list), 'second_atom_count': len(b.atoms_list),
        'first_composition': dict(sorted(Counter(a.atoms_list).items())),
        'second_composition': dict(sorted(Counter(b.atoms_list).items())),
        'atom_labels': a.atoms_list if connectivity['same_atom_order'] else None,
        **connectivity,
        'chemical_identity': {
            'first': identities[0], 'second': identities[1],
            'canonical_smiles_match': (identities[0]['canonical_smiles'] ==
                                       identities[1]['canonical_smiles']) if successful else None,
        },
        'connectivity_match': None if geometry is None else geometry.metadata.get('connectivity_match'),
        'rmsd_angstrom': None if geometry is None else geometry.distance,
        'rmsd_atom_mode': effective_mode,
        'rmsd_threshold_angstrom': rmsd_threshold,
        'geometry_equivalent_under_threshold': None if geometry is None else geometry.equivalent,
        'comparison_complete': None if geometry is None else geometry.metadata.get('comparison_complete', True),
        'mappings_evaluated': 0 if geometry is None else geometry.metadata.get('isomorphisms_evaluated', 0),
        'coordinate_model': 'covalent-radii', 'bond_scale': bond_scale,
    }
