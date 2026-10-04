"""Side-effect-free pairwise structure inspection command."""

import argparse
import json

from pyar.structure_inspection import compare_structures


def build_parser(prog=None):
    parser = argparse.ArgumentParser(
        prog=prog or 'pyar compare',
        description='Compare XYZ coordinates and stored energies. ΔE = E(B) - E(A).')
    parser.add_argument('first', metavar='A.xyz')
    parser.add_argument('second', metavar='B.xyz')
    parser.add_argument('--atom-mode', choices=('heavy', 'all'), default='heavy')
    parser.add_argument('--bond-scale', type=float, default=1.15,
                        help='Covalent-radius adjacency scale (default: 1.15)')
    parser.add_argument('--maximum-mappings', type=int, default=10000)
    parser.add_argument('--rmsd-threshold', type=float,
                        help='Explicit geometry equivalence threshold in Å (strictly below)')
    parser.add_argument('--json', action='store_true', help='Print machine-readable JSON')
    return parser


def _yes(value):
    return 'not applicable' if value is None else 'yes' if value else 'no'


def render_comparison(result):
    """Render the same structured result used for JSON."""
    r = result
    print(f"Comparison: {r['first_file']} ↔ {r['second_file']}")
    print('\nComposition')
    print(f"  Atoms: {r['first_atom_count']} / {r['second_atom_count']}")
    for side in ('first', 'second'):
        formula = ''.join(element + (str(count) if count != 1 else '')
                          for element, count in r[f'{side}_composition'].items())
        print(f"  {'A' if side == 'first' else 'B'}: {formula}")
    print(f"  Same composition: {_yes(r['same_composition'])}")
    print(f"  Same atom order (element sequence): {_yes(r['same_atom_order'])}")
    print('\nEnergy')
    for label, key in [('A', 'first_energy_hartree'), ('B', 'second_energy_hartree')]:
        energy = r[key]
        print(f"  {label}: " + ('unavailable' if energy is None else f'{energy:.6f} Eh'))
    delta = r['delta_energy_kcal_mol']
    print('  ΔE (B - A): ' + ('unavailable' if delta is None else f'{delta:+.2f} kcal/mol'))
    if delta is not None:
        print('  Lower energy: ' + ('A' if delta > 0 else 'B' if delta < 0 else 'equal'))
    print('\nGeometry')
    print(f"  Connectivity match: {_yes(r['connectivity_match'])}")
    rmsd = r['rmsd_angstrom']
    missing = ('incomplete' if r['comparison_complete'] is False else
               'not applicable under graph-equivalent comparison')
    print(f"  {r['rmsd_atom_mode'].capitalize()}-atom RMSD: " +
          (missing if rmsd is None else f'{rmsd:.6f} Å'))
    print(f"  Comparison complete: {_yes(r['comparison_complete'])}")
    print(f"  Mappings evaluated: {r['mappings_evaluated']}")
    if r['rmsd_threshold_angstrom'] is not None:
        print(f"  Equivalent below {r['rmsd_threshold_angstrom']} Å: "
              f"{_yes(r['geometry_equivalent_under_threshold'])}")
    print('\nConnectivity (inferred adjacency; indices are 0-based)')
    print(f"  Components: {r['input_component_count']} → {r['output_component_count']}")
    print(f"  Inferred edges: {r['input_edge_count']} → {r['output_edge_count']}")
    if r['same_atom_order']:
        for key in ('added_edges', 'removed_edges'):
            edges = r[key]
            print(f"  {key.split('_')[0].capitalize()} inferred contacts/edges:")
            if not edges:
                print('    none')
            for left, right in edges:
                labels = r['atom_labels']
                print(f'    {left}({labels[left]}) — {right}({labels[right]})')
    else:
        print('  Specific added/removed edge indices were not reported because atom ordering differs.')
    print('\nLimitations')
    for limitation in r['limitations']:
        print('  ' + limitation)


def main(argv=None, *, prog=None):
    parser = build_parser(prog)
    args = parser.parse_args(argv)
    try:
        result = compare_structures(args.first, args.second, atom_mode=args.atom_mode,
                                   bond_scale=args.bond_scale,
                                   maximum_mappings=args.maximum_mappings,
                                   rmsd_threshold=args.rmsd_threshold)
    except (ValueError, OSError) as exc:
        parser.error(str(exc))
    if args.json:
        print(json.dumps(result, allow_nan=False))
    else:
        render_comparison(result)
