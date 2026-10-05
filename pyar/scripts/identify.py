"""Inspect XYZ composition, coordinate topology and optional perceived identity."""
import argparse
import json

from pyar.structure_inspection import identify_structure
from pyar.utility_io import expand_charges, load_structures


def build_parser(prog=None):
    parser = argparse.ArgumentParser(prog=prog or 'pyar identify', description=__doc__)
    parser.add_argument('inputs', nargs='+', metavar='xyz')
    parser.add_argument('--charge', nargs='+', type=int, help='One charge for all files or one per file')
    parser.add_argument('--bond-scale', type=float, default=1.15)
    parser.add_argument('--json', action='store_true')
    return parser


def main(argv=None, *, prog=None):
    parser = build_parser(prog)
    args = parser.parse_args(argv)
    try:
        charges = expand_charges(args.charge, len(args.inputs))
        molecules = load_structures(args.inputs)
        results = [identify_structure(molecule, charge=charge, bond_scale=args.bond_scale)
                   for molecule, charge in zip(molecules, charges)]
    except (ValueError, OSError) as exc:
        parser.error(str(exc))
    if args.json:
        print(json.dumps({'structures': results}, allow_nan=False))
        return
    for result in results:
        identity = result['chemical_identity']
        print(f"\n{result['file']}\n  Atoms: {result['atom_count']}\n  Formula: {result['formula']}")
        energy = result['energy_hartree']
        print('  Energy: ' + ('unavailable' if energy is None else f'{energy:.6f} Eh'))
        print(f"  Coordinate components: {result['component_count']} (geometric adjacency only)")
        print('  Canonical SMILES: ' + (identity['canonical_smiles'] or 'unavailable'))
        print(f"  Charge: {identity['charge_used']} ({identity['charge_source']})")
        print(f"  Perception: RDKit DetermineBonds (xyz2mol), {identity['status']}")
        if identity['reason']:
            print('  Reason: ' + identity['reason'])
        if identity['installation_hint']:
            print('  Enable perception: ' + identity['installation_hint'])
        print('  XYZ does not establish charge, spin, formal bond orders or complete chemical identity.')
