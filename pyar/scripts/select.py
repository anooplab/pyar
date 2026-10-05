"""Select stored-energy structures, with optional graph-first deduplication."""
import argparse
import sys

from pyar.selection.energy_window import select_structures
from pyar.selection.reports import print_energy_table
from pyar.utility_io import copy_structures, load_structures


def build_parser(prog=None):
    parser = argparse.ArgumentParser(prog=prog or 'pyar select', description=__doc__)
    parser.add_argument('inputs', nargs='+', metavar='xyz')
    parser.add_argument('--within', type=float, help='Inclusive energy window in kcal/mol')
    parser.add_argument('--top', type=int)
    parser.add_argument('--unique', action='store_true')
    parser.add_argument('--rmsd-threshold', type=float)
    parser.add_argument('--output')
    return parser


def main(argv=None, *, prog=None):
    parser = build_parser(prog)
    args = parser.parse_args(argv)
    if args.rmsd_threshold is not None and not args.unique:
        parser.error('--rmsd-threshold requires --unique')
    try:
        result = select_structures(load_structures(args.inputs, require_energy=True),
                                   within=args.within, top=args.top, unique=args.unique,
                                   threshold=args.rmsd_threshold)
        if args.output:
            copy_structures(result['kept'], args.output)
    except (ValueError, OSError) as exc:
        parser.error(str(exc))
    print_energy_table(result['kept'], stream=sys.stdout, title='Selected structures (absolute energy: Eh):')
    print('Selected paths:')
    for molecule in result['kept']:
        print('  ' + molecule.relative_path)
    if result['deduplication'] is not None:
        for diagnostic in result['deduplication']['diagnostics']:
            print('  Retention diagnostic: ' + diagnostic)
