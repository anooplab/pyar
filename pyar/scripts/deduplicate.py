"""Inspect unique geometries with PyAR's conservative graph-first policy."""
import argparse

from pyar.selection.deduplication import deduplicate_structures
from pyar.utility_io import copy_structures, load_structures


def build_parser(prog=None):
    parser = argparse.ArgumentParser(prog=prog or 'pyar deduplicate', description=__doc__)
    parser.add_argument('inputs', nargs='+', metavar='xyz')
    parser.add_argument('--rmsd-threshold', type=float)
    parser.add_argument('--output')
    return parser


def main(argv=None, *, prog=None):
    parser = build_parser(prog)
    args = parser.parse_args(argv)
    try:
        result = deduplicate_structures(load_structures(args.inputs), threshold=args.rmsd_threshold)
        if args.output:
            copy_structures(result['kept'], args.output)
    except (ValueError, OSError) as exc:
        parser.error(str(exc))
    print(f"Deduplication\n  Input structures: {result['input_count']}\n"
          f"  Unique structures: {len(result['kept'])}\n  Duplicates removed: {len(result['removed'])}")
    print(f"  RMSD threshold: {result['rmsd_threshold']:.6f} Å\nKept")
    for molecule in result['kept']:
        print('  ' + molecule.relative_path)
    print('Removed')
    for candidate, kept, distance in result['removed']:
        print(f'  {candidate} duplicate of {kept}; RMSD {distance:.6f} Å')
    for diagnostic in result['diagnostics']:
        print('  Retention diagnostic: ' + diagnostic)
