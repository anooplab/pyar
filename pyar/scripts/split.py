"""Split coordinate-connected components; no bond orders or electronic states assigned."""
import argparse

from pyar.geometry_utilities import split_structure


def build_parser(prog=None):
    parser = argparse.ArgumentParser(prog=prog or 'pyar split', description=__doc__)
    parser.add_argument('input', metavar='complex.xyz')
    parser.add_argument('--output')
    parser.add_argument('--bond-scale', type=float, default=1.15)
    return parser


def main(argv=None, *, prog=None):
    parser = build_parser(prog)
    args = parser.parse_args(argv)
    try:
        result = split_structure(args.input, output=args.output, bond_scale=args.bond_scale)
    except (ValueError, OSError) as exc:
        parser.error(str(exc))
    if result['component_count'] == 1:
        print(f'{args.input} contains one connected component; nothing to split.')
    else:
        print(f"Extracted {result['component_count']} coordinate components (not formal chemical fragments):")
        for filename in result['files']:
            print('  ' + filename)
