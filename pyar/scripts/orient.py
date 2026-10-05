"""Geometry-only encounter orientation command."""
import argparse

from pyar.geometry_utilities import DEFAULT_ORIENTATIONS, orient_structures


def build_parser(prog=None):
    parser = argparse.ArgumentParser(prog=prog or 'pyar orient', description=__doc__)
    parser.add_argument('first', metavar='A.xyz')
    parser.add_argument('second', metavar='B.xyz')
    parser.add_argument('--orientations', '-N', type=int, default=DEFAULT_ORIENTATIONS)
    parser.add_argument('--distance-scale', type=float)
    parser.add_argument('--sequence-offset', type=int, default=0)
    parser.add_argument('--output', default='orientations')
    return parser


def main(argv=None, *, prog=None):
    parser = build_parser(prog)
    args = parser.parse_args(argv)
    try:
        result = orient_structures(args.first, args.second, orientations=args.orientations,
                                   distance_scale=args.distance_scale, sequence_offset=args.sequence_offset,
                                   output=args.output)
    except (ValueError, OSError) as exc:
        parser.error(str(exc))
    print(f"Generated {len(result['files'])} encounter geometries in {args.output}")
    print(f"Trial vectors: {result['trial_vectors']}")
