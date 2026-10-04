#!/usr/bin/env python3
"""Importable ``pyar-clustering`` entrypoint.

This utility clusters or filters XYZ pools using the selection algorithms
implemented in :mod:`pyar.selection.clustering`. It prints the energy
table for the input pool and then emits the selected geometries.
"""

import argparse
import csv
import json
import sys

from pyar.selection import clustering
from pyar.selection import reports as selection_reports
from pyar.core.molecule import Molecule
from pyar.selection.clusterers import CLUSTERING_ALGORITHMS, cluster_molecules, validate_clustering_options
from pyar.selection.distances import DISTANCE_METRICS
from pyar.selection.features import FEATURES
from pyar.selection.policy import SYSTEM_TYPES, classify_system_pool, resolve_clustering_policy
from pyar.structure_comparison.coordinate_graph import analyze_coordinate_structure


def build_parser(prog=None):
    """Build the shared parser for the standalone and task-oriented interfaces."""
    parser = argparse.ArgumentParser(prog=prog)
    parser.add_argument('input_files', type=str, nargs='+',
                        help="input xyz files for analysis")
    parser.add_argument('-m', '--mode', choices=['filter', 'cluster', 'labels', 'analyze'],
                        default='cluster')
    parser.add_argument(
        '-a', '--algorithm',
        choices=[*CLUSTERING_ALGORITHMS, 'maxmin'],
        default='auto',
        help="cluster-label algorithm; maxmin uses auto clustering and trims cluster minima",
    )
    parser.add_argument(
        '-n', '--maximum-number-of-seeds',
        type=int,
        default=12,
        help="maximum number of geometries to keep",
    )
    parser.add_argument('--feature', choices=FEATURES, default='auto')
    parser.add_argument(
        '--system-type', choices=SYSTEM_TYPES, default='auto',
        help='System class used by automatic feature selection; auto infers cautiously from XYZ geometry',
    )
    parser.add_argument('--distance', choices=DISTANCE_METRICS, default='euclidean')
    parser.add_argument('--distance-atom-mode', choices=['heavy', 'all'],
                        help='RMSD atoms: graph RMSD defaults to heavy; fragment RMSD to all')
    parser.add_argument('--soap-cutoff', type=float, help='Local SOAP cutoff in Angstrom (default: 5)')
    parser.add_argument('--rematch-alpha', type=float, help='Positive REMatch entropy parameter (default: 1)')
    parser.add_argument('--maximum-mappings', type=int, help='Bounded RMSD matching budget (default: 10000)')
    parser.add_argument('--min-samples', type=int, default=2)
    parser.add_argument('--min-cluster-size', type=int, default=2)
    parser.add_argument('--eps', type=float, help='DBSCAN radius in the units of the distance actually used')
    parser.add_argument('--xi', type=float, default=0.05, help='OPTICS steepness parameter')
    parser.add_argument('--labels-output', help='Write structure labels to this CSV file (labels mode)')
    parser.add_argument('--report-output', help='Write clustering labels and fallback provenance as JSON')
    parser.add_argument(
        '--coordinate-model',
        choices=['covalent-radii', 'distance-cutoff', 'none'],
        default='covalent-radii',
        help='Adjacency model for analyze mode; uses coordinates only and assigns no bond orders or charges',
    )
    parser.add_argument(
        '--bond-scale', type=float, default=1.15,
        help='Covalent-radius sum multiplier for coordinate-only adjacency (analyze mode)',
    )
    parser.add_argument(
        '--bond-cutoff', type=float,
        help='Explicit adjacency cutoff in Angstrom for --coordinate-model distance-cutoff',
    )
    parser.add_argument('--structure-report', help='Write coordinate-only structure analysis as JSON')
    return parser


def main(argv=None, *, prog=None):
    """Cluster or filter the provided XYZ pool and print the selected files."""
    parser = build_parser(prog=prog)
    args = parser.parse_args(argv)
    if args.mode == 'labels' and args.algorithm == 'maxmin':
        parser.error('maxmin selects a subset and does not assign cluster labels')
    if args.report_output and args.mode not in {'cluster', 'labels'}:
        parser.error('--report-output is available in cluster and labels modes')
    if args.labels_output and args.mode != 'labels':
        parser.error('--labels-output is available in labels mode')
    if args.structure_report and args.mode != 'analyze':
        parser.error('--structure-report is available in analyze mode')
    algorithm_options = {
        'min_samples': args.min_samples,
        'min_cluster_size': args.min_cluster_size,
        'eps': args.eps,
        'xi': args.xi,
    }
    distance_options = {key: value for key, value in {
        'atom_mode': args.distance_atom_mode, 'soap_cutoff': args.soap_cutoff,
        'rematch_alpha': args.rematch_alpha, 'max_mappings': args.maximum_mappings,
    }.items() if value is not None}
    if args.mode in {'cluster', 'labels'}:
        try:
            validate_clustering_options(
                'auto' if args.algorithm == 'maxmin' else args.algorithm,
                args.maximum_number_of_seeds, args.distance, algorithm_options,
            )
            from pyar.selection.structural_distances import validate_distance_options
            validate_distance_options(distance_options)
        except ValueError as exc:
            parser.error(str(exc))
    input_files = args.input_files
    if (
        args.mode == 'analyze'
        and args.coordinate_model == 'distance-cutoff'
        and args.bond_cutoff is None
    ):
        parser.error('--coordinate-model distance-cutoff requires --bond-cutoff')
    if (
        args.mode == 'analyze'
        and args.coordinate_model != 'distance-cutoff'
        and args.bond_cutoff is not None
    ):
        parser.error('--bond-cutoff requires --coordinate-model distance-cutoff')
    if len(input_files) < 2 and args.mode not in {'filter', 'analyze'}:
        parser.error('cluster and labels modes require at least two input files')

    mols = []
    for each_file in input_files:
        mol = Molecule.from_xyz(each_file)
        mol.energy = selection_reports.read_energy_from_xyz_file(each_file)
        mol.relative_path = each_file
        mols.append(mol)

    if args.mode == 'analyze':
        try:
            classification = classify_system_pool(
                mols, args.system_type, coordinate_model=args.coordinate_model,
                bond_scale=args.bond_scale, bond_cutoff=args.bond_cutoff,
            )
            report = {
                'method': 'coordinate-only-adjacency',
                'input_files': input_files,
                'classification': classification.to_dict(),
                'resolved_policy': resolve_clustering_policy(
                    classification.system_type, feature=args.feature, algorithm=args.algorithm
                ),
                'structures': [
                    analyze_coordinate_structure(
                        molecule,
                        model=args.coordinate_model,
                        scale=args.bond_scale,
                        cutoff=args.bond_cutoff,
                    )
                    for molecule in mols
                ],
                'limitations': [
                    'XYZ coordinates do not provide bond orders, charge, spin, or fragment intent.',
                    'Adjacency is geometric and should not be interpreted as definitive chemical identity.',
                ],
            }
        except (KeyError, ValueError) as exc:
            parser.error(str(exc))
        if args.structure_report:
            with open(args.structure_report, 'w', encoding='utf-8') as stream:
                json.dump(report, stream, indent=2)
                stream.write('\n')
        else:
            json.dump(report, sys.stdout, indent=2)
            sys.stdout.write('\n')
        return

    selection_reports.print_energy_table(
        mols,
        stream=sys.stdout,
        title="Input pool energies:",
    )
    selected = []
    cluster_result = None
    report = {}
    if args.mode == 'cluster':
        selected = clustering.choose_geometries(
            mols,
            maximum_number_of_seeds=args.maximum_number_of_seeds,
            algorithm=args.algorithm,
            feature=args.feature,
            distance_metric=args.distance,
            algorithm_options=algorithm_options,
            system_type=args.system_type,
            diagnostics=report,
            distance_options=distance_options,
        )
    if args.mode == 'labels':
        cluster_result = cluster_molecules(
            mols,
            feature=args.feature,
            distance_metric=args.distance,
            algorithm=args.algorithm,
            maximum_number_of_clusters=args.maximum_number_of_seeds,
            algorithm_options=algorithm_options,
            system_type=args.system_type,
            distance_options=distance_options,
        )
        selected = mols
        if args.labels_output:
            with open(args.labels_output, 'w', newline='', encoding='utf-8') as stream:
                writer = csv.writer(stream)
                writer.writerow(['input', 'name', 'energy', 'label'])
                for path, molecule, label in zip(input_files, mols, cluster_result.labels):
                    writer.writerow([path, molecule.name, molecule.energy, int(label)])
    if args.mode == 'filter':
        selected = clustering.remove_similar(mols)
    if cluster_result is not None:
        report = cluster_result.to_dict()
    if args.report_output:
        report['input_files'] = input_files
        with open(args.report_output, 'w', encoding='utf-8') as stream:
            json.dump(report, stream, indent=2)
            stream.write('\n')
    selection_reports.print_energy_table(
        selected,
        stream=sys.stdout,
        title="Selected pool energies:",
    )
    print(' '.join(one.name + '.xyz' for one in selected))


if __name__ == '__main__':
    main()
