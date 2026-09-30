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
from pyar.selection.clusterers import CLUSTERING_ALGORITHMS, cluster_molecules
from pyar.selection.distances import DISTANCE_METRICS
from pyar.selection.features import FEATURES
from pyar.structure_comparison.coordinate_graph import analyze_coordinate_structure


def main():
    """Cluster or filter the provided XYZ pool and print the selected files."""
    parser = argparse.ArgumentParser()
    parser.add_argument('input_files', type=str, nargs='+',
                        help="input xyz files for analysis")
    parser.add_argument('-m', '--mode', choices=['filter', 'cluster', 'labels', 'analyze'],
                        default='cluster')
    parser.add_argument(
        '-a', '--algorithm',
        choices=[*CLUSTERING_ALGORITHMS, 'maxmin'],
        default='hybrid',
        help="selection algorithm used in cluster mode",
    )
    parser.add_argument(
        '-n', '--maximum-number-of-seeds',
        type=int,
        default=12,
        help="maximum number of geometries to keep",
    )
    parser.add_argument('--feature', choices=FEATURES, default='mbtr')
    parser.add_argument('--distance', choices=DISTANCE_METRICS, default='euclidean')
    parser.add_argument('--min-samples', type=int, default=2)
    parser.add_argument('--min-cluster-size', type=int, default=2)
    parser.add_argument('--eps', type=float, help='DBSCAN radius in standardized feature-vector units')
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
    args = parser.parse_args()
    if args.mode == 'labels' and args.algorithm == 'maxmin':
        parser.error('maxmin selects a subset and does not assign cluster labels')
    if args.report_output and args.mode != 'labels':
        parser.error('--report-output is available in labels mode')
    if args.structure_report and args.mode != 'analyze':
        parser.error('--structure-report is available in analyze mode')
    algorithm_options = {
        'min_samples': args.min_samples,
        'min_cluster_size': args.min_cluster_size,
        'eps': args.eps,
        'xi': args.xi,
    }
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
            report = {
                'method': 'coordinate-only-adjacency',
                'input_files': input_files,
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
    if args.mode == 'cluster':
        selected = clustering.choose_geometries(
            mols,
            maximum_number_of_seeds=args.maximum_number_of_seeds,
            algorithm=args.algorithm,
            feature=args.feature,
            distance_metric=args.distance,
            algorithm_options=algorithm_options,
        )
    if args.mode == 'labels':
        cluster_result = cluster_molecules(
            mols,
            feature=args.feature,
            distance_metric=args.distance,
            algorithm=args.algorithm,
            maximum_number_of_clusters=args.maximum_number_of_seeds,
            algorithm_options=algorithm_options,
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
    if args.report_output and cluster_result is not None:
        with open(args.report_output, 'w', encoding='utf-8') as stream:
            json.dump(cluster_result.to_dict(), stream, indent=2)
            stream.write('\n')
    selection_reports.print_energy_table(
        selected,
        stream=sys.stdout,
        title="Selected pool energies:",
    )
    print(' '.join(one.name + '.xyz' for one in selected))


if __name__ == '__main__':
    main()
