"""Public selection entrypoints for PyAR 2.0."""

from pyar.selection import reports
from pyar.selection.clustering import choose_geometries
from pyar.selection.clusterers import ClusteringResult, cluster_molecules
from pyar.selection.features import FEATURES, FeatureMatrix, compute_feature_matrix
from pyar.selection.distances import DISTANCE_METRICS, pairwise_distances
from pyar.selection.reports import print_energy_table, read_energy_from_xyz_file

__all__ = [
    "choose_geometries",
    "cluster_molecules",
    "ClusteringResult",
    "FEATURES",
    "FeatureMatrix",
    "compute_feature_matrix",
    "DISTANCE_METRICS",
    "pairwise_distances",
    "reports",
    "print_energy_table",
    "read_energy_from_xyz_file",
]
