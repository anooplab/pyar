"""Conformer CLI options and request translation."""
from __future__ import annotations
import argparse

def build_parser(prog=None):
    """Parse conformer-search CLI arguments."""
    parser = argparse.ArgumentParser(prog=prog or 'pyar conformer', description='Generate and optionally refine RDKit conformers.')
    parser.add_argument('input', help='SMILES string, SDF/MOL file, or XYZ file')
    parser.add_argument('--input-format', choices=['auto', 'smiles', 'sdf', 'mol', 'xyz'], default='auto', help='Input format; auto treats existing .xyz/.sdf/.mol paths as files and other input as SMILES.')
    parser.add_argument('--num-conformers', type=int, default=150)
    parser.add_argument('--top-n', type=int, default=10)
    parser.add_argument('--backend-top-n', type=int)
    advanced = parser.add_argument_group('Advanced generation controls')
    advanced.add_argument('--num-seeds', type=int, default=5)
    advanced.add_argument('--diversity-fraction', type=float, default=0.2)
    advanced.add_argument('--compactness-fraction', type=float, default=0.2, help='Protected contact-rich folded-basin quota; matched by an open-basin quota for diversity.')
    advanced.add_argument('--rms-threshold', '--prune-rms-threshold', dest='rms_threshold', type=float, default=0.25, help='RDKit greedy prune RMS threshold; lower values keep more embedded conformers.')
    advanced.add_argument('--use-random-coords', dest='use_random_coords', action=argparse.BooleanOptionalAction, default=True, help='Start RDKit embedding from random coordinates instead of distance geometry eigenvectors.')
    advanced.add_argument('--torsion-kicks', dest='torsion_kicks', action=argparse.BooleanOptionalAction, default=True, help='Generate local torsion-perturbed conformers before backend refinement.')
    advanced.add_argument('--torsion-mode', choices=['random'], default='random', help='Use the stratified random torsion-kick sampler.')
    advanced.add_argument('--torsion-rounds', type=int, default=2)
    advanced.add_argument('--torsion-kicks-per-conformer', type=int, default=6)
    advanced.add_argument('--torsion-max-bonds', type=int, default=3)
    advanced.add_argument('--torsion-dedup-rms', type=float, default=0.5)
    advanced.add_argument('--dedup-atom-mode', choices=['heavy', 'all'], default='heavy', help='Atoms scored by shared graph RMSD; full connectivity is always checked.')
    parser.add_argument('--force-field', choices=['auto', 'mmff', 'uff'], default='auto')
    advanced.add_argument('--seed', type=int, default=1)
    advanced.add_argument('--num-threads', type=int, default=0)
    advanced.add_argument('--max-iterations', type=int, default=200)
    molecule_group = parser.add_argument_group('molecule')
    molecule_group.add_argument('-c', '--charge', type=int)
    molecule_group.add_argument('-m', '--multiplicity', type=int, default=None)
    molecule_group.add_argument('--scftype', type=str, default=None)
    from pyar.cli_inputs import add_backend_arguments
    add_backend_arguments(parser)
    parser.add_argument('--check', action='store_true', help='Validate without embedding or creating a run directory')
    return parser

def workflow_options(args, qc_params):
    """Translate canonical conformer options into a workflow request."""
    names = ('input_format', 'num_conformers', 'top_n', 'backend_top_n', 'num_seeds', 'diversity_fraction', 'compactness_fraction', 'rms_threshold', 'use_random_coords', 'torsion_kicks', 'torsion_mode', 'torsion_rounds', 'torsion_kicks_per_conformer', 'torsion_max_bonds', 'torsion_dedup_rms', 'dedup_atom_mode', 'force_field', 'seed', 'num_threads', 'max_iterations', 'charge', 'multiplicity', 'scftype')
    return {**{name: getattr(args, name) for name in names}, 'qc_params': qc_params}
