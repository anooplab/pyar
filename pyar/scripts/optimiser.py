#!/usr/bin/env python3
"""Importable `pyar-optimiser` entrypoint."""

import argparse
import datetime
import logging
import os
import sys
from collections import defaultdict

from pyar.Molecule import Molecule
from pyar import optimiser as optimiser_mod
from pyar.data import defualt_parameters

logger = logging.getLogger('pyar')
handler = logging.FileHandler('pyar-optimiser.log', 'w')


def main():
    parser = argparse.ArgumentParser(prog='pyar', description='pyar is a program to predict aggregation, reaction, clustering.')
    parser.add_argument('-v', '--verbosity', choices=[0, 1, 2, 3, 4], type=int)
    parser.add_argument("input_files", metavar='files', type=str, nargs='+')
    reactor_group = parser.add_mutually_exclusive_group(required=False)
    reactor_group.add_argument('-g', '--gamma', type=float)
    reactor_group.add_argument('--site', type=int, nargs=2)
    molecule_group = parser.add_argument_group('molecule', 'Options related to the electronic structure of the molecule')
    molecule_group.add_argument("-c", "--charge", type=int, required=True)
    molecule_group.add_argument("-m", "--multiplicity", type=int, required=True)
    molecule_group.add_argument("--scftype", type=str, choices=['rhf', 'uhf'], default='rhf')
    quantum_chemistry_group = parser.add_argument_group('calculation', 'Calculation specific options')
    quantum_chemistry_group.add_argument("--software", type=str, choices=['gaussian', 'mopac', 'obabel', 'orca', 'psi4', 'turbomole', 'xtb', 'xtb_turbo', 'mlatom_aiqm1', 'aimnet_2', 'aiqm1_mlatom', 'xtb-aimnet2', 'xtb-aiqm1'], required=True)
    quantum_chemistry_group.add_argument('-nprocs', '--nprocs', type=int, nargs=1)
    quantum_chemistry_group.add_argument('-basis', '--basis', type=str)
    quantum_chemistry_group.add_argument('-method', '--method', type=str)
    quantum_chemistry_group.add_argument('--opt-threshold', type=str, default='normal', choices=['loose', 'normal', 'tight'])
    quantum_chemistry_group.add_argument('--opt-cycles', type=int, default=100)
    quantum_chemistry_group.add_argument('--scf-threshold', type=str, default='normal', choices=['loose', 'normal', 'tight'])
    quantum_chemistry_group.add_argument('--scf-cycles', type=int, default=1000)
    quantum_chemistry_group.add_argument('--custom-keywords', type=str)
    args = vars(parser.parse_args())
    run_parameters = defaultdict(lambda: None, defualt_parameters.values)
    for key, value in args.items():
        if value is not None and run_parameters[key] != value:
            run_parameters[key] = value
    logger.setLevel(logging.INFO)
    handler.setFormatter(logging.Formatter('%(message)s'))
    logger.addHandler(handler)
    logger.info(f'Starting pyar-Optimiser at {datetime.datetime.now().strftime("%d %b %Y, %H:%M:%S")}')
    input_molecules = []
    for each_file in run_parameters['input_files']:
        try:
            mol = Molecule.from_xyz(each_file)
            for prop in ['charge', 'multiplicity', 'scftype']:
                vars(mol)[prop] = run_parameters[prop]
            input_molecules.append(mol)
        except IOError:
            logger.critical(f"File {each_file} does not exist")
            sys.exit()
    quantum_chemistry_parameters = {
        'basis': run_parameters['basis'],
        'method': run_parameters['method'],
        'software': run_parameters['software'],
        'opt_cycles': run_parameters['opt_cycles'],
        'opt_threshold': run_parameters['opt_threshold'],
        'scf_cycles': run_parameters['scf_cycles'],
        'scf_threshold': run_parameters['scf_threshold'],
        'nprocs': run_parameters['nprocs'],
        'gamma': run_parameters['gamma'],
        'custom_keywords': run_parameters['custom_keywords'],
        'custom_keyword': run_parameters['custom_keywords']
    }
    try:
        optimiser_mod.bulk_optimize(input_molecules, quantum_chemistry_parameters)
    except FileNotFoundError as exc:
        logger.critical(str(exc))
        sys.exit(str(exc))


if __name__ == '__main__':
    main()
