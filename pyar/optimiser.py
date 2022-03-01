# encoding: utf-8
"""
Optimiser Module

Functions
---------

optimise(molecule, qc_params)
write_csv(csv_file, energy_dit)
bulk_optimize(input_files, qc_params)
"""

import logging
import os

from pyar import file_manager
from pyar.Molecule import Molecule

optimiser_logger = logging.getLogger('pyar.optimiser')


def optimise(molecule, qc_params):
    opt_options = {
        option: qc_params[option]
        for option in ['gamma', 'opt_cycles', 'opt_threshold']
    }

    gamma = qc_params['gamma']
    custom_keyword = qc_params['custom_keyword']

    cwd = os.getcwd()
    if molecule.name == '':
        molecule.name = 'Opt job'
    job_dir = 'job_' + molecule.name
    if not os.path.exists(job_dir):
        file_manager.make_directories(job_dir)
    os.chdir(job_dir)

    if os.path.exists(f'result_{molecule.name}.xyz'):
        read_molecule = Molecule.from_xyz(f'result_{molecule.name}.xyz')
        molecule.energy = read_molecule.energy
        molecule.optimized_coordinates = read_molecule.coordinates
        optimiser_logger.info(
            f'     {molecule.name:35s}: {molecule.energy:15.6f}')
        os.chdir(cwd)
        return True

    software = qc_params['software']
    if software == 'xtb':
        from pyar.interface import xtb
        geometry = xtb.Xtb(molecule, qc_params)
    elif software == 'xtb_turbo':
        if gamma == 0.0:
            from pyar.interface import xtb
            geometry = xtb.Xtb(molecule, qc_params)
        else:
            from pyar.interface import xtbturbo
            geometry = xtbturbo.XtbTurbo(molecule, qc_params)
    elif software == 'turbomole':
        from pyar.interface import turbomole
        geometry = turbomole.Turbomole(molecule, qc_params)
    elif software == "mopac":
        from pyar.interface import mopac
        geometry = mopac.Mopac(molecule, qc_params)
    elif software == "orca":
        from pyar.interface import orca
        geometry = orca.Orca(molecule, qc_params, custom_keyword=custom_keyword)
    elif software == 'obabel':
        from pyar.interface import babel
        geometry = babel.OBabel(molecule)
    elif software == 'psi4':
        from pyar.interface import psi4
        geometry = psi4.Psi4(molecule, qc_params)
    elif software == 'gaussian':
        from pyar.interface import gaussian
        geometry = gaussian.Gaussian(molecule, qc_params)
    else:
        optimiser_logger.error(software, "is not implemented yet")
        return NotImplementedError

    optimize_status = geometry.optimize(opt_options)

    if optimize_status is True \
            or optimize_status == 'converged' \
            or optimize_status == 'CycleExceeded':
        molecule.energy = geometry.energy
        molecule.coordinates = geometry.optimized_coordinates
        optimiser_logger.info(f'     {molecule.name:35s}: {geometry.energy:15.6f}')
    elif optimize_status == 'SCFFailed':
        from numpy.random import uniform
        molecule.coordinates += uniform(-0.1, 0.1, (molecule.number_of_atoms, 3))
        os.chdir(cwd)
        optimize_status = optimise(molecule, qc_params)
    else:
        molecule.energy = None
        molecule.coordinates = None

    os.chdir(cwd)
    return optimize_status


def write_csv_file(csv_filename, energy_dict):
    import csv
    with open(csv_filename, 'w') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(["Name", "Energy"])
        writer.writerows(energy_dict.items())


def bulk_optimize(input_molecules, qc_params):
    status_list = [optimise(each_mol, qc_params) for each_mol in input_molecules]
    return [
        n
        for n, s in zip(input_molecules, status_list)
        if s is True or s == 'cycleexceeded'
    ]


def main():
    pass


if __name__ == '__main__':
    main()
