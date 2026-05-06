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

SUCCESS_STATUSES = {True, 'converged'}
CYCLE_EXCEEDED_STATUSES = {'CycleExceeded', 'cycle_exceeded', 'cycleexceeded'}


def is_success(status):
    return status in SUCCESS_STATUSES


def is_cycle_exceeded(status):
    return status in CYCLE_EXCEEDED_STATUSES


def is_usable(status):
    return is_success(status) or is_cycle_exceeded(status)


def apply_geometry_result(molecule, geometry):
    molecule.energy = geometry.energy
    if hasattr(geometry, 'optimized_coordinates'):
        molecule.optimized_coordinates = geometry.optimized_coordinates
        molecule.coordinates = geometry.optimized_coordinates


def optimise(molecule, qc_params):
    cwd = os.getcwd()
    if molecule.name == '':
        molecule.name = 'Opt job'
    job_dir = f'job_{molecule.name}'
    if not os.path.exists(job_dir):
        file_manager.make_directories(job_dir)
    os.chdir(job_dir)
    if os.path.exists(f'result_{molecule.name}.xyz'):
        read_molecule = Molecule.from_xyz(f'result_{molecule.name}.xyz')
        molecule.energy = read_molecule.energy
        molecule.optimized_coordinates = read_molecule.coordinates
        molecule.coordinates = read_molecule.coordinates
        optimiser_logger.info(f'     {molecule.name:35s}: {molecule.energy:15.6f}')
        os.chdir(cwd)
        return True
    software = qc_params['software']
    gamma = qc_params.get('gamma', None)
    print(software)
    if software == 'mlatom_aiqm1':
        from pyar.interface import mlatom_aiqm1
        geometry = mlatom_aiqm1.MlatomAiqm1(molecule, qc_params)
        
    elif software == "orca":
        from pyar.interface import orca
        geometry = orca.Orca(molecule, qc_params)
        
    elif software == "orca-aiqm1":
        from pyar.interface import orca_aiqm1
        geometry = orca_aiqm1.OrcaAIQM1(molecule, qc_params)
        if gamma is not None:
            geometry.set_gamma(gamma)
        
    elif software == "xtb":
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
    elif software == "aimnet_2":
        from pyar.interface import aimnet_2
        geometry = aimnet_2.Aimnet2(molecule, qc_params)
    elif software == "aiqm1_mlatom":
         from pyar.interface import aiqm1_mlatom
         geometry = aiqm1_mlatom.AIQM1(molecule, qc_params) # noqa: F401
    elif software == "xtb-aimnet2":
        from pyar.interface import xtb_aimnet2
        geometry = xtb_aimnet2.XtbAimnet2(molecule, qc_params)
    elif software == "xtb-aiqm1":
        from pyar.interface import xtb_aiqm1
        geometry = xtb_aiqm1.XtbAIQM1(molecule, qc_params)
    elif software == 'obabel': 
        from pyar.interface import babel
        geometry = babel.OBabel(molecule)
    
    optimize_status = geometry.optimize()
    if is_usable(optimize_status):
        apply_geometry_result(molecule, geometry)
        optimiser_logger.info(f'     {molecule.name:35s}: {float(molecule.energy):15.6f}')
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
        if is_usable(s)
    ]


def main():
    pass


if __name__ == '__main__':
    main()
