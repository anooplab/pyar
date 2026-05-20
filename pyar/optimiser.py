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


def build_geometry(molecule, qc_params):
    """Create the interface wrapper for a configured software backend."""
    software = qc_params['software']
    gamma = qc_params.get('gamma', None)

    if software == 'mlatom_aiqm1':
        from pyar.interface import mlatom_aiqm1
        return mlatom_aiqm1.MlatomAiqm1(molecule, qc_params)
    if software == "gaussian":
        from pyar.interface import gaussian
        return gaussian.Gaussian(molecule, qc_params)
    if software == "orca":
        from pyar.interface import orca
        return orca.Orca(molecule, qc_params)
    if software == "orca-aiqm1":
        from pyar.interface import orca_aiqm1
        geometry = orca_aiqm1.OrcaAIQM1(molecule, qc_params)
        if gamma is not None:
            geometry.set_gamma(gamma)
        return geometry
    if software == "psi4":
        from pyar.interface import psi4
        return psi4.Psi4(molecule, qc_params)
    if software == "xtb":
        from pyar.interface import xtb
        return xtb.Xtb(molecule, qc_params)
    if software == 'xtb_turbo':
        if gamma == 0.0:
            from pyar.interface import xtb
            return xtb.Xtb(molecule, qc_params)
        from pyar.interface import xtbturbo
        return xtbturbo.XtbTurbo(molecule, qc_params)
    if software == 'turbomole':
        from pyar.interface import turbomole
        return turbomole.Turbomole(molecule, qc_params)
    if software == "mopac":
        from pyar.interface import mopac
        return mopac.Mopac(molecule, qc_params)
    if software == "aimnet_2":
        from pyar.interface import aimnet_2
        return aimnet_2.Aimnet2(molecule, qc_params)
    if software == "aiqm1_mlatom":
        from pyar.interface import aiqm1_mlatom
        return aiqm1_mlatom.AIQM1(molecule, qc_params)
    if software == "xtb-aimnet2":
        from pyar.interface import xtb_aimnet2
        return xtb_aimnet2.XtbAimnet2(molecule, qc_params)
    if software == "xtb-aiqm1":
        from pyar.interface import xtb_aiqm1
        return xtb_aiqm1.XtbAIQM1(molecule, qc_params)
    if software == 'obabel':
        from pyar.interface import babel
        return babel.OBabel(molecule)

    raise ValueError(f"Unsupported software backend: {software}")


def optimise(molecule, qc_params):
    cwd = os.getcwd()
    if molecule.name == '':
        molecule.name = 'Opt job'
    job_dir = f'job_{molecule.name}'
    if not os.path.exists(job_dir):
        file_manager.make_directories(job_dir)
    os.chdir(job_dir)
    try:
        if os.path.exists(f'result_{molecule.name}.xyz'):
            read_molecule = Molecule.from_xyz(f'result_{molecule.name}.xyz')
            molecule.energy = read_molecule.energy
            molecule.optimized_coordinates = read_molecule.coordinates
            molecule.coordinates = read_molecule.coordinates
            optimiser_logger.info(f'     {molecule.name:35s}: {molecule.energy:15.6f}')
            return True

        geometry = build_geometry(molecule, qc_params)
        optimize_status = geometry.optimize()
        if is_usable(optimize_status):
            apply_geometry_result(molecule, geometry)
            optimiser_logger.info(f'     {molecule.name:35s}: {float(molecule.energy):15.6f}')
        else:
            molecule.energy = None
            molecule.coordinates = None
        return optimize_status
    finally:
        os.chdir(cwd)


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
