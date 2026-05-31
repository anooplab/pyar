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
from pyar.core.molecule import Molecule
from pyar.backend_capabilities import backend_supports_geometry_optimization

optimiser_logger = logging.getLogger('pyar.optimiser')

SUCCESS_STATUSES = {True, 'converged'}
CYCLE_EXCEEDED_STATUSES = {'CycleExceeded', 'cycle_exceeded', 'cycleexceeded'}


def is_success(status):
    return status in SUCCESS_STATUSES


def is_cycle_exceeded(status):
    return status in CYCLE_EXCEEDED_STATUSES


def is_usable(status):
    return is_success(status) or is_cycle_exceeded(status)


def _status_label(status):
    """Return a normalized human-readable optimization status label."""
    if is_success(status):
        return "success"
    if is_cycle_exceeded(status):
        return "cycle_exceeded_usable"
    if status is None:
        return "failed_none_status"
    return f"failed:{status}"


def apply_geometry_result(molecule, geometry):
    molecule.energy = geometry.energy
    if hasattr(geometry, 'optimized_coordinates'):
        molecule.optimized_coordinates = geometry.optimized_coordinates
        molecule.coordinates = geometry.optimized_coordinates


def build_geometry(molecule, qc_params):
    """Create the interface wrapper for a configured software backend."""
    software = qc_params['software']
    geometry_optimizer = qc_params.get('geometry_optimizer', 'native')
    gamma = qc_params.get('gamma', None)

    if geometry_optimizer == 'geometric':
        if not backend_supports_geometry_optimization(software):
            raise ValueError(
                f"Backend {software!r} cannot be used with geomeTRIC AFIR optimisation "
                "because it does not expose Cartesian energy and gradients."
            )
        from pyar.backends import geometric
        return geometric.Geometric(molecule, qc_params)

    if software == 'mlatom_aiqm1':
        from pyar.backends import mlatom_aiqm1
        return mlatom_aiqm1.MlatomAiqm1(molecule, qc_params)
    if software == "gaussian":
        from pyar.backends import gaussian
        return gaussian.Gaussian(molecule, qc_params)
    if software == "orca":
        from pyar.backends import orca
        return orca.Orca(molecule, qc_params)
    if software == "orca-aiqm1":
        from pyar.backends import orca_aiqm1
        geometry = orca_aiqm1.OrcaAIQM1(molecule, qc_params)
        if gamma is not None:
            geometry.set_gamma(gamma)
        return geometry
    if software == "psi4":
        from pyar.backends import psi4
        return psi4.Psi4(molecule, qc_params)
    if software == "xtb":
        from pyar.backends import xtb
        return xtb.Xtb(molecule, qc_params)
    if software == 'xtb_turbo':
        if gamma == 0.0:
            from pyar.backends import xtb
            return xtb.Xtb(molecule, qc_params)
        from pyar.backends import xtb_turbo
        return xtb_turbo.XtbTurbo(molecule, qc_params)
    if software == 'turbomole':
        from pyar.backends import turbomole
        return turbomole.Turbomole(molecule, qc_params)
    if software == "mopac":
        from pyar.backends import mopac
        return mopac.Mopac(molecule, qc_params)
    if software == "aimnet_2":
        from pyar.backends import aimnet_2
        return aimnet_2.Aimnet2(molecule, qc_params)
    if software == "aiqm1_mlatom":
        from pyar.backends import aiqm1_mlatom
        return aiqm1_mlatom.AIQM1(molecule, qc_params)
    if software == "xtb-aimnet2":
        from pyar.backends import xtb_aimnet2
        return xtb_aimnet2.XtbAimnet2(molecule, qc_params)
    if software == "xtb-aiqm1":
        from pyar.backends import xtb_aiqm1
        return xtb_aiqm1.XtbAIQM1(molecule, qc_params)
    if software == 'obabel':
        from pyar.backends import babel
        return babel.OBabel(molecule)

    raise ValueError(f"Unsupported software backend: {software}")


def optimise(molecule, qc_params):
    """Run one optimization job and return the backend status object."""
    cwd = os.getcwd()
    if molecule.name == '':
        molecule.name = 'Opt job'
    job_dir = f'job_{molecule.name}'
    if not os.path.exists(job_dir):
        file_manager.make_directories(job_dir)
    os.chdir(job_dir)
    try:
        software = qc_params.get('software', 'unknown')
        gamma = qc_params.get('gamma', None)
        optimiser_logger.debug(
            "Optimization job start: name=%s software=%s geometry_optimizer=%s gamma=%s dir=%s",
            molecule.name,
            software,
            qc_params.get('geometry_optimizer', 'native'),
            gamma,
            os.path.abspath(job_dir),
        )
        if os.path.exists(f'result_{molecule.name}.xyz'):
            read_molecule = Molecule.from_xyz(f'result_{molecule.name}.xyz')
            molecule.energy = read_molecule.energy
            molecule.optimized_coordinates = read_molecule.coordinates
            molecule.coordinates = read_molecule.coordinates
            optimiser_logger.info(
                "Optimization reused cached result: name=%s status=success energy=%15.6f",
                molecule.name,
                float(molecule.energy),
            )
            optimiser_logger.info(f'     {molecule.name:35s}: {molecule.energy:15.6f}')
            return True

        geometry = build_geometry(molecule, qc_params)
        optimize_status = geometry.optimize()
        normalized_status = _status_label(optimize_status)
        if is_usable(optimize_status):
            apply_geometry_result(molecule, geometry)
            optimiser_logger.info(
                "Optimization completed: name=%s status=%s energy=%15.6f",
                molecule.name,
                normalized_status,
                float(molecule.energy),
            )
            optimiser_logger.info(f'     {molecule.name:35s}: {float(molecule.energy):15.6f}')
        else:
            molecule.energy = None
            molecule.coordinates = None
            optimiser_logger.warning(
                "Optimization failed: name=%s status=%s software=%s",
                molecule.name,
                normalized_status,
                software,
            )
        return optimize_status
    except Exception:
        optimiser_logger.exception("Optimization crashed: name=%s", molecule.name)
        molecule.energy = None
        molecule.coordinates = None
        return None
    finally:
        os.chdir(cwd)


def write_csv_file(csv_filename, energy_dict):
    import csv
    with open(csv_filename, 'w') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(["Name", "Energy"])
        writer.writerows(energy_dict.items())


def bulk_optimize(input_molecules, qc_params):
    """Optimize molecules and keep only usable results."""
    status_list = [optimise(each_mol, qc_params) for each_mol in input_molecules]
    status_counts = {}
    for status in status_list:
        key = _status_label(status)
        status_counts[key] = status_counts.get(key, 0) + 1
    optimiser_logger.info(
        "Bulk optimization summary: total=%d usable=%d status_counts=%s",
        len(input_molecules),
        sum(1 for s in status_list if is_usable(s)),
        status_counts,
    )
    return [
        n
        for n, s in zip(input_molecules, status_list)
        if is_usable(s)
    ]


def main():
    pass


if __name__ == '__main__':
    main()
