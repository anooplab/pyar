import logging
import os

from pyar import file_manager

optimiser_logger = logging.getLogger('pyar.optimiser')


def optimise(molecule, method, gamma=0.0, max_cycles=350, convergence='normal', restart='False', custom_keyword=None):
    cwd = os.getcwd()
    if molecule.name == '':
        molecule.name = 'opt'
    job_dir = 'job_' + molecule.name
    file_manager.make_directories(job_dir)
    os.chdir(job_dir)

    software = method['software']
    if software == 'xtb':
        from pyar.interface import xtb
        geometry = xtb.Xtb(molecule, method)
    elif software == 'xtb_turbo':
        if gamma == 0.0:
            from pyar.interface import xtb
            geometry = xtb.Xtb(molecule, method)
        else:
            from pyar.interface import xtbturbo
            geometry = xtbturbo.XtbTurbo(molecule, method)
    elif software == 'turbomole':
        from pyar.interface import turbomole
        geometry = turbomole.Turbomole(molecule, method)
    elif software == "mopac":
        from pyar.interface import mopac
        geometry = mopac.Mopac(molecule, method)
    elif software == "orca":
        from pyar.interface import orca
        geometry = orca.Orca(molecule, method, custom_keyword=custom_keyword)
    elif software == 'obabel':
        from pyar.interface import babel
        geometry = babel.OBabel(molecule)
    elif software == 'psi4':
        from pyar.interface import psi4
        geometry = psi4.Psi4(molecule, method)
    elif software == 'gaussian':
        from pyar.interface import gaussian
        geometry = gaussian.Gaussian(molecule, method)
    else:
        print(software, "is not implemented yet")
        return NotImplementedError

    optimize_status = geometry.optimize(gamma=gamma,
                                        max_cycles=max_cycles,
                                        convergence=convergence)

    if optimize_status is True \
            or optimize_status == 'converged' \
            or optimize_status == 'CycleExceeded':
        molecule.energy = geometry.energy
        molecule.coordinates = geometry.optimized_coordinates
        optimiser_logger.info('Energy ({:9s}): {:15.6f}'.format(molecule.name, geometry.energy))
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


def bulk_optimize(input_molecules, method_args, gamma):
    status_list = [optimise(each_mol, method_args, gamma=gamma) for each_mol in input_molecules]
    optimized_molecules = [n for n, s in zip(input_molecules, status_list) if s is True or s == 'cycleexceeded']
    return optimized_molecules


def main():
    pass


if __name__ == '__main__':
    main()
