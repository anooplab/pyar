import os

import file_manager

import logging
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
        from interface import xtb
        geometry = xtb.Xtb(molecule, method)
    elif software == 'xtb_turbo':
        if gamma == 0.0:
            from interface import xtb
            geometry = xtb.Xtb(molecule, method)
        else:
            from interface import xtbturbo
            geometry = xtbturbo.XtbTurbo(molecule, method)
    elif software == 'turbomole':
        from interface import turbomole
        geometry = turbomole.Turbomole(molecule, method)
    elif software == "mopac":
        from interface import mopac
        geometry = mopac.Mopac(molecule, method)
    elif software == "orca":
        from interface import orca
        geometry = orca.Orca(molecule, method, custom_keyword=custom_keyword)
    elif software == 'obabel':
        from interface import babel
        geometry = babel.OBabel(molecule)
    elif software == 'psi4':
        from interface import psi4
        geometry = psi4.Psi4(molecule, method)
    else:
        print(software, "is not implemented yet")
        return NotImplementedError

    optimize_status = geometry.optimize(gamma=gamma, max_cycles=max_cycles, convergence=convergence)
    if optimize_status is True\
            or optimize_status == 'converged'\
            or optimize_status == 'CycleExceeded':
        molecule.energy = geometry.energy
        molecule.coordinates = geometry.optimized_coordinates
        optimiser_logger.info("Energy: %15.6f", geometry.energy)
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
    import sys
    input_files = sys.argv[1:]
    from Molecule import Molecule
    method = {'charge': 0, 'multiplicity': 1, 'scftype': 'rhf', 'software': 'orca'}
    gamma = 0.0

    for m in input_files:
        mol = Molecule.from_xyz(m)
        if optimise(mol, method=method, gamma=gamma):
            print(mol.energy)


if __name__ == '__main__':
    main()