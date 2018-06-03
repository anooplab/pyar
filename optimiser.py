import os

import file_manager


def optimise(molecule, method, gamma=0.0, custom_keyword=None):
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
    elif software == 'OBabel':
        from interface import babel
        geometry = babel.OBabel(molecule)
    else:
        print(software, "is not implemented yet")
        return NotImplementedError

    optimize_status = geometry.optimize(gamma=gamma)
    if optimize_status is True or optimize_status == 'converged' or optimize_status == 'cycle_exceeded':
        molecule.energy = geometry.energy
        molecule.coordinates = geometry.optimized_coordinates
    else:
        molecule.energy = None
        molecule.coordinates = None
 
    os.chdir(cwd)
    return optimize_status
    return True


def bulk_optimize(input_molecules, method_args, gamma):
    status_list = [optimise(each_mol, method_args, gamma=gamma) for each_mol in input_molecules]
    energy_dict = {n.name: n.energy for n, s in zip(input_molecules, status_list) if s}
    write_csv_file('energy.csv', energy_dict)
    optimized_molecules = [n for n, s in zip(input_molecules, status_list) if s is True or s == 'cycleexceeded']
    return optimized_molecules


def write_csv_file(csv_filename, energy_dict):
    import csv
    with open(csv_filename, 'w') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(["Name", "Energy"])
        writer.writerows(energy_dict.items())


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
