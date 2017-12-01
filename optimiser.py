def optimise(molecule, method, gamma=0.0):
    job_name = molecule.name

    charge = method['charge']
    scf_type = method['scftype']
    multiplicity = method['multiplicity']

    software = method['software']
    if software == 'xtb':
        from interface import xtb
        geometry = xtb.Xtb(molecule, charge, multiplicity, scf_type)
    elif software == 'xtb_turbo':
        from interface import xtbturbo
        geometry = xtbturbo.XtbTurbo(molecule)
    elif software == 'turbomole':
        from interface import turbomole
        geometry = turbomole.Turbomole(molecule)
    elif software == "mopac":
        from interface import mopac
        geometry = mopac.Mopac(molecule)
    elif software == 'obabel':
        from interface import babel
        geometry = babel.obabel(molecule)
    else:
        print(software, "is not implemented yet")
        return NotImplementedError
    optimize_status = geometry.optimize(gamma=gamma)
    if optimize_status:
        molecule.energy = geometry.energy
        molecule.coordinates = geometry.optimized_coordinates
        return optimize_status
    else:
        return optimize_status


def main():
    import sys
    input_files = sys.argv[1:]
    from Molecule import Molecule
    method = {'charge':0, 'multiplicity':1, 'scftype':'rhf', 'software'
    :'xtb'}
    gamma = 0.0

    for m in input_files:
        mol = Molecule.from_xyz(m)
        if optimise(mol, method=method, gamma=gamma): print(mol.energy)


if __name__ == '__main__':
    main()