import os
from pyar.mlatom.data import molecule  # noqa: F401


def which(program):
    import os

    def is_exe(exec_path):
        return os.path.isfile(exec_path) and os.access(exec_path, os.X_OK)

    file_path, file_name = os.path.split(program)
    if file_path:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


# from pyar.interface.mlatom_aiqm1 import MlatomAiqm1

class SF(object):
    def __init__(self, molecule):  # noqa: F811
        self.job_name = molecule.name
        self.start_xyz_file = 'trial_' + self.job_name + '.xyz'
        self.result_xyz_file = 'result_' + self.job_name + '.xyz'
        self.atoms_list = molecule.atoms_list
        self.number_of_atoms = molecule.number_of_atoms
        self.title = 'opt {}'.format(molecule.title)
        self.charge = molecule.charge
        self.multiplicity = molecule.multiplicity
        self.scftype = molecule.scftype
        if not os.path.isfile(self.start_xyz_file):
            write_xyz(self.atoms_list,
                      molecule.coordinates, self.start_xyz_file,
                      job_name=self.job_name)


def write_xyz(atoms_list, coordinates, filename, job_name='no_name', energy=0.0):
    with open(filename, 'w') as fp:
        fp.write("%3d\n" % len(coordinates))
        fp.write(job_name + ':' + str(energy) + '\n')
        for a, c in zip(atoms_list, coordinates):
            fp.write("{:<2}{:12.5f}{:12.5f}{:12.5f}\n".format(a, c[0], c[1], c[2]))


# import torchani  # noqa: E402
import ase.optimize  # noqa: E402
from ase.calculators.calculator import Calculator  # noqa: E402


class ANI(Calculator):
    def __init__(self, species, model='ANI-1x'):
        self.species = species
        self.model = torchani.models.__dict__[model]()

    def calculate(self, atoms=None, properties=['energy'], system_changes=[]):
        positions = atoms.get_positions()
        species = atoms.get_chemical_symbols()
        energy = self.model((species, positions)).energy
        self.results = {'energy': energy}

    def optimize(atoms):
        dyn = ase.optimize.BFGS(atoms, trajectory='optimization.traj')
        dyn.run(fmax=0.001)
        return atoms.get_potential_energy()

    def singlepoint(atoms):
        calc = ANI(atoms.get_chemical_symbols())
        atoms.set_calculator(calc)
        energy = atoms.get_potential_energy()
        return energy
