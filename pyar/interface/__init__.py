import os


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


class SF(object):
    def __init__(self, molecule):
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
