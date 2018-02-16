"""
turbomole.py - interface to turbomole program

Copyright (C) 2016 by Surajit Nandi, Anoop Ayyappan, and Mark P. Waller
Indian Institute of Technology Kharagpur, India and Westfaelische Wilhelms
Universitaet Muenster, Germany

This file is part of the PyAR project.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
"""
import os
import shutil
import subprocess as subp
import sys

import fortranformat as ff
import numpy as np

import file_manager
import interface.babel
from afir import restraints
from units import angstrom2bohr, bohr2angstrom


class Turbomole(object):
    def __init__(self, molecule, charge=0, multiplicity=1, scftype='rhf'):
        self.job_name = molecule.name
        self.charge = charge
        self.multiplicity = multiplicity
        self.scftype = scftype
        self.start_xyz_file = 'trial_' + self.job_name + '.xyz'
        self.result_xyz_file = 'result_' + self.job_name + '.xyz'
        self.number_of_atoms = molecule.number_of_atoms
        self.title = molecule.title
        self.atoms_in_fragments = molecule.fragments
        self.atoms_list = molecule.atoms_list
        self.start_coords = angstrom2bohr(molecule.coordinates)
        self.energy = None
        self.optimized_coordinates = None

    def optimize(self, max_cycles=35, gamma=0.0):

        cwd = os.getcwd()
        job_dir = 'job_' + self.job_name
        file_manager.make_directories(job_dir)

        if gamma > 0.0:
            if os.path.isfile('fragment'):
                shutil.copy('fragment', job_dir)
            else:
                print("file, {}, not found".format('fragment'))
                sys.exit()

        os.chdir(job_dir)

        make_coord(self.atoms_list, self.start_coords)
        prepare_control()

        for cycle in range(max_cycles):
            # Calculate energy
            status, message, energy = self.calc_energy
            if status is False:
                print('Energy evaluation failed')
                os.chdir(cwd)
                return False

            # Calculate Gradients
            status, message, gradients = self.calc_gradients
            if status is False:
                print('Gradient evaluation failed')
                os.chdir(cwd)
                return False

            # Calculate afir gradient if gamma is greater than zero
            if gamma > 0.0:
                re, trg, rg = restraints.isotropic(self.atoms_in_fragments,self.atoms_list,get_coords(), gamma)
                rewrite_turbomole_energy_and_gradient_files(self.number_of_atoms, rg, re, trg)

            # Update coordinates and check convergence.
            status, msg = update_coord()
            if status is True:
                print('converged at', cycle)
                self.energy = get_energy()
                self.optimized_coordinates = bohr2angstrom(get_coords())
                interface.babel.write_xyz(self.atoms_list, self.optimized_coordinates, self.result_xyz_file,
                                          self.job_name,
                                          energy=self.energy)
                shutil.copy(self.result_xyz_file, cwd)
                os.chdir(cwd)
                return True
            with open('energy.dat','a') as fe:
                if gamma > 0.0:
                    fe.writelines("{:3d} {:15.8f} {:15.8f}\n".format(cycle, energy, energy+re))
                else:
                    fe.writelines("{:3d} {:15.8f}\n".format(cycle, energy))
        else:
            print("         cycle exceeded")
            status = 'cycle_exceeded'
            self.energy = get_energy()
            self.optimized_coordinates = bohr2angstrom(get_coords())
            os.chdir(cwd)
            return status

    @property
    def calc_energy(self):
        with open('job.last', 'a') as fp, open('ridft.out', 'w') as fc:
            subp.check_call(['ridft'], stdout=fp, stderr=fc)
        msg = [line for line in open('ridft.out').readlines() if 'ended' in line]
        if os.path.isfile('dscf_problem'):
            msg = "SCF Failure. Check files in"+os.getcwd()
            return False, msg, None
        if 'abnormally' in msg:
            return False, msg, None
        else:
            return True, msg, get_energy()

    @property
    def calc_gradients(self):
        with open('job.last','a') as fp, open('rdgrad.out', 'w') as fc:
            subp.check_call(['rdgrad'], stdout=fp, stderr=fc)
        msg = [line for line in open('rdgrad.out').readlines() if 'ended' in line]
        if 'abnormally' in msg:
            return False, msg, None
        else:
            gradients = get_gradients(self.number_of_atoms)
            return True, msg, gradients


def make_coord(atoms_list, coordinates, outfile='coord'):
    """
    runs turbomole program x2t to create turbomole coord file.
    """

    with open(outfile, 'w') as fcoord:
        fcoord.write('$coord\n')
        for s, c in zip(atoms_list, coordinates):
            fcoord.write("{:20.13} {:20.13} {:20.13} {:3s}\n".format(c[0], c[1], c[2], s))
        fcoord.write('$user-defined bonds\m')
        fcoord.write('$end')
    return


def prepare_control(basis="def2-SVP", func="b-p", ri="on",
                    memory=8000, grid="m4", scf_conv=7, scf_maxiter=500,
                    dftd3="yes", charge=0, multiplicity=1,
                    read_define_input="no", define_input_file="define.inp"):

    """This function will prepare the control file. Most of the arguments have
       default values(all have in its current form). The read_define_input var
       if "yes", then define_input_file should provide. The define will run
       from this input file."""

    # Turbomole 'coord' file should exist
    if not os.path.isfile('coord'):
        print("Turbomole coordinate file, 'coord', should exist")
        sys.exit()

    if read_define_input == "no":
        define_input_file = "define.inp"
        with open(define_input_file, "w") as fdefine:
            fdefine.write("\n\n")
            fdefine.write("a coord\n*\nno\n")
            fdefine.write("bb all %s\n" % basis)
            fdefine.write("*\neht\ny \n")
            fdefine.write("%d\n" % charge)
            fdefine.write("y\n\n")
            fdefine.write("\n\n\n")
            fdefine.write("dft\n")
            fdefine.write("on\n")
            fdefine.write("func %s\n" % func)  # here put condition for advance use
            fdefine.write("grid %s\n" % grid)
            fdefine.write("\n")
            if ri == "on":
                fdefine.write("ri\n on\n m\n")
                fdefine.write("%d\n" % memory)
                fdefine.write("\n")
            if dftd3=="yes":
                fdefine.write("dsp\n on\n")
                fdefine.write("\n")
            fdefine.write("scf\nconv\n")
            fdefine.write("%d\n" % scf_conv)
            fdefine.write("iter\n")
            fdefine.write("%d\n" % scf_maxiter)
            fdefine.write("\n*")

    define_log_file = "define.log"
    with open('tmp_bash_command', 'w') as f:
        f.write("define<%s>&define.log" % define_input_file)
    with open(define_log_file, 'w') as fdout:
        subp.check_call(["bash", "tmp_bash_command"], stdout=fdout)
    os.remove('tmp_bash_command')
    if 'abnormally' in open(define_log_file).read():
        print('define ended abnormally')
        return False
    return True


def get_energy():
    with open('energy') as fp:
        return float(fp.readlines()[-2].split()[1])


def get_gradients(natoms):
    start =natoms+1
    grad = [line.replace('D', 'E').split() for line in open('gradient').readlines()[-start:-1]]
    return grad


def get_coords():
    return np.loadtxt('coord', comments='$', usecols=(0, 1, 2))


def update_coord():
    with open('job.last', 'a') as fp, open('statpt.out', 'w') as fc:
        subp.check_call(['relax'], stdout=fp, stderr=fc)
    msg = [line for line in open('statpt.out').readlines() if 'ended' in line]
    if 'abnormally' in msg:
        return False, msg
    else:
        return check_geometry_convergence('statpt.out'), msg


def check_geometry_convergence(outfile):
    convergence_status = False
    with open("not.converged") as fp:
        all_lines = fp.readlines()
    with open(outfile, 'a') as fout:
        for lines in all_lines:
            if 'convergence reached' in lines:
                fout.write("THE OPTIMIZATION IS CONVERGED.\n")
                convergence_status = True
                os.remove("not.converged")
    return convergence_status


def rewrite_turbomole_gradient(number_of_atoms, total_restraint_energy, total_restraint_gradients,
                               restraint_gradient):
    """Rewrite new gradient file in turbomole format"""
    import re
    with open('gradient', 'r') as gradient_file:
        gradient_contents = gradient_file.read()
        gradients = re.split("cycle", gradient_contents)
        prelog = "cycle".join(gradients[:-1])
        last_gradient = gradients[-1]
        list_from_last_gradient = last_gradient.split('\n')

        cycle_number = int(list_from_last_gradient[0].split()[1])
        scf_energy = float(list_from_last_gradient[0].split()[5])
        total_gradients = float(list_from_last_gradient[0].split(' =')[-1])
        new_energy = scf_energy + total_restraint_energy
        new_total_gradients = total_gradients + total_restraint_gradients
        new_gradients = "cycle = %6d    SCF energy =%20.10f   |dE/dxyz| =%10.6f \n" \
                        % (cycle_number, new_energy, new_total_gradients)

        for line in list_from_last_gradient[1:number_of_atoms + 1]:
            new_gradients = new_gradients + line + "\n"

        i = 0
        for line in list_from_last_gradient[number_of_atoms + 1:number_of_atoms * 2 + 1]:
            dx, dy, dz = line.split()[0], line.split()[1], line.split()[2]
            dx = float(re.sub('D', 'E', dx)) - restraint_gradient[i][0]
            dy = float(re.sub('D', 'E', dy)) - restraint_gradient[i][1]
            dz = float(re.sub('D', 'E', dz)) - restraint_gradient[i][2]
            formatted = ff.FortranRecordWriter('3D22.13')
            temp_line = formatted.write([dx, dy, dz])
            new_gradients = new_gradients + str(temp_line) + "\n"
            i += 1
    with open('gradient', 'w') as g:
        g.write(prelog + new_gradients + '$end')


def rewrite_turbomole_energy(total_restraint_energy):
    import re
    with open('energy', 'r') as energy_file:
        energy = energy_file.readlines()
        till_last_energy = energy[:-2]
        old_energy = energy[-2].split()[1]
        new_energy = float(old_energy) + total_restraint_energy
    with open('energy', 'w') as new_energy_file:
        new_energy_file.write(''.join(till_last_energy) + re.sub(old_energy, str(new_energy), energy[-2]) + '$end\n')


def rewrite_turbomole_energy_and_gradient_files(number_of_atoms,
                                                restraint_gradients,
                                                total_restraint_energy,
                                                total_restraint_gradients):

    rewrite_turbomole_gradient(number_of_atoms, total_restraint_energy,
                               total_restraint_gradients,
                               restraint_gradients)
    rewrite_turbomole_energy(total_restraint_energy)


def main():
    from Molecule import Molecule
    mol = Molecule.from_xyz(sys.argv[1])
    geometry = Turbomole(mol)
    if len(sys.argv) == 2:
        geometry.optimize()
    elif len(sys.argv) == 3:
        gamma_force = sys.argv[2]
        geometry.optimize(gamma=gamma_force)
    else:
        sys.exit()
    return

if __name__ == "__main__":
    from Molecule import Molecule


    mol = Molecule.from_xyz(sys.argv[1])
    geometry = Turbomole(mol)
    geometry.optimize(1000)
