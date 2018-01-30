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

import numpy as np

import file_manager
import interface.babel
from afir import restraints


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
        self.atoms_list = molecule.atoms_list
        self.start_coords = molecule.coordinates
        self.prepare_input()
        self.energy = None
        self.optimized_coordinates = None

    def make_coord(self, inp_xyz_file=None, outfile='coord'):
        """
        runs turbomole program x2t to create turbomole coord file.
        :param inp_xyz_file: xyz file containing coordinates
        :param outfile: coordinate file in turbomole format
        """
        if inp_xyz_file is None:
            # print "xyz file will be searched with job_name"
            xyz_file_name = self.start_xyz_file
        else:
            xyz_file_name = inp_xyz_file
        if not os.path.isfile(xyz_file_name):
            interface.babel.write_xyz(self.atoms_list, self.start_coords, self.start_xyz_file, self.job_name)

        with open(outfile, 'w') as fcoord:
            try:
                subp.check_call(["x2t", str(xyz_file_name), ">"], stdout=fcoord)
            except OSError:
                print("\nx2t not found.\nCheck your turbomole installation")
                sys.exit()
        return

    def prepare_control(self, basis="def2-SVP", func="b-p", ri="on",
                        memory=8000, grid="m4", scf_conv=7, scf_maxiter=500,
                        dftd3="yes", charge=0, multiplicity=1,
                        read_define_input="no", define_input_file="define.inp"):
        """This function will prepare the control file. Most of the arguments have
           default values(all have in its current form). The read_define_input var
           if "yes", then define_input_file should provide. The define will run
           from this input file."""
        define_log_file = "define.log"
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
                fdefine.write("dsp\n on\n")
                fdefine.write("\n")
                fdefine.write("scf\n conv\n")
                fdefine.write("%d\n" % scf_conv)
                fdefine.write("iter\n")
                fdefine.write("%d\n" % scf_maxiter)
                fdefine.write("\n*")
        with open('tmp_bash_command', 'w') as f:
            f.write("define<%s>&define.log" % define_input_file)
        with open(define_log_file, 'w') as fdout:
            subp.check_call(["bash", "tmp_bash_command"], stdout=fdout)
        os.remove('tmp_bash_command')
        define_run_status = self.check_status_from_log(define_log_file, 'define')
        return define_run_status

    def prepare_input(self):
        self.make_coord()
        define_status = self.prepare_control()
        if define_status:
            print("Error in define.")
            sys.exit()
        return define_status

    @staticmethod
    def check_status_from_log(log_file, program):
        run_status = 0
        with open(log_file) as fp:
            for this_line in fp.readlines():
                program = ''.join(program)
                if program + " ended abnormally" in this_line:
                    print(this_line)
                    run_status = 1
        return run_status

    @staticmethod
    def run_turbomole_module(program, outfile):
        with open(outfile, 'a') as fp:
            if ' ' in program:
                program_object = subp.Popen(program.split(), stdout=fp, stderr=fp)
            else:
                program_object = subp.Popen([program], stdout=fp, stderr=fp)
            program_object.communicate()
            program_object.poll()
            return program_object.returncode

    @staticmethod
    def run_convgrep(outfile):
        with open(outfile, 'a') as fp:
            run_status = subp.Popen(["convgrep"], stdout=fp, stderr=fp)
            out, error = run_status.communicate()
            poll = run_status.poll()
            convgrep_exit_status = run_status.returncode
        return convgrep_exit_status

    @staticmethod
    def check_dscf():
        dscf_conv_status = 0
        if os.path.isfile('dscf_problem'):
            dscf_conv_status = 1
            print("\n ========= :( !SCF Failure! :( ==========")
            print("Check files in", os.getcwd())
        return dscf_conv_status

    @staticmethod
    def check_convergence(outfile):
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

    def get_energy(self):
        with open('energy') as fp:
            return float(fp.readlines()[-2].split()[1])

    def get_coords(self):
        return np.loadtxt('coord', comments='$', usecols=(0, 1, 2)) * 0.52917726

    def optimize(self, ri="on", cycle=10000, gamma=0.0):
        cwd = os.getcwd()
        job_dir = 'job_' + self.job_name
        file_manager.make_directories(job_dir)
        for f in ['auxbasis', 'basis', 'control', 'coord', 'mos', 'alpha', 'beta', 'define.inp', 'define.log',
                  'fragment']:
            if os.path.exists(f):
                shutil.move(f, job_dir)
        os.chdir(job_dir)
        if not os.path.exists('control'):
            print("TURBOMOLE input file (control) not found, stopping")
            sys.exit()
        if gamma > 0.0:
            print("      gamma", gamma)
            if not os.path.isfile('fragment'):
                print("fragment file is not found")
                sys.exit()

        c = 0
        if ri == "on":
            energy_program = "ridft"
            gradient_program = "rdgrad"
        else:
            energy_program = "dscf"
            gradient_program = "grad"
        update_coord = "statpt"
        status = True
        outfile = 'job.start'
        # initial energy
        if self.run_turbomole_module(energy_program, outfile) \
                or self.check_status_from_log(outfile, energy_program) \
                or self.check_dscf():
            os.chdir(cwd)
            return False
        converged = False

        while c <= cycle and not converged:
            sys.stdout.flush()
            outfile = "job.last"
            if self.run_turbomole_module(gradient_program, outfile) \
                    or self.check_status_from_log(outfile, gradient_program):
                os.chdir(cwd)
                return False

            if gamma > 0.0:
                if restraints.isotropic(force=gamma) is False:
                    print("problem with afir restraints")
                    os.chdir(cwd)
                    return False

            if self.run_turbomole_module(update_coord, outfile):
                print("problem with statpt")
                os.chdir(cwd)
                return False

            self.run_convgrep(outfile)
            converged = self.check_convergence(outfile)

            if converged:
                print("\nConverged, at cycle", c)

            if self.run_turbomole_module(energy_program, outfile) \
                    or self.check_status_from_log(outfile, energy_program) \
                    or self.check_dscf():
                os.chdir(cwd)
                return False
            c += 1
        print()

        if converged:
            self.energy = self.get_energy()
            self.optimized_coordinates = self.get_coords()
            interface.babel.write_xyz(self.atoms_list, self.optimized_coordinates, self.result_xyz_file, self.job_name,
                                      energy=self.energy)
            shutil.copy(self.result_xyz_file, cwd)
            status = True
        elif c > cycle:
            print("cycle exceeded")
            status = 'cycle_exceeded'
        os.chdir(cwd)
        print(os.getcwd())
        sys.stdout.flush()
        return status


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
    main()
