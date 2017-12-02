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
from afir import restraints
import interface.babel



class XtbTurbo(object):
    def __init__(self, molecule, charge=0, multiplicity=1, scftype='rhf'):
        self.job_name = molecule.name
        self.charge = charge
        self.multiplicity = multiplicity
        self.scftype = scftype
        self.atoms_list = molecule.atoms_list
        self.start_xyz_file = 'trial_' + self.job_name + '.xyz'
        self.xyz_to_coord()
        self.result_xyz_file = 'result_' + self.job_name + '.xyz'
        self.number_of_atoms = molecule.number_of_atoms
        self.job_dir = '{}/job_{}'.format(os.getcwd(), self.job_name)
        self.coord_file = '{}/coord'.format(self.job_dir)
        self.energy_file = '{}/energy'.format(self.job_dir)
        self.prepare_control()

    def prepare_control(self, charge=0, multiplicity=1,
                        read_define_input="no", define_input_file="define.inp"):
        """This function will prepare the control file. Most of the arguments have
           default values(all have in its current form). The read_define_input var
           if "yes", then define_input_file should provide. The define will run
           from this input file."""
        define_log_file = "define.log"
        if read_define_input == "no":
            define_input_file = "define.inp"
            with open(define_input_file, "w") as fdefine:
                fdefine.write("\n\n\n\n")
                fdefine.write("a coord\n*\nno\n")
                fdefine.write("*\n eht\n y \n")
                fdefine.write("{:d}\n".format(charge))
                if multiplicity != 1:
                    return NotImplemented
                fdefine.write("y\n\n")
                fdefine.write("\n\n\n")
                fdefine.write("\n")
                fdefine.write("\n")
                fdefine.write("\n*")
        with open('tmp_bash_command', 'w') as f:
            f.write("define<%s>&define.log" % define_input_file)
        with open(define_log_file, 'w') as fdout:
            subp.check_call(["bash", "tmp_bash_command"], stdout=fdout)
        os.remove('tmp_bash_command')
        if self.check_status_from_log('define.log', 'define') is False:
            print("Error in define.")
            sys.exit()
        return

    @staticmethod
    def run_turbomole_convgrep(outfile):
        with open(outfile, 'a') as fp:
            run_status = subp.Popen(["convgrep"], stdout=fp, stderr=fp)
            out, error = run_status.communicate()
            poll = run_status.poll()
            convgrep_exit_status = run_status.returncode
        return convgrep_exit_status

    @staticmethod
    def check_scf_convergance():
        dscf_conv_status = True
        if os.path.isfile('.sccnotconverged'):
            dscf_conv_status = False
            print("\n :( !SCF Failure! :(")
            print("Check files in", os.getcwd())
        return dscf_conv_status

    @staticmethod
    def check_status_from_log(log_file, program):
        run_status = True
        if program + " ended abnormally" in open(log_file).read():
            run_status = False
        return run_status

    @staticmethod
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

    def run_xtb(self, program, outfile):
        with open(outfile, 'a') as fp:
            if ' ' in program:
                program_object = subp.Popen(program.split(), stdout=fp, stderr=fp)
            else:
                program_object = subp.Popen([program], stdout=fp, stderr=fp)
            program_object.communicate()
            program_object.poll()
            if self.check_status_from_log(outfile, program) is False:
                return False
            return program_object.returncode

    def xyz_to_coord(self, inp_xyz_file=None, outfile='coord'):
        """
        runs turbomole program x2t to create turbomole coord file.
        :param inp_xyz_file: inp_xyz_file will be converted to turbomole coord file.
        :param outfile: coordinate file in turbomole format
        """
        if inp_xyz_file is None:
            # print "xyz file will be searched with job_name"
            xyz_file_name = self.start_xyz_file
        else:
            xyz_file_name = inp_xyz_file
        if not os.path.isfile(xyz_file_name):
            print(xyz_file_name, "not found")
            sys.exit()

        with open(outfile, 'w') as fcoord:
            try:
                subp.check_call(["x2t", str(xyz_file_name), ">"], stdout=fcoord)
            except OSError:
                print("\nx2t not found.\nCheck your turbomole installation")
                sys.exit()


    @property
    def energy(self):
        return float(open(self.energy_file).readlines()[-2].split()[1])

    @property
    def optimized_coordinates(self):
        return np.loadtxt(self.coord_file, comments='$', usecols=(0, 1, 2)) * 0.52917726

    def optimize(self, cycle=10000, gamma=0.0):
        cwd = os.getcwd()
        file_manager.make_directories(self.job_dir)
        for f in ['auxbasis', 'basis', 'control', 'coord', 'mos', 'alpha', 'beta', 'define.inp', 'define.log',
                  'fragment']:
            if os.path.exists(f):
                shutil.move(f, self.job_dir)
        os.chdir(self.job_dir)
        if not os.path.exists('control'):
            print("TURBOMOLE input file (control) not found, stoping")
            sys.exit()
        if gamma > 0.0:
            print("      gamma", gamma)
            if not os.path.isfile('fragment'):
                print("fragment file is not found")
                sys.exit()

        c = 0
        energy_program = "xtb coord"
        gradient_program = "xtb coord -grad"
        update_coord = "statpt"
        status = True
        outfile = 'job.start'
        # initial energy
        self.run_xtb(energy_program, outfile)
        if self.check_scf_convergance() is False:
            print('SF Failure. Check files in', os.getcwd())
            os.chdir(cwd)
            return False
        converged = False
        outfile = "xtb.log"

        while c <= cycle and not converged:
            print("{:4d} {:15.6f}".format(c, self.energy))
            sys.stdout.flush()

            self.run_xtb(gradient_program, outfile)

            if gamma > 0.0:
                status = restraints.isotropic(force=gamma)
                if status is False:
                    print("problem with afir restraints")
                    os.chdir(cwd)
                    return False

            self.run_xtb(update_coord, outfile)
            self.run_turbomole_convgrep(outfile)
            converged = self.check_geometry_convergence(outfile)

            if converged:
                print("\nConverged, at cycle", c)

            self.run_xtb(energy_program, outfile)
            if self.check_scf_convergance() is False:
                print('SF Failure. Check files in', os.getcwd())
                os.chdir(cwd)
                return False
            c += 1
        print()

        if converged:
            interface.babel.write_xyz(self.optimized_coordinates, self.result_xyz_file)
            shutil.copy(self.result_xyz_file, cwd)
            status = True
        elif c > cycle:
            print("cycle exceeded")
            status = 'cycle_exceeded'
        os.chdir(cwd)
        print(os.getcwd())
        sys.stdout.flush()
        return status


def usage():
    print("This usage is for when it will run from the command line.")
    print("otherwise import the module and run prepare_input and optimize.")
    print("\n\n")
    print("USAGE:", sys.argv[0], "start_xyz_file and/or gamma_value")
    print("         gamma_value should be a float")
    print()


def main():
    from Molecule import Molecule
    mol = Molecule.from_xyz(sys.argv[1])
    geometry = XtbTurbo(mol)
    if len(sys.argv) == 2:
        geometry.optimize()
    elif len(sys.argv) == 3:
        gamma_force = sys.argv[2]
        geometry.optimize(gamma=gamma_force)
    else:
        usage()
        sys.exit()
    return


if __name__ == "__main__":
    main()
