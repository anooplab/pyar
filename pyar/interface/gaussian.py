"""
gaussian.py - interface to gaussian program

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
import subprocess as subp

import numpy as np


class gaussian(object):
    def __init__(self, molecule, charge=0, multiplicity=1, scftype='rhf'):
        self.job_name = molecule.name
        self.charge = charge
        self.multiplicity = multiplicity
        self.scftype = scftype
        self.start_xyz_file = 'trial_' + self.job_name + '.xyz'
        self.result_xyz_file = 'result_' + self.job_name + '.xyz'
        self.inp_file = 'trial_' + self.job_name + '.com'
        self.out_file = 'trial_' + self.job_name + '.log'
        self.atoms_list = molecule.atoms_list
        self.start_coords = molecule.coordinates
        self.optimized_coordinates = []
        self.number_of_atoms = len(self.atoms_list)
        self.energy = 0.0
        keyword="%nprocshared=5\n%chk=molecule.chk\n%mem=5GB\n# opt=(maxcycles=1000) mp2/3-21g"
       
        if sum(molecule.atomic_number)%2 != 0:
            self.multiplicity = 2   
        self.prepare_input(keyword=keyword)


    def prepare_input(self,keyword=""):
        coords=self.start_coords
        f1=open(self.inp_file,"w")
        f1.write(keyword+"\n\n")
        f1.write(self.job_name+"\n\n")
        f1.write(str(self.charge)+" "+str(self.multiplicity)+"\n")
        for i in range(self.number_of_atoms):
            f1.write("%3s  %10.7f  %10.7f %10.7f\n" % (self.atoms_list[i], coords[i][0], coords[i][1], coords[i][2]))
        f1.write("\n")
        f1.close()

    def optimize(self, max_cycles=350, gamma=0.0, restart=False, convergence='normal'):
        """
        :return:This object will return the optimization status. It will
        optimize a structure.
        """
        # TODO: Add a return 'CycleExceeded'
        logfile = "trial_{}.out".format(self.job_name)

        with open(self.out_file, 'w') as fopt:
            out = subp.Popen(["g09", self.inp_file], stdout=fopt, stderr=fopt)
        out.communicate()
        out.poll()
        exit_status = out.returncode
        if exit_status == 0:
            f=open(self.out_file,"r")
            l=f.readlines()
            check_1=0
            check_2=0

            for j in l:
                if "Optimization completed" in j:
                    check_1=1
                if "SCF Done" in j:
                    check_2=1

            if ("Normal termination of Gaussian 09" in l[-1]) and check_1==1 and check_2==1:
                self.energy = self.get_energy()
                self.optimized_coordinates = self.get_coords()
                interface.write_xyz(self.atoms_list, self.optimized_coordinates, self.result_xyz_file, self.job_name, energy=self.energy)
                f.close()
                return True
            else:
                print("Error: OPTIMIZATION PROBABLY FAILED.")
                print("Location: {}".format(os.getcwd()))
                return False

    def get_coords(self):
        """
        :return: coords It will return coordinates
        """
        opt_status = False
        coordinates = []
        with open(self.out_file) as v:
            t = v.readlines()
        for i, lines in enumerate(t):
            if 'Stationary point found.' in lines:
                opt_status = True
            if opt_status is True and 'Standard orientation' in lines:
                pos = i
                coords_lines = t[pos+5:pos+5+self.number_of_atoms]
                for ilines in coords_lines:
                    coordinates.append(ilines.split()[3:6])
                return np.array(coordinates, dtype=float)
        if opt_status is False:
            return None

    def get_energy(self):
        """
        :return:This object will return energy from an orca calculation. It will return Hartree units.
        """
        try:
            with open(self.out_file,"r") as out:
                l = out.readlines()
                en_steps=[]
                for i in range(len(l)):
                    if "SCF Done" in l[i]:
                        en_steps.append( l[i])
                en_Eh=float((en_steps[-1].strip().split())[4])
            return en_Eh
        except IOError:
            print("Warning: File ", self.out_file, "was not found.")


def main():
    pass


if __name__ == "__main__":
    main()
