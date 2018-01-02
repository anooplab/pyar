"""
orca.py - interface to mopac program

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
import sys

import numpy as np


class Orca(object):
    def __init__(self, molecule, charge=0, multiplicity=1, scftype='rhf'):
        self.job_name = molecule.name
        self.charge = charge
        if (sum(molecule.atomic_number) - self.charge) % 2 == 1 and multiplicity == 1:
            self.multiplicity = 2
        else:
            self.multiplicity = multiplicity
        if self.multiplicity % 2 == 0 and scftype is 'rhf':
            self.scftype = 'uhf'
        else:
            self.scftype = scftype
        self.start_xyz_file = 'trial_' + self.job_name + '.xyz'
        self.result_xyz_file = 'result_' + self.job_name + '.xyz'
        self.inp_file = 'trial_' + self.job_name + '.inp'
        self.out_file = 'trial_' + self.job_name + '.out'
        self.atoms_list = molecule.atoms_list
        self.start_coords = molecule.coordinates
        self.optimized_coordinates = []
        self.number_of_atoms = len(self.atoms_list)
        self.energy = 0.0
        keyword="!BP RI opt def2-SVP def2-SVP/J KDIIS"
        if any(x >=21 for x in molecule.atomic_number):
            keyword += 'def2-ECP'
        self.prepare_input(keyword=keyword)


    def prepare_input(self,keyword=""):
        coords=self.start_coords
        f1=open(self.inp_file,"w")
        if self.scftype is 'uks':
            keyword += 'UKS'
        f1.write(keyword+"\n")
        f1.write("*xyz {0} {1}\n".format(str(self.charge), str(self.multiplicity)))
        for i in range(self.number_of_atoms):
            f1.write(" "+"%3s  %10.7f  %10.7f %10.7f\n" % (self.atoms_list[i], coords[i][0], coords[i][1], coords[i][2]))
        f1.write("*")
        f1.close()


    def optimize(self, gamma=None):
        """
        :return:This object will return the optimization status. It will
        optimize a structure.
        """
        logfile = "trial_{}.log".format(self.job_name)
	
        with open(self.out_file, 'w') as fopt:
            out = subp.Popen(["orca", self.inp_file], stdout=fopt, stderr=fopt)
        out.communicate()
        out.poll()
        exit_status = out.returncode
        if exit_status == 0:
            f=open(self.out_file,"r")
            l=f.readlines()
            if ("****ORCA TERMINATED NORMALLY****" in l[-2]):
                self.energy = self.get_energy()
                self.optimized_coordinates = np.loadtxt(self.inp_file[:-4]+".xyz", dtype=float, skiprows=2, usecols=(1, 2, 3))
                self.write_xyz(self.optimized_coordinates, self.result_xyz_file)
                f.close()
                return True
            else:
                print("Error: OPTIMIZATION PROBABLY FAILED. CHECK THE .out FILE FOR PARTIAL OPTIMIZTION ")
                print("Check for partial optimization.")
                print("Location: {}".format(os.getcwd()))
                return False

    def get_energy(self):
        """
        :return:This object will return energy from an orca calculation. It will return Hartree units.
        """
        try:
            with open(self.out_file,"r") as out:
                l = out.readlines()
                en_steps=[]
                for i in range(len(l)):
                    if "FINAL SINGLE POINT ENERGY" in l[i]:
                        en_steps.append( l[i])
                en_Eh=float((en_steps[-1].strip().split())[-1])
            return en_Eh
        except IOError:
            print("Warning: File ", self.out_file, "was not found.")

    def write_xyz(self, coords, filename):
        with open(filename, 'w') as fp:
            fp.write(str(self.number_of_atoms) + '\n')
            fp.write(self.job_name + '\n')
            for i in range(self.number_of_atoms):
                fp.write("%3s  %10.7f  %10.7f %10.7f\n"
                         % (self.atoms_list[i], coords[i][0], coords[i][1], coords[i][2]))



def main():
    pass


if __name__ == "__main__":
    main()
