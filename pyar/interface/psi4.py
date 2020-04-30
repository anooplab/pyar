"""
psi4.py - interface to PSI4 program

Copyright (C) 2016 by Anoop Ayyappan
Indian Institute of Technology Kharagpur, India

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

from pyar.interface import SF, write_xyz


class Psi4(SF):
    def __init__(self, molecule, method, custom_keyword=None):

        super(Psi4, self).__init__(molecule)

        self.start_coords = molecule.coordinates
        self.inp_file = 'trial_' + self.job_name + '.in'
        self.out_file = 'trial_' + self.job_name + '.out'
        self.optimized_coordinates = []
        self.number_of_atoms = len(self.atoms_list)
        self.energy = 0.0
        keyword = "set basis def2-SVP\n"
        if any(x >= 21 for x in molecule.atomic_number):
            keyword += ' def2-ECP'
        if custom_keyword is not None:
            keyword += custom_keyword
        self.prepare_input(keyword=keyword)

    def prepare_input(self, keyword=""):
        coords = self.start_coords
        f1 = open(self.inp_file, "w")
        if self.scftype is 'uks':
            keyword += 'UKS'
        f1.write(keyword + "\n")
        f1.write("molecule {\n")
        f1.write("{0} {1}\n".format(str(self.charge), str(self.multiplicity)))
        for i in range(self.number_of_atoms):
            f1.write(
                " " + "%3s  %10.7f  %10.7f %10.7f\n" % (self.atoms_list[i], coords[i][0], coords[i][1], coords[i][2]))
        f1.write("}\n")
        f1.write("optimize(\"B97-D\")")
        f1.close()

    def optimize(self, max_cycles=350, gamma=0.0, restart=False, convergence='normal'):
        """
        :return:This object will return the optimization status. It will
        optimize a structure.
        """
        # TODO: Add a return 'CycleExceeded'

        with open(self.out_file, 'w') as fopt:
            out = subp.Popen(["psi4", self.inp_file], stdout=fopt, stderr=fopt)
        out.communicate()
        out.poll()
        exit_status = out.returncode
        if exit_status == 0:
            f = open(self.out_file, "r")
            if "  **** Optimization is complete!" in f.read():
                print("Optimized")
                self.energy = self.get_energy()
                self.optimized_coordinates = self.get_coordinates()
                write_xyz(self.atoms_list,
                          self.optimized_coordinates,
                          self.result_xyz_file, energy=self.energy)
                f.close()
                return True
            else:
                print("Error: OPTIMIZATION PROBABLY FAILED. "
                      "CHECK THE .out FILE FOR PARTIAL OPTIMIZTION ")
                print("Check for partial optimization.")
                print("Location: {}".format(os.getcwd()))
                return False

    def get_energy(self):
        """
        :return:This object will return energy from an psi4 calculation. It will return Hartree units.
        """
        try:
            energy_in_hartrees = None
            with open(self.out_file, "r") as out:
                lines = out.readlines()
                for line in lines:
                    if "Final energy is" in line:
                        energy_in_hartrees = float((line.strip().split())[-1])
            return energy_in_hartrees
        except IOError:
            print("Warning: File ", self.out_file, "was not found.")

    def get_coordinates(self):
        """
        :return:This object will return energy from an psi4 calculation. It will return Hartree units.
        """
        try:
            with open(self.out_file, "r") as out:
                for line in out:
                    if "Final optimized geometry and variables:" in line:
                        break
                coords = []
                for line in out.readlines()[5:5 + self.number_of_atoms + 1]:
                    lc = line.split()
                    if len(lc) == 4:
                        try:
                            _, x, y, z = lc
                            coords.append([float(x), float(y), float(z)])
                        except Exception as e:
                            print(e)
                            print("some problem", line)
            return np.array(coords)
        except IOError:
            print("Warning: File ", self.out_file, "was not found.")


def main():
    from pyar.Molecule import Molecule
    import sys
    mol = Molecule.from_xyz(sys.argv[1])
    method = {'charge': 0, 'multiplicity': 1, 'scftype': 'rhf'}
    geometry = Psi4(mol, method)
    geometry.optimize()


if __name__ == "__main__":
    main()
