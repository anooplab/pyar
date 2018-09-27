"""
orca.py - interface to ORCA program

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

from interface import SF, write_xyz, which


class Orca(SF):
    def __init__(self, molecule, method, custom_keyword=None):

        super(Orca, self).__init__(molecule)

        self.charge = method['charge']
        self.multiplicity = method['multiplicity']
        self.scftype = method['scftype']

        if (sum(molecule.atomic_number) - self.charge) % 2 == 1 and self.multiplicity == 1:
            self.multiplicity = 2
        else:
            self.multiplicity = method['multiplicity']
        if self.multiplicity % 2 == 0 and self.scftype is 'rhf':
            self.scftype = 'uhf'
        else:
            self.scftype = method['scftype']

        self.start_coords = molecule.coordinates
        self.inp_file = 'trial_' + self.job_name + '.inp'
        self.out_file = 'trial_' + self.job_name + '.out'
        self.optimized_coordinates = []
        self.number_of_atoms = len(self.atoms_list)
        self.energy = 0.0
        keyword = "!BP RI Opt def2-SVP D3BJ KDIIS def2/J PAL3"
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
        f1.write("*xyz {0} {1}\n".format(str(self.charge), str(self.multiplicity)))
        for i in range(self.number_of_atoms):
            f1.write(
                " " + "%3s  %10.7f  %10.7f %10.7f\n" % (self.atoms_list[i], coords[i][0], coords[i][1], coords[i][2]))
        f1.write("*")
        f1.close()

    def optimize(self, max_cycles=350, gamma=0.0, restart=False):
        """
        :return:This object will return the optimization status. It will
        optimize a structure.
        """

        with open(self.out_file, 'w') as fopt:
            out = subp.Popen([which("orca"), self.inp_file], stdout=fopt, stderr=fopt)
        out.communicate()
        out.poll()
        exit_status = out.returncode
        if exit_status == 0:
            f = open(self.out_file, "r")
            line = f.readlines()
            if "****ORCA TERMINATED NORMALLY****" in line[-2]:
                self.energy = self.get_energy()
                self.optimized_coordinates = np.loadtxt(self.inp_file[:-4] + ".xyz", dtype=float, skiprows=2,
                                                        usecols=(1, 2, 3))
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
        :return:This object will return energy from an orca calculation. It will return Hartree units.
        """
        try:
            with open(self.out_file, "r") as out:
                line = out.readlines()
                en_steps = []
                for i in range(len(line)):
                    if "FINAL SINGLE POINT ENERGY" in line[i]:
                        en_steps.append(line[i])
                energy_in_hartrees = float((en_steps[-1].strip().split())[-1])
            return energy_in_hartrees
        except IOError:
            print("Warning: File ", self.out_file, "was not found.")



def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', type=int, required=True, nargs=2,
                        help='scan between two atoms')
    parser.add_argument('-c', '--charge', type=int, default=0,
                        help='charge')
    parser.add_argument('-m', '--multiplicity', type=int, default=1,
                        help='multiplicity')
    parser.add_argument('--scftype', type=str, default= 'rhf', choices= ['rhf', 'uhf'],
                        help='SCF type (rhf/uhf)')
    parser.add_argument("input_file", metavar='file',
                        type=str,
                        help='input coordinate file')
    parser.add_argument('-o', '--opt', action='store_true',
                        help='optimize')
    args = parser.parse_args()
    from Molecule import Molecule
    mol = Molecule.from_xyz(args.input_file)

    method_args = {
        'charge': args.charge,
        'multiplicity': args.multiplicity,
        'scftype': args.scftype,
        'software': 'orca'
    }
    a, b = args.s
    coordinates = mol.coordinates
    import numpy as np
    start_dist = np.linalg.norm(coordinates[a] - coordinates[b])
    final_distance = mol.covalent_radius[a] + mol.covalent_radius[b]
    step = int(abs(final_distance - start_dist)*10)
    c_k = '\n!ScanTS\n% geom scan B '+str(a)+' '+str(b)+ '= '+ str(start_dist)\
          + ', ' + str(final_distance) + ', ' + str(step) + ' end end\n'
    geometry = Orca(mol, method_args, custom_keyword=c_k)
    if args.opt:
        import optimiser
        print('optimising')
        optimiser.optimise(mol, method_args, 0.0, custom_keyword=c_k)
    else:
        print('created input file')


if __name__ == "__main__":
    main()
