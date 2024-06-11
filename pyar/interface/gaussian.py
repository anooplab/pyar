"""
gaussian.py - interface to gaussian program

Copyright (C) 2016 by Surajit Nandi, Anoop Ayyappan, and Mark P. Waller
Indian Institute of Technology Kharagpur, India and Westfaelische Wilhelms
Universitaet Muenster, Germany

This file is part of the pyar project.

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

from pyar import interface
from pyar.interface import SF


class Gaussian(SF):
    def __init__(self, molecule, qc_params):

        super(Gaussian, self).__init__(molecule)

        self.start_coords = molecule.coordinates
        self.inp_file = f'trial_{self.job_name}.com'
        self.out_file = f'trial_{self.job_name}.log'
        self.optimized_coordinates = []
        self.energy = 0.0
        basis = qc_params['basis']
        if basis.lower() == 'def2-svp':
            basis = 'def2SVP'
        self.keyword = f"%nprocshared={qc_params['nprocs']}3\n" \
                       f"%chk=trial_{self.job_name}.chk\n" \
                       f"%mem=2GB\n" \
                       f"# {qc_params['method']} {basis} " \
                       f"SCF=(MaxCycle={qc_params['scf_cycles']})"

    def prepare_input(self):
        coords = self.start_coords
        f1 = open(self.inp_file, "w")
        f1.write(f"{self.keyword}\n\n")
        f1.write(f"{self.job_name}\n\n")
        f1.write(f"{str(self.charge)} {str(self.multiplicity)}\n")
        for i in range(self.number_of_atoms):
            f1.write(f"{self.atoms_list[i]:>3}  {coords[i][0]:10.7f}  {coords[i][1]:10.7f} {coords[i][2]:10.7f}\n")
        f1.write(f"\n")
        f1.close()

    def optimize(self):
        """
        :return:This object will return the optimization status. It will
        optimize a structure.
        """
        # TODO: Add a return 'CycleExceeded'
        logfile = "trial_{}.out".format(self.job_name)

        # max_opt_cycles = options['opt_cycles']
        # gamma = options['gamma']
        # opt_keywords = f"maxcycles={max_opt_cycles}"
        # convergence = options['opt_threshold']
        # if convergence != 'normal':
        #     opt_keywords += f", {convergence}"
        # self.keyword += f" opt=({opt_keywords})"
        self.prepare_input()
        with open(self.out_file, 'w') as fopt:
            out = subp.Popen(["g16", self.inp_file], stdout=fopt, stderr=fopt)
        out.communicate()
        out.poll()
        exit_status = out.returncode
        if exit_status == 0:
            file_pointer = open(self.out_file, "r")
            this_line = file_pointer.readlines()
            check_1 = 0
            check_2 = 0

            for j in this_line:
                if "Optimization completed" in j:
                    check_1 = 1
                if "SCF Done" in j:
                    check_2 = 1

            if ("Normal termination" in this_line[-1]) and check_1 == 1 and check_2 == 1:
                self.energy = self.get_energy()
                self.optimized_coordinates = self.get_coords()
                interface.write_xyz(self.atoms_list, self.optimized_coordinates, self.result_xyz_file, self.job_name,
                                    energy=self.energy)
                file_pointer.close()
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
            if opt_status and 'Standard orientation' in lines:
                pos = i
                coords_lines = t[pos + 5:pos + 5 + self.number_of_atoms]
                for ilines in coords_lines:
                    coordinates.append(ilines.split()[3:6])
                return np.array(coordinates, dtype=float)
        if not opt_status:
            return None

    def get_energy(self):
        """
        :return:This object will return energy from an orca calculation. It will return Hartree units.
        """
        try:
            with open(self.out_file, "r") as out:
                lines_in_file = out.readlines()
                en_steps = []
                for item in lines_in_file:
                    if "SCF Done" in item:
                        en_steps.append(item)
                en_Eh = float((en_steps[-1].strip().split())[4])
            return en_Eh
        except IOError:
            print("Warning: File ", self.out_file, "was not found.")


def main():
    pass


if __name__ == "__main__":
    main()
