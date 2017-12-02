"""
mopac.py - interface to mopac program

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
import interface.babel

class Mopac(object):
    def __init__(self, molecule, charge=0, multiplicity=1, scftype='rhf'):
        self.job_name = molecule.name
        self.charge = charge
        self.multiplicity = multiplicity
        self.scftype = scftype
        self.start_xyz_file = 'trial_' + self.job_name + '.xyz'
        self.result_xyz_file = 'result_' + self.job_name + '.xyz'
        self.inp_file = 'trial_' + self.job_name + '.mop'
        self.arc_file = 'trial_' + self.job_name + '.arc'
        self.atoms_list = molecule.atoms_list
        self.start_coords = molecule.coordinates
        self.optimized_coordinates = []
        self.number_of_atoms = len(self.atoms_list)
        self.energy = 0.0
        keyword = "PM7 PRECISE LET DDMIN=0.0 " \
                  "CYCLES=10000 charge={}".format(charge)
        self.prepare_input(keyword=keyword)

    def prepare_input(self, keyword=""):
        """
        :param keyword: this is the keyword for optimizations. This parameter
        should be a strings of characters which are mopac keywords
        :return: It will not return anything. It will prepare the input file for
        the purpose given in the keyword. Note that babel will be used to prepare
        the input(.mop) file.
        """
        if not keyword:
            keyword_line = "-xkPM7"
        elif keyword:
            keyword_line = "-xk" + keyword
        else:
            print("keyword:", keyword, "not recognized.")
            print("program will stop")
            sys.exit()
        with open('tmp.log', 'w') as fminp:
            out = subp.Popen(["babel", "-ixyz", self.start_xyz_file, "-omop", self.inp_file, keyword_line],
                             stdout=fminp, stderr=fminp)
        output, error = out.communicate()
        poll = out.poll()
        exit_status = out.returncode
        os.remove('tmp.log')
        return exit_status

    def optimize(self, gamma=None):
        """
        :return:This object will return the optimization status. It will
        optimize a structure.
        """
        logfile = "trial_{}.log".format(self.job_name)
        with open(logfile, 'w') as fopt:
            out = subp.Popen(["mopac", self.inp_file], stdout=fopt, stderr=fopt)
        out.communicate()
        out.poll()
        exit_status = out.returncode
        if exit_status == 0:
            if os.path.exists(self.arc_file):
                self.energy = self.get_energy()
                self.optimized_coordinates = self.get_coords()
                interface.babel.write_xyz(self.optimized_coordinates, self.result_xyz_file, energy=self.energy)
                return True
            else:
                print("Error: File ", self.arc_file, "was not found.")
                print("Check for partial optimization.")
                print("Location: {}".format(os.getcwd()))
                return False

    def get_energy(self):
        """
        :return:This object will return energy from a mopac calculation. It will return both the kj/mol and
        kcal/mol units.
        """
        en_kcal = 0.0
        en_kj = 0.0
        try:
            with open(self.arc_file, 'r') as arc_out:
                arc_cont = arc_out.readlines()
            for lines in arc_cont:
                if "HEAT OF FORMATION" in lines:
                    line_cont = lines.split('=')
                    en_kcal = float(line_cont[1].split()[0])
                    en_kj = float(line_cont[2].split()[0])
                    break
        except IOError:
            print("Warning: File ", self.arc_file, "was not found.")
        return en_kcal / 627.51

    def get_coords(self):
        """
        :param out_file: This is the output file in which the final xyz coordinates will be
        written
        :return: It will return coordinates
        """
        number_of_atoms = None
        with open(self.arc_file) as arc_out:
            arc_cont = arc_out.readlines()
        for lines in arc_cont:
            if 'Empirical Formula:' in lines:
                number_of_atoms = int(lines.split()[-2])
        coordinates = arc_cont[-(number_of_atoms + 1):-1]
        coords = []
        atoms_list = []
        for i in coordinates:
            c = i.split()
            try:
                coords.append(np.array([c[1], c[3], c[5]], dtype=str).astype(np.float))
            except ValueError:
                return
            atoms_list.append(c[0])
        return np.array(coords)


def main():
    pass


if __name__ == "__main__":
    main()
