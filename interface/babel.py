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


class obabel(object):
    def __init__(self, molecule, forcefield=None):
        self.job_name = molecule.name
        self.start_xyz_file = 'trial_'+self.job_name+'.xyz'
        if not os.path.isfile(self.start_xyz_file):
            molecule.mol_to_xyz(self.start_xyz_file)
        self.result_xyz_file = 'result_'+self.job_name+'.xyz'
        self.optimized_coordinates = []
        self.energy = 0.0

    def optimize(self, gamma=0.0):
        """
        """
        with open('tmp.log', 'w') as logfile, open('tmp.xyz', 'w') as xyzfile:
            out = subp.Popen(["obminimize", "-ff", "uff", self.start_xyz_file], stdout=xyzfile, stderr=logfile)
        output, error = out.communicate()
        poll = out.poll()
        exit_status = out.returncode
        if exit_status == 1:
            with open('tmp.xyz') as xyzfile, open(self.result_xyz_file, 'w') as result_xyz_file:
                for line in xyzfile:
                    if not 'WARNING' in line:
                        result_xyz_file.write(line)
            self.energy = self.get_energy()
            self.optimized_coordinates = self.get_coords()
            return True
        os.remove('tmp.log')
        return exit_status

    def get_coords(self):
        """
        :param out_file: This is the output file in which the final xyz coordinates will be
        written
        :return: It will return coordinates
        """
        return np.loadtxt(self.result_xyz_file, dtype=float, skiprows=2, usecols=(1, 2, 3))


    def get_energy(self):
        """
        """
        with open(self.job_name+'.ene', 'w') as energy_file:
            out = subp.Popen(["obenergy", "-ff", "uff", self.result_xyz_file], stdout=energy_file, stderr=energy_file)
        output, error = out.communicate()
        poll = out.poll()
        exit_status = out.returncode
        if exit_status == 1:
            with open(self.job_name+'.ene', 'r') as energy_file:
                energy = float(energy_file.readlines()[-1].split()[3])
                return energy


def xyz_to_mopac_input(xyzfile, mopac_input_file, keyword=None):
    """

    :param xyzfile: input xyz file
    :param mopac_input_file: Mopac input file to be written
    :param keyword: this is the keyword for optimizations. This parameter
    should be a strings of characters which are mopac keywords
    :return: It will not return anything. It will prepare the input file for
    the purpose given in the keyword. Note that babel will be used to prepare
    the input(.mop) file.
    """
    if keyword is None:
        keyword_line = "-xkPM7"
    else:
        keyword_line = "-xk" + keyword
    with open('tmp.log', 'w') as fminp:
        subp.call(["babel", "-ixyz", xyzfile, "-omop", mopac_input_file, keyword_line], stderr=fminp, stdout=fminp)
    print(open('tmp.log').readlines()[0])
    os.remove('tmp.log')


def xyz_to_sdf_file(xyz_input_files, sdf_output_file):
    print(xyz_input_files)
    with open('tmp.log', 'w') as fminp:
        subp.call(["babel", "-ixyz"]+ xyz_input_files+ ["-osdf", sdf_output_file] , stderr=fminp, stdout=fminp)
    os.remove('tmp.log')

def make_inchi_string_from_xyz(xyzfile):
    """This function will make a inchi string from a xyz file with
       babel as the tools
    """
    if os.path.isfile(xyzfile):
        inchi = subp.check_output(["babel", "-ixyz", str(xyzfile), "-oinchi"])
        return inchi
    else:
        raise IOError("file %s does not exists" % xyzfile)


def make_smile_string_from_xyz(xyzfile):
    """This function will make smile string from a xyz file.
       if more than one xyz file provide, it will return smiles
       in a list. if one xyz file supplied, it will return the
       string
    """
    if os.path.isfile(xyzfile):
        pre_smile = subp.check_output(["babel", "-ixyz", str(xyzfile), "-osmi"])
        smile = pre_smile.split()[0]
        # print "smile string inside make_smile_string_from_xyz() is: ", smile
        return smile
    else:
        raise IOError("file %s does not exists" % xyzfile)


def min_all(input_sdf_file):
    pass

def write_xyz(job_name, atoms_list, coordinates, filename, energy=0.0):
    with open(filename, 'w') as fp:
        fp.write("%3d\n" % len(coordinates))
        fp.write(job_name + ':' + str(energy) + '\n')
        for a, c in zip(atoms_list, coordinates):
            fp.write("{:<2}{:12.5f}{:12.5f}{:12.5f}\n".format(a, c[0], c[1], c[2]))


def main(input_files):
    import Molecule
    for f in input_files:
        mol = Molecule.Molecule.from_xyz(f)
        g = obabel(mol)
        g.optimize()


if __name__ == "__main__":
    main(sys.argv[1:])