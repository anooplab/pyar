"""
mopac.py - interface to mopac program

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
import sys

import numpy as np

from pyar.interface import SF, which


class OBabel(SF):
    def __init__(self, molecule, forcefield=None):

        if which('obabel') is None:
            print('setup OBabel path')
            sys.exit()

        super(OBabel, self).__init__(molecule)

        if not os.path.isfile(self.start_xyz_file):
            molecule.mol_to_xyz(self.start_xyz_file)
        self.optimized_coordinates = []
        self.energy = 0.0

    def optimize(self, max_cycles=350, gamma=0.0, restart=False, convergence='normal'):
        """
        """
        # TODO: Add a return 'CycleExceeded'

        with open('tmp.log', 'w') as logfile, open('tmp.xyz', 'w') as xyzfile:
            try:
                subp.check_call(["obminimize", "-ff", "uff", '-n',
                                 max_cycles, self.start_xyz_file],
                                stdout=xyzfile, stderr=logfile)
            except subp.CalledProcessError as e:
                print('done')
        with open('tmp.xyz') as xyzfile, open(self.result_xyz_file, 'w') as result_xyz_file:
            for line in xyzfile:
                if 'WARNING' not in line:
                    result_xyz_file.write(line)
        self.energy = self.get_energy()
        self.optimized_coordinates = self.get_coords()
        return True

    def get_coords(self):
        """
        :return: It will return coordinates
        """
        return np.loadtxt(self.result_xyz_file, dtype=float, skiprows=2, usecols=(1, 2, 3))

    def get_energy(self):
        """
        """
        with open(self.job_name + '.ene', 'w') as energy_file:
            out = subp.Popen(["obenergy", "-ff", "uff", self.result_xyz_file], stdout=energy_file, stderr=energy_file)
        output, error = out.communicate()
        poll = out.poll()
        exit_status = out.returncode
        if exit_status is None:
            with open(self.job_name + '.ene', 'r') as energy_file:
                return float(energy_file.readlines()[-1].split()[3])


def xyz_to_mopac_input(xyzfile, mopac_input_file, keyword=None):
    """
    Convert xyz file to mopac input

    :param xyzfile: input xyz file
    :param mopac_input_file: Mopac input file to be written
    :param keyword: this is the keyword for optimizations. This parameter
        should be a strings of characters which are mopac keywords
    :return: It will not return anything. It will prepare the input file for
        the purpose given in the keyword. Note that babel will be used to prepare
        the input(.mop) file.
    """
    keyword_line = '-xkPM7' if keyword is None else '-xk' + keyword
    with open('tmp.log', 'w') as fminp:
        subp.call(["babel", "-ixyz", xyzfile, "-omop", mopac_input_file, keyword_line], stderr=fminp, stdout=fminp)
    print(open('tmp.log').readlines()[0])
    os.remove('tmp.log')


def xyz_to_sdf_file(xyz_input_files, sdf_output_file):
    print(xyz_input_files)
    with open('tmp.log', 'w') as fminp:
        subp.call(["babel", "-ixyz"] + xyz_input_files + ["-osdf", sdf_output_file], stderr=fminp, stdout=fminp)
    os.remove('tmp.log')


def make_inchi_string_from_xyz(xyzfile):
    """This function will make a inchi string from a xyz file with
       babel as the tools
    """
    if os.path.isfile(xyzfile):
        with open('OBabel.log', 'w') as ferr:
            inchi = subp.check_output(["babel", "-ixyz", str(xyzfile), "-oinchi"], stderr=ferr)
        return inchi.decode("utf-8").strip()
    else:
        raise IOError("file %s does not exists" % xyzfile)


def make_smile_string_from_xyz(xyzfile):
    """This function will make smile string from a xyz file.
       if more than one xyz file provide, it will return smiles
       in a list. if one xyz file supplied, it will return the
       string
    """
    if os.path.isfile(xyzfile):
        with open('OBabel.log', 'w') as ferr:
            try:
                pre_smile = subp.check_output(["babel", "-ixyz", str(xyzfile), "-osmi", "-xn"], stderr=ferr)
                smile = pre_smile.decode("utf-8").strip()
            except Exception as e:
                ferr.write(str(e))
                smile = ''
            return smile
    else:
        raise IOError("file %s does not exists" % xyzfile)


def main(input_files):
    from pyar import Molecule
    for f in input_files:
        mol = Molecule.Molecule.from_xyz(f)
        g = OBabel(mol)
        g.optimize()


if __name__ == "__main__":
    main(sys.argv[1:])
