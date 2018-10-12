"""
xtb.py - interface to mopac program
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

from interface import SF, which, write_xyz


class Xtb(SF):

    def __init__(self, molecule, method):
        if which('xtb') is None:
            print('set XTB path')
            sys.exit()

        super(Xtb, self).__init__(molecule)

        charge = method['charge']
        scftype = method['scftype']
        multiplicity = method['multiplicity']

        self.cmd = "xtb {} -opt vtight".format(self.start_xyz_file)
        if charge != 0:
            self.cmd = "{} -chrg {}".format(self.cmd, charge)
        if multiplicity != 1:
            self.cmd = "{} -uhf {}".format(self.cmd, multiplicity)
        if multiplicity == 1 and scftype is not 'rhf':
            self.cmd = "{} -{}".format(self.cmd, scftype)

        self.trajectory_xyz_file = 'traj_' + self.job_name + '.xyz'

    def optimize(self, max_cycles=350, gamma=0.0, restart=False, convergence='normal'):
        """
        :returns: True,
                  'SCFFailed',
                  'GradFailed',
                  'UpdateFailed',
                  'CycleExceeded',
                  False
        """
        if gamma > 0.0:
            print('not implemented in this module. Use xtb_turbo')

        with open('xtb.out', 'w') as fout:
            try:
                out = subp.check_call(self.cmd.split(), stdout=fout, stderr=fout)
            except:
                print('Optimization failed')
                return False

        if os.path.isfile('.xtboptok'):

            write_xyz(self.atoms_list, self.optimized_coordinates, self.result_xyz_file,
                      job_name=self.job_name,
                      energy=self.energy)
            os.rename('xtbopt.log', self.trajectory_xyz_file)
            os.remove('.xtboptok')
            return True
        elif os.path.isfile('.sccnotconverged'):
            print('SCF Convergence failure in {} run in {}'.format(self.start_xyz_file, os.getcwd()))
            return 'SCFFailed'
        else:
            print('Something went wrong with {} run in {}'.format(self.start_xyz_file, os.getcwd()))
            return False

    @property
    def optimized_coordinates(self):
        """"""
        return np.loadtxt('xtbopt.xyz', dtype=float, skiprows=2, usecols=(1, 2, 3))

    @property
    def energy(self):
        with open('energy') as fp:
            return float(fp.readlines()[-2].split()[1])


def main():

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--charge', type=int, default=0,
                        help='charge')
    parser.add_argument('-m', '--multiplicity', type=int, default=1,
                        help='multiplicity')
    parser.add_argument('--scftype', type=str, default='rhf',
                        choices=['rhf', 'uhf'],
                        help='SCF type (rhf/uhf)')
    parser.add_argument("input_file", metavar='file',
                        type=str,
                        help='input coordinate file')
    parser.add_argument('--scan', type=int, nargs=2,
                        help='scan between two atoms')
    parser.add_argument('--cycle', type=int, default=350,
                        help='maximum number of optimization cycles')
    args = parser.parse_args()

    from Molecule import Molecule
    mol = Molecule.from_xyz(args.input_file)
    method_args = {
        'charge': args.charge,
        'multiplicity': args.multiplicity,
        'scftype': args.scftype,
        'software': 'xtb'
    }
    Xtb(mol, method_args)

    import optimiser

    print('optimising')
    optimiser.optimise(mol, method_args)

    if args.scan:
        pass


if __name__ == "__main__":
    main()
