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
import logging
import os
import subprocess as subp
import sys

import numpy as np

from pyar.interface import SF, which, write_xyz

xtb_logger = logging.getLogger('pyar.xtb')


class Xtb(SF):

    def __init__(self, molecule, method):
        if which('xtb') is None:
            xtb_logger.error('set XTB path')
            sys.exit()

        super(Xtb, self).__init__(molecule)

        self.cmd = f"xtb {self.start_xyz_file} -opt {method['opt_threshold']}"

        if self.charge != 0:
            self.cmd = "{} -chrg {}".format(self.cmd, self.charge)
        if self.multiplicity != 1:
            self.cmd = "{} -uhf {}".format(self.cmd, self.multiplicity)
        if self.multiplicity == 1 and self.scftype is not 'rhf':  # noqa: F632
            self.cmd = "{} -{}".format(self.cmd, self.scftype)

        self.trajectory_xyz_file = 'traj_' + self.job_name + '.xyz'

    def optimize(self, max_cycles=350, gamma=None, restart=False, convergence='normal'):
        """
        :returns: True,
                  'SCFFailed',
                  'GradFailed',
                  'UpdateFailed',
                  'CycleExceeded',
                  False
        """
        if gamma is not None:
            xtb_logger.error('not implemented in this module. Use xtb_turbo')

        with open('xtb.out', 'w') as output_file_pointer:
            try:
                out = subp.check_call(self.cmd.split(), stdout=output_file_pointer, stderr=output_file_pointer)  # noqa: F841
            except Exception as e:
                xtb_logger.info('    Optimization failed')
                xtb_logger.error(f"      {e}")
                return False

        if os.path.isfile('.xtboptok'):

            write_xyz(self.atoms_list, self.optimized_coordinates, self.result_xyz_file,
                      job_name=self.job_name,
                      energy=self.energy)
            os.rename('xtbopt.log', self.trajectory_xyz_file)
            os.remove('.xtboptok')
            return True
        elif os.path.isfile('.sccnotconverged') or os.path.isfile('NOT_CONVERGED'):
            xtb_logger.info('      SCF Convergence failure in {} run in {}'.format(self.start_xyz_file, os.getcwd()))
            return 'SCFFailed'
        else:
            xtb_logger.info('      Something went wrong with {} run in {}'.format(self.start_xyz_file, os.getcwd()))
            return False

    @property
    def optimized_coordinates(self):
        """"""
        return np.loadtxt('xtbopt.xyz', dtype=float, skiprows=2, usecols=(1, 2, 3))

    @property
    def energy(self):
        if os.path.exists('energy'):
            with open('energy') as fp:
                return float(fp.readlines()[-2].split()[1])
        elif os.path.exists('xtb.out'):
            with open('xtb.out') as fp:
                for line in fp.readlines():
                    if 'total E' in line:
                        energy = float(line.split()[-1])
                    if 'TOTAL ENERGY' in line:
                        energy = float(line.split()[3])
                return energy
        else:
            return None


def main():
    pass


if __name__ == "__main__":
    main()
