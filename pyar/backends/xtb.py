"""Pure xTB geometry optimization backend."""
import logging
import os
import subprocess as subp
import sys

import numpy as np

from pyar.backends import SF, require_executable, write_xyz
from pyar.backends.xtb_utils import build_xtb_command

xtb_logger = logging.getLogger('pyar.xtb')


def xtb_supports_gxtb(xtb_executable):
    try:
        help_text = subp.check_output(
            [xtb_executable, '--help'],
            stderr=subp.STDOUT,
            text=True,
        )
        return '--gxtb' in help_text
    except Exception as e:
        xtb_logger.warning(f'Could not check xTB --gxtb support: {e}')
        return False


class Xtb(SF):
    """Run a standalone xTB optimization for a single molecule."""

    def __init__(self, molecule, method):
        """Build the xTB command for a molecule and QC parameter set."""
        self.xtb_executable = require_executable('xtb', 'xTB')

        super(Xtb, self).__init__(molecule)

        self.cmd = build_xtb_command(
            self.xtb_executable,
            self.start_xyz_file,
            method,
            opt_threshold=method['opt_threshold'],
        )

        if xtb_supports_gxtb(self.xtb_executable):
            self.cmd.append('--gxtb')
            xtb_logger.info('Using xTB with --gxtb flag')
        else:
            xtb_logger.info('Installed xTB does not support --gxtb; running without it')

        self.trajectory_xyz_file = 'traj_' + self.job_name + '.xyz'

    def optimize(self, max_cycles=350, gamma=None, restart=False, convergence='normal'):
        """Execute xTB and return a success flag or failure status string."""
        if gamma is not None:
            xtb_logger.error('not implemented in this module. Use xtb_turbo')

        with open('xtb.out', 'w') as output_file_pointer:
            try:
                out = subp.check_call(
                    self.cmd,
                    stdout=output_file_pointer,
                    stderr=output_file_pointer,
                )
            except Exception as e:
                xtb_logger.info('    Optimization failed')
                xtb_logger.error(f'      {e}')
                return False

        if os.path.isfile('.xtboptok'):
            write_xyz(
                self.atoms_list,
                self.optimized_coordinates,
                self.result_xyz_file,
                job_name=self.job_name,
                energy=self.energy,
            )

            if os.path.isfile('xtbopt.log'):
                os.rename('xtbopt.log', self.trajectory_xyz_file)

            os.remove('.xtboptok')
            return True

        elif os.path.isfile('.sccnotconverged') or os.path.isfile('NOT_CONVERGED'):
            xtb_logger.info(
                '      SCF Convergence failure in {} run in {}'.format(
                    self.start_xyz_file,
                    os.getcwd(),
                )
            )
            return 'SCFFailed'

        else:
            xtb_logger.info(
                '      Something went wrong with {} run in {}'.format(
                    self.start_xyz_file,
                    os.getcwd(),
                )
            )
            return False

    @property
    def optimized_coordinates(self):
        """Read the optimized xTB coordinates from ``xtbopt.xyz``."""
        return np.loadtxt(
            'xtbopt.xyz',
            dtype=float,
            skiprows=2,
            usecols=(1, 2, 3),
        )

    @property
    def energy(self):
        """Read the final total energy from xTB output files."""
        if os.path.exists('energy'):
            with open('energy') as fp:
                return float(fp.readlines()[-2].split()[1])

        elif os.path.exists('xtb.out'):
            energy = None
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


if __name__ == '__main__':
    main()
