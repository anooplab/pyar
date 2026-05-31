"""xTB followed by AIQM1 refinement backend."""

import logging
import os
import subprocess as subp
import sys
from importlib import resources

import numpy as np

from pyar.backends import SF, require_executable, write_xyz
from pyar.backends.xtb_utils import build_xtb_command

xtb_aiqm1_logger = logging.getLogger('pyar.xtb_aiqm1')

aiqm1_opt = str(resources.files('pyar').joinpath('backends/mlopt.py'))

class XtbAIQM1(SF):
    """Run xTB, then refine the xTB minimum with AIQM1."""

    def __init__(self, molecule, method):
        """Build the two-stage xTB + AIQM1 command sequence."""
        self.xtb_executable = require_executable('xtb', 'xTB')

        super(XtbAIQM1, self).__init__(molecule)

        self.xtb_cmd = build_xtb_command(
            self.xtb_executable,
            self.start_xyz_file,
            {
                **method,
                "charge": self.charge,
                "multiplicity": self.multiplicity,
                "scftype": self.scftype,
            },
            opt_threshold=method['opt_threshold'],
        )

        self.aiqm1_cmd = f"{sys.executable} {aiqm1_opt} {self.xtb_optimized_xyz_file} -c {self.charge} -m {self.multiplicity}  {self.aiqm1_optimized_xyz_file}"

        self.trajectory_xyz_file = 'traj_' + self.job_name + '.xyz'

    def optimize(self, max_cycles=350, gamma=None, restart=False, convergence='normal'):
        # XTB optimization
        with open('xtb.out', 'w') as output_file_pointer:
            try:
                out = subp.check_call(self.xtb_cmd, stdout=output_file_pointer, stderr=output_file_pointer)
            except Exception as e:
                xtb_aiqm1_logger.info('    XTB optimization failed')
                xtb_aiqm1_logger.error(f"      {e}")
                return False

        if os.path.isfile('.xtboptok'):
            os.rename('xtbopt.xyz', self.xtb_optimized_xyz_file)
            os.remove('.xtboptok')
        else:
            xtb_aiqm1_logger.info('      XTB optimization failed')
            return False

        # AIQM1 optimization
        with open('aiqm1.out', 'w') as output_file_pointer:
            try:
                out = subp.check_call(self.aiqm1_cmd.split(), stdout=output_file_pointer, stderr=output_file_pointer)
            except Exception as e:
                xtb_aiqm1_logger.info('    AIQM1 optimization failed')
                xtb_aiqm1_logger.error(f"      {e}")
                return False

        write_xyz(self.atoms_list, self.optimized_coordinates, self.result_xyz_file,
                  job_name=self.job_name,
                  energy=self.energy)

        return True

    @property
    def xtb_optimized_xyz_file(self):
        return 'xtb_optimized_' + self.job_name + '.xyz'

    @property
    def aiqm1_optimized_xyz_file(self):
        return 'aiqm1_optimized_' + self.job_name + '.xyz'

    @property
    def optimized_coordinates(self):
        return np.loadtxt(self.aiqm1_optimized_xyz_file, dtype=float, skiprows=2, usecols=(1, 2, 3))

    @property
    def energy(self):
        if os.path.exists('aiqm1.out'):
            with open('aiqm1.out') as fp:
                lines = fp.readlines()
                last_line = lines[-1].strip()
                energy_tokens = last_line.split()
                if len(energy_tokens) >= 3:
                    energy = float(energy_tokens[2])
                    return energy
                else:
                    print("Error: Unexpected format of the last line")
                    return None
        else:
            print("Error: File 'aiqm1.out' not found")
            return None


def main():
    pass


if __name__ == "__main__":
    main()
