"""xTB followed by AIMNet2 refinement backend."""

import logging
import os
import subprocess as subp
import sys
from importlib import resources
from pathlib import Path

import numpy as np

from pyar.backends import SF, require_executable, write_xyz
from pyar.backends.xtb_utils import build_xtb_command

xtb_aimnet2_logger = logging.getLogger('pyar.xtb_aimnet2')

MODEL_PATH = str(resources.files("pyar").joinpath("AIMNet2/models/aimnet2_wb97m-d3_0.jpt"))
SCRIPT_PATH = str(resources.files("pyar").joinpath("AIMNet2/calculators/aimnet2_ase_opt.py"))

class XtbAimnet2(SF):
    """Run xTB, then refine the xTB minimum with AIMNet2."""

    def __init__(self, molecule, method):
        """Build the two-stage xTB + AIMNet2 command sequence."""
        self.xtb_executable = require_executable('xtb', 'xTB')

        super(XtbAimnet2, self).__init__(molecule)

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

        self._validate_runtime_files()
        self.aimnet2_cmd = [
            sys.executable,
            SCRIPT_PATH,
            MODEL_PATH,
            '--traj',
            'result.traj',
            self.xtb_optimized_xyz_file,
            self.aimnet2_optimized_xyz_file,
        ]

        if self.charge != 0:
            self.aimnet2_cmd.extend(["-c", str(self.charge)])

        self.trajectory_xyz_file = 'traj_' + self.job_name + '.xyz'

    @staticmethod
    def _validate_runtime_files():
        """Validate AIMNet2 runtime assets before running jobs."""
        missing = []
        if not Path(MODEL_PATH).is_file():
            missing.append(MODEL_PATH)
        if not Path(SCRIPT_PATH).is_file():
            missing.append(SCRIPT_PATH)
        if missing:
            raise FileNotFoundError(
                "AIMNet2 runtime assets are not bundled in the pyar-chem wheel. "
                "Provide the model files separately or install AIMNet2 from a source checkout. "
                "Missing files: "
                + ", ".join(missing)
            )

    def optimize(self, max_cycles=350, gamma=None, restart=False, convergence='normal'):
        # XTB optimization
        with open('xtb.out', 'w') as output_file_pointer:
            try:
                out = subp.check_call(self.xtb_cmd, stdout=output_file_pointer, stderr=output_file_pointer)
            except Exception as e:
                xtb_aimnet2_logger.info('    XTB optimization failed')
                xtb_aimnet2_logger.error(f"      {e}")
                return False

        if os.path.isfile('.xtboptok'):
            os.rename('xtbopt.xyz', self.xtb_optimized_xyz_file)
            os.remove('.xtboptok')
        else:
            xtb_aimnet2_logger.info('      XTB optimization failed')
            return False

        # AIMNet2 optimization
        with open('aimnet2.out', 'w') as output_file_pointer:
            try:
                subp.check_call(self.aimnet2_cmd, stdout=output_file_pointer, stderr=output_file_pointer)
            except Exception as e:
                xtb_aimnet2_logger.info('    AIMNet2 optimization failed')
                xtb_aimnet2_logger.error(f"      {e}")
                return False

        write_xyz(self.atoms_list, self.optimized_coordinates, self.result_xyz_file,
                  job_name=self.job_name,
                  energy=self.energy)

        return True

    @property
    def xtb_optimized_xyz_file(self):
        return 'xtb_optimized_' + self.job_name + '.xyz'

    @property
    def aimnet2_optimized_xyz_file(self):
        return 'aimnet2_optimized_' + self.job_name + '.xyz'

    @property
    def optimized_coordinates(self):
        return np.loadtxt(self.aimnet2_optimized_xyz_file, dtype=float, skiprows=2, usecols=(1, 2, 3))

    @property
    def energy(self):
        if os.path.exists('aimnet2.out'):
            with open('aimnet2.out') as fp:
                lines = fp.readlines()
                last_line = lines[-1].strip()
                energy_tokens = last_line.split()
                if len(energy_tokens) >= 4:
                    energy = energy_tokens[3]
                    try:
                        energy = float(energy) * 0.0367493  # Convert energy from eV to Hartree
                        return energy
                    except ValueError:
                        print("Error: Unable to convert energy value to float")
                        return None
                else:
                    print("Error: Unexpected format of the last line")
                    return None
        else:
            print("Error: File 'aimnet2.out' not found")
            return None


def main():
    pass


if __name__ == "__main__":
    main()
