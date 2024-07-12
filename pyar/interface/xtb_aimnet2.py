import logging
import os
import subprocess as subp
import sys

import numpy as np
import torch

from pyar.interface import SF, which, write_xyz
import pkg_resources
import os

xtb_aimnet2_logger = logging.getLogger('pyar.xtb_aimnet2')

device = torch.device('cpu')
print(device)

model_path = pkg_resources.resource_filename('pyar', 'AIMNet2/models/aimnet2_wb97m-d3_0.jpt')
aimnet2_script = pkg_resources.resource_filename('pyar', 'AIMNet2/calculators/aimnet2_ase_opt.py')
# Load the model
aimnet2 = torch.jit.load(model_path, map_location=device)

class XtbAimnet2(SF):
    def __init__(self, molecule, method):
        if which('xtb') is None:
            xtb_aimnet2_logger.error('set XTB path')
            sys.exit()

        if which('python') is None:
            xtb_aimnet2_logger.error('set python3 path')
            sys.exit()

        super(XtbAimnet2, self).__init__(molecule)

        self.xtb_cmd = f"xtb {self.start_xyz_file} -opt {method['opt_threshold']}"

        if self.charge != 0:
            self.xtb_cmd = "{} -chrg {}".format(self.xtb_cmd, self.charge)
        if self.multiplicity != 1:
            self.xtb_cmd = "{} -uhf {}".format(self.xtb_cmd, self.multiplicity)
        if self.multiplicity == 1 and self.scftype is not 'rhf':
            self.xtb_cmd = "{} -{}".format(self.xtb_cmd, self.scftype)

        self.aimnet2_cmd = f"python {aimnet2_script} {model_path} --traj result.traj {self.xtb_optimized_xyz_file} {self.aimnet2_optimized_xyz_file}"

        if self.charge != 0:
            self.aimnet2_cmd = "{} -c {}".format(self.aimnet2_cmd, self.charge)

        self.trajectory_xyz_file = 'traj_' + self.job_name + '.xyz'

    def optimize(self, max_cycles=350, gamma=None, restart=False, convergence='normal'):
        # XTB optimization
        with open('xtb.out', 'w') as output_file_pointer:
            try:
                out = subp.check_call(self.xtb_cmd.split(), stdout=output_file_pointer, stderr=output_file_pointer)
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
                out = subp.check_call(self.aimnet2_cmd.split(), stdout=output_file_pointer, stderr=output_file_pointer)
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