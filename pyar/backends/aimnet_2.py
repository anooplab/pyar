"""AIMNet2 backend implementation."""

import logging
import os
import subprocess as subp
import sys
from importlib import resources
from pathlib import Path

import numpy as np
import pyar  # noqa: F401
from pyar.backends import SF, write_xyz  # noqa: F401
from pyar.optional_dependencies import optional_dependency_error

Aimnet2_logger = logging.getLogger('pyar.aimnet-2')

try:
    import torch  # noqa: E402
    from pyar.AIMNet2.calculators import aimnet2_ase_opt  # noqa: F401
    from pyar.AIMNet2.calculators import aimnet2ase  # noqa: F401
except ImportError as exc:
    module_name = getattr(exc, "name", None) or "torch"
    raise optional_dependency_error(module_name, feature="AIMNet2 backend", extra="aimnet2") from exc

torch.backends.cuda.matmul.allow_tf32 = False
torch.backends.cudnn.allow_tf32 = False



device = torch.device('cpu')
Aimnet2_logger.debug("AIMNet2 torch device: %s", device)

MODEL_PATH = str(resources.files("pyar").joinpath("AIMNet2/models/aimnet2_wb97m-d3_0.jpt"))
SCRIPT_PATH = str(resources.files("pyar").joinpath("AIMNet2/calculators/aimnet2_ase_opt.py"))

class Aimnet2(SF):
    def __init__(self, molecule, qc_params):
        super(Aimnet2, self).__init__(molecule)

        self.trajectory_xyz_file = 'traj_' + self.job_name + '.traj'
        self.molecule = molecule
        self.qc_params = qc_params
        # self.energy = 0.0
        self.number_of_atoms = molecule.number_of_atoms
        self.job_name = molecule.name
        self.start_coords = molecule.coordinates
        # self.optimized_coordinates = []
        self.inp_file = 'trial_' + self.job_name + '.xyz'
        self.inp_min_file = 'trial_' + self.job_name + '_min.xyz'
        self.out_file = 'trial_' + self.job_name + '.out'

        self._validate_runtime_files()
        self.cmd = [
            sys.executable,
            SCRIPT_PATH,
            MODEL_PATH,
            '--traj',
            'result.traj',
            self.inp_file,
            self.inp_min_file,
        ]
        if self.charge != 0:
            self.cmd.extend(["-c", str(self.charge)])

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

    def prepare_input(self):
        coords = self.start_coords
        with open(self.inp_file, "w") as f1:
            f1.write(str(self.number_of_atoms) + "\n")
            f1.write("trial_" + self.job_name + "\n")
            for i in range(self.number_of_atoms):
                f1.write(
                    " " + "%3s  %10.7f  %10.7f %10.7f\n" % (self.atoms_list[i], coords[i][0], coords[i][1], coords[i][2]))

    def optimize(self):
        """
        :returns: True,
                  'SCFFailed',
                  'GradFailed',
                  'UpdateFailed',
                  'CycleExceeded',
                  False
        """


        with open('aimnet2.out', 'w') as output_file_pointer:
            try:
                subp.check_call(self.cmd, stdout=output_file_pointer, stderr=output_file_pointer)
            except Exception as e:
                Aimnet2_logger.info('    Optimization failed')
                Aimnet2_logger.error(f"      {e}")
                return False

            write_xyz(self.atoms_list, self.optimized_coordinates, self.result_xyz_file,
                      job_name=self.job_name,
                      energy=self.energy)

            return True



    @property
    def optimized_coordinates(self):
        """"""
        return np.loadtxt(f'{self.inp_min_file}', dtype=float, skiprows=2, usecols=(1, 2, 3))

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
