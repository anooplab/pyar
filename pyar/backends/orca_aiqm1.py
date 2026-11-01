"""ORCA followed by AIQM1 backend implementation."""
from pyar.backend_errors import backend_operation
from pyar.backend_errors import BackendExecutionError

import os
import subprocess as subp

import numpy as np

from pyar.backends import SF, require_executable, write_xyz

EXTERNAL_PROGRAM_ENV = "PYAR_ORCA_AIQM1_PROGEXT"


def _external_program(qc_params):
    """Return the ORCA ``ProgExt`` wrapper supplied by the user, failing fast if absent."""
    program = qc_params.get('orca_aiqm1_progext') or os.environ.get(EXTERNAL_PROGRAM_ENV)
    if not program:
        raise ValueError(
            "The orca-aiqm1 backend needs the ORCA external-program wrapper that combines "
            f"AIQM1 with the AFIR force. Set {EXTERNAL_PROGRAM_ENV} (or the "
            "'orca_aiqm1_progext' backend setting) to its path."
        )
    if not os.path.isfile(program):
        raise FileNotFoundError(f"orca-aiqm1 external program not found: {program}")
    return os.path.abspath(program)


class OrcaAIQM1(SF):
    def __init__(self, molecule, qc_params):
        super(OrcaAIQM1, self).__init__(molecule)
        self.orca_executable = require_executable("orca", "ORCA")

        self.start_coords = molecule.coordinates
        self.inp_file = f'trial_{self.job_name}.inp'
        self.out_file = f'trial_{self.job_name}.out'
        self.optimized_coordinates = []
        self.energy = 0.0
        self.gamma = qc_params.get('gamma', 0.0)
        self.index = qc_params.get('index', len(molecule.atoms_list))
        self.external_program = _external_program(qc_params)

        # Prepare the ORCA input
        self.prepare_keyword(qc_params)

    def prepare_keyword(self, qc_params):
        keyword = "! Extopt opt xyzfile"
        keyword += "\n%base \"external\""
        keyword += "\n%method"
        keyword += f"\n  ProgExt \"{self.external_program}\""
        keyword += f"\n  Ext_Params \"--force {self.gamma} --index {self.index}\""
        keyword += "\nend"

        self.keyword = keyword

    def set_gamma(self, gamma):
        self.gamma = gamma
        self.prepare_keyword({'gamma': gamma, 'index': self.index})

    def prepare_input(self):
        with open(self.inp_file, "w") as f:
            f.write(self.keyword + "\n\n")
            f.write(f"*xyz {self.charge} {self.multiplicity}\n")
            for i, coord in enumerate(self.start_coords):
                f.write(f"{self.atoms_list[i]:3s} {coord[0]:10.7f} {coord[1]:10.7f} {coord[2]:10.7f}\n")
            f.write("*")

    @backend_operation
    def optimize(self):
        self.prepare_input()

        with open(self.out_file, 'w') as fopt:
            out = subp.Popen([self.orca_executable, self.inp_file], stdout=fopt, stderr=fopt)
        out.communicate()
        exit_status = out.returncode

        if exit_status == 0:
            with open(self.out_file, "r") as f:
                lines = f.readlines()
                if "****ORCA TERMINATED NORMALLY****" in lines[-2]:
                    self.energy = self.get_energy()
                    self.optimized_coordinates = np.loadtxt(f"{self.inp_file[:-4]}.xyz", dtype=float, skiprows=2, usecols=(1, 2, 3))
                    write_xyz(self.atoms_list, self.optimized_coordinates, self.result_xyz_file, energy=self.energy)
                    return True
                else:
                    print("Error: OPTIMIZATION PROBABLY FAILED. CHECK THE .out FILE FOR PARTIAL OPTIMIZATION")
                    print(f"Location: {os.getcwd()}")
                    raise BackendExecutionError("orca_aiqm1 optimization failed; inspect the retained backend output")
        else:
            print(f"Error: ORCA calculation failed with exit status {exit_status}")
            raise BackendExecutionError("orca_aiqm1 optimization failed; inspect the retained backend output")

    def get_energy(self):
        try:
            with open(self.out_file, "r") as out:
                lines = out.readlines()
                en_steps = [item for item in lines if "FINAL SINGLE POINT ENERGY" in item]
                if en_steps:
                    energy_in_hartrees = float(en_steps[-1].strip().split()[-1])
                else:
                    energy_in_hartrees = 0.0
            return energy_in_hartrees
        except IOError:
            print(f"Warning: File {self.out_file} was not found.")
            return 0.0

def main():
    pass

if __name__ == "__main__":
    main()
