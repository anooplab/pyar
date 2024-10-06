import os
import subprocess as subp
import numpy as np
from pyar.interface import SF, write_xyz, which

class OrcaAIQM1(SF):
    def __init__(self, molecule, qc_params):
        super(OrcaAIQM1, self).__init__(molecule)

        self.start_coords = molecule.coordinates
        self.inp_file = f'trial_{self.job_name}.inp'
        self.out_file = f'trial_{self.job_name}.out'
        self.optimized_coordinates = []
        self.energy = 0.0
        self.gamma = qc_params.get('gamma', 0.0)
        self.index = qc_params.get('index', len(molecule.atoms_list))

        # Prepare the ORCA input
        self.prepare_keyword(qc_params)

    def prepare_keyword(self, qc_params):
        keyword = "! Extopt opt xyzfile"
        keyword += "\n%base \"external\""
        keyword += "\n%method"
        keyword += "\n  ProgExt \"/scratch/20cy91r19/app/orca_6_0_0/combine.py\""
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

    def optimize(self):
        self.prepare_input()

        with open(self.out_file, 'w') as fopt:
            out = subp.Popen([which("orca"), self.inp_file], stdout=fopt, stderr=fopt)
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
                    return False
        else:
            print(f"Error: ORCA calculation failed with exit status {exit_status}")
            return False

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