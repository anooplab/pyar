import logging
import pyar  # noqa: F401
# from time import sleep
from pyar.interface import SF, write_xyz, which  # noqa: F401
import os
import subprocess as subp
import numpy as np
import sys
import pkg_resources
import os

mlatom_logger = logging.getLogger('pyar.mlatom')

aiqm1_opt = pkg_resources.resource_filename('pyar', 'interface/mlopt.py')

class AIQM1(SF):
    def __init__(self, molecule, qc_params):
        if which('python') is None:
            mlatom_logger.error('set python3 path')
            sys.exit()


        super(AIQM1, self).__init__(molecule)

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
        
        self.cmd = f"python {aiqm1_opt}  {self.inp_file} -c {self.charge} -m {self.multiplicity}  {self.inp_min_file}"
        if self.charge != 0:
            self.cmd = "{} -c {}".format(self.cmd, self.charge)

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
        

        with open('mlatom.out', 'w') as output_file_pointer:
            try:
                out = subp.check_call(self.cmd.split(), stdout=output_file_pointer, stderr=output_file_pointer)  # noqa: F841
            except Exception as e:
                mlatom_logger.info('    Optimization failed')
                mlatom_logger.error(f"      {e}")
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
        if os.path.exists('mlatom.out'):
            with open('mlatom.out') as fp:
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
            print("Error: File 'mlatom.out' not found")
            return None
    
        
    
def main():
    pass


if __name__ == "__main__":
    main()
