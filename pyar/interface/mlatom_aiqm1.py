import logging  # noqa: F401
import pyar.mlatom as ml
from time import sleep  # noqa: F401
from pyar.mlatom.data import molecule  # noqa: F401 
import numpy as np
from pyar import interface
from pyar.interface import SF  # noqa: F811
import os 
# import time
import sys  # noqa: F401

mlatom_logger = logging.getLogger('pyar.mlatom')


class MlatomAiqm1(SF):
    def __init__(self, molecule, qc_params=None):  # noqa: F811
        super(MlatomAiqm1, self).__init__(molecule)
        self.molecule = molecule
        self.qc_params = qc_params
        self.energy = 0.0
        # self.number_of_atoms = molecule.number_of_atoms
        # self.job_name = molecule.name
        self.start_coords = molecule.coordinates
        self.optimized_coordinates = []

        self.inp_file = 'trial_' + self.job_name + '.xyz'
        # self.out_file = 'trial_' + self.job_name + '.log'

        self.aiqm1 = ml.models.methods(method='AIQM1')

    def prepare_input(self):
        coords = self.start_coords      
        with open(self.inp_file, "w") as f1:
            f1.write(str(self.number_of_atoms) + "\n")
            f1.write("trial_" + self.job_name + "\n")
            for i in range(self.number_of_atoms):
                f1.write(
                    " " + "%3s  %10.7f  %10.7f %10.7f\n" % (self.atoms_list[i], coords[i][0], coords[i][1], coords[i][2]))

    def read_molecule_from_file(self, filename):
        mol = ml.data.molecule()
        return mol.read_from_xyz_file(filename)

    
    
    def optimize(self):
        self.prepare_input()
        try:
            initial_molecule = self.read_molecule_from_file(self.inp_file)
            # initial_molecule = ml.data.molecule.read_from_xyz_file(self, self.inp_file)
            initial_molecule.charge = self.charge
            initial_molecule.multiplicity = self.multiplicity
            # future = ml.optimize_geometry(model=self.aiqm1, initial_molecule=initial_molecule)
            # self.optimized_molecule = future.result().optimized_molecule
            self.optimized_molecule =  ml.optimize_geometry(model=self.aiqm1, initial_molecule=initial_molecule).optimized_molecule
            sleep(90)
    
            self.out_file = 'gaussian.log'
            
            # Check if the file exists
            if not os.path.isfile(self.out_file):
                mlatom_logger.info("Error: File does not exist.")
                return None
    
            file_pointer = open(self.out_file, "r")
            this_line = file_pointer.readlines()
            check_1 = 0 
            check_2 = 0
    
            for j in this_line:
                if "Optimization completed" in j:
                    check_1 = 1
                if "SCF Done" in j: 
                   check_2 = 1
    
            if ("Normal termination" in this_line[-1]) and check_1 == 1 and check_2 == 1:
                self.energy = self.get_energy()
                self.optimized_coordinates = self.get_coords()
                interface.write_xyz(self.atoms_list, self.optimized_coordinates, self.result_xyz_file, self.job_name,
                                    energy=self.energy)
                file_pointer.close()
                return True
            else:
                mlatom_logger.info("Error: OPTIMIZATION PROBABLY FAILED.")
                mlatom_logger.info("Last line in log file: {}".format(this_line[-1]))
                mlatom_logger.info("Location: {}".format(os.getcwd()))
                return None
                
        except Exception as e:
            mlatom_logger.info('    Optimization failed:', f"      {e}")
            return None
               
        
        
    def get_coords(self):
        """
        :return: coords It will return coordinates
        """
        opt_status = False
        coordinates = []
        with open(self.out_file) as v:
            t = v.readlines()
        for i, lines in enumerate(t):
            if 'Stationary point found.' in lines:
                opt_status = True
            if opt_status and 'Standard orientation' in lines:
                pos = i
                coords_lines = t[pos + 5:pos + 5 + self.number_of_atoms]
                for ilines in coords_lines:
                    coordinates.append(ilines.split()[3:6])
                return np.array(coordinates, dtype=float)
        if not opt_status:
            return None
        
    def get_energy(self):
        """
        :return: This object will return energy from an gaussian calculation. It will return energy in Hartree units.
        """
        try:
            with open(self.out_file, "r") as out:
                lines_in_file = out.readlines()
                last_energy_line = None
                for line in reversed(lines_in_file):
                    if "Recovered energy" in line:
                        last_energy_line = line.strip()
                        break
                if last_energy_line:
                    en_Eh = float(last_energy_line.split("=")[1].split()[0])
                    return en_Eh
                else:
                    print("Warning: No 'Recovered energy' line found in file.")
                    return None
        except IOError:
            print("Warning: File", self.out_file, "was not found.")

    

def main():
    pass

if __name__ == '__main__':
    main()


