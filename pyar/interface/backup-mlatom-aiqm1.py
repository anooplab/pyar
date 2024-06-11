import logging
import pyar  # noqa: F401
import pyar.mlatom as ml
from time import sleep
from pyar.mlatom.data import molecule  # noqa: F401
from pyar.interface import SF, write_xyz, which  # noqa: F401

mlatom_logger = logging.getLogger('pyar.mlatom')
# optimiser_logger = logging.getLogger('pyar.optimiser')


class MlatomAiqm1(SF):
    def __init__(self, molecule, qc_params=None):  # noqa: F811
        super(MlatomAiqm1, self).__init__(molecule)
        self.molecule = molecule
        self.qc_params = qc_params
        # self.multiplicity = molecule.multiplicity
        # self.charge = molecule.charge
        self.energy = 0.0
        self.number_of_atoms = molecule.number_of_atoms
        self.job_name = molecule.name
        self.start_coords = molecule.coordinates
        self.optimized_coordinates = []
        self.inp_file = 'trial_' + self.job_name + '.xyz'
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
           initial_molecule.charge = self.charge
           initial_molecule.multiplicity = self.multiplicity
           self.optimized_molecule = ml.optimize_geometry(model=self.aiqm1, initial_molecule=initial_molecule).optimized_molecule
           self.energy = self.optimized_molecule.energy
        #    coordinates = self.optimized_molecule.get_xyz_coordinates
           mlatom_logger.info(f'     {self.job_name:35s}: {self.energy:15.6f}')
           #    print(coordinates)
               
           with open(self.result_xyz_file, "w") as f2:
            f2.write(self.optimized_molecule.get_xyz_string().split('\n')[0] + "\n" +
             "result_" + self.job_name + ':' + str(self.energy) + "\n")
            f2.write('\n'.join(self.optimized_molecule.get_xyz_string().split('\n')[2:]))
            
            sleep(5)

            # f2.write(str(self.number_of_atoms) + "\n")
            # f2.write("result_" + self.job_name + ':' + str(self.energy) + "\n")
            # for i in range(self.number_of_atoms):
            #     f2.write(
            #         " " + "%3s  %10.7f  %10.7f %10.7f\n" % (self.atoms_list[i], coordinates[i][0], coordinates[i][1], coordinates[i][2]))
            # f2.write(self.optimized_molecule.get_xyz_string())
            return
            # return self.optimized_molecule
        except Exception as e:
            print(f"An error occurred: {e}")
            return None

    # def run_frequency_analysis(self):
    #     next_molecule = self.read_molecule_from_file('predict1.xyz')
    #     ml.freq(model=self.aiqm1, molecule=next_molecule, program='ASE')

    # def calculate_thermochemical_properties(self):
    #     next_molecule = self.read_molecule_from_file('predict1.xyz')
    #     ml.thermochemistry(model=self.aiqm1, molecule=next_molecule, program='ASE')

    


    # def run(self):
    #     # self.prepare_input()
    #     self.optimize()
    #     # self.run_frequency_analysis()
    #     # self.calculate_thermochemical_properties()
    # def optimize(self, options):
    #     """
    #     :return:This object will return the optimization status. It will
    #     optimize a structure.
    #     """
    #     self.run()
    #     self.energy = self.optimized_molecule.energy
    #     self.optimized_coordinates = self.optimized_molecule.coordinates
    #     return self.optimized_molecule.energy
    
    # def get_energy(self):
    #     return self.energy
    



def main():
    pass

if __name__ == '__main__':
    main()




# this is another old version of the mlatom_aiqm1.py file
    import os  # noqa: F401
import argparse  # noqa: E402
import logging  # noqa: E402
from pyar.mlatom.data import molecule  # noqa: E402, F811, F401
from pyar.interface import SF, write_xyz, which  # noqa: E402, F401, F811

mlatom_logger = logging.getLogger('pyar.mlatom')

def read_molecule_from_file(filename):
    mol = ml.data.molecule()
    return mol.read_from_xyz_file(filename)

def optimize_geometry(model, initial_molecule):
    return ml.optimize_geometry(model=model, initial_molecule=initial_molecule).optimized_molecule
    

def run_frequency_analysis(model, molecule):  # noqa: F811
    return ml.freq(model=model, molecule=molecule, program='ASE')

def calculate_thermochemical_properties(model, molecule):  # noqa: F811
    return ml.thermochemistry(model=model, molecule=molecule, program='ASE')

def main(input_files):
    for input_file in input_files:
        # Initialize molecule
        init_mol = read_molecule_from_file(input_file)
        # Initialize methods
        aiqm1 = ml.models.methods(method='AIQM1') 
        # Run geometry optimization
        opt_mol = optimize_geometry(aiqm1, init_mol)
        
        # Write optimized geometry to file
        next_mol = opt_mol.get_xyz_string()
        # Run frequency analysis
        run_frequency_analysis(aiqm1, next_mol)
        # Run thermochemical properties calculation
        calculate_thermochemical_properties(aiqm1, next_mol)
    # def print_results(molecule, energy, gradient, std, frequencies, ZPE, Hof):
    #     print(f"Final energy: {energy:.4f} hartree")
    #     print(f"Standard deviation: {std:.4f} hartree")
    #     print(f"Zero-point energy: {ZPE:.4f} hartree")
    #     print(f"Heat of formation: {Hof:.4f} kcal/mol")
    #     print(f"Frequencies (cm^-1): {', '.join([f'{freq:.2f}' for freq in frequencies])}")

    
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run geometry optimization, frequency analysis, and thermochemical properties calculation on a molecule.')
    parser.add_argument('input_file', type=str, help='Input XYZ file path')
    # parser.add_argument('output_file', type=str, help='Output XYZ file path')
    args = parser.parse_args()
    main(args.input_files)
