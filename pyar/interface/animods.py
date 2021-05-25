
# import logging
# import os
# import subprocess as subp
# # import numpy as np
# import torchani
# from pyar.interface import SF, write_xyz, which
# from ase import Atoms
# from ase.optimize import BFGS

# ani2x_logger = logging.getLogger('pyar.ani2x')


# class ANI2X(SF):
#     def __init__(self, molecule, qc_params, custom_keyword=None):

#         super(ANI2X, self).__init__(molecule)
#         self.asemol = Atoms(self.atoms_list,positions=molecule.coordinates)
#         if qc_params['gamma'] == 0.0:
#             self.asemol.calc = torchani.models.ANI2x().ase()
#         else:
#             self.asemol.calc = torchani.models.ANI2x().ase() #Add restraint energy
#         self.optimized_coordinates = self.asemol.get_positions()
#         self.energy = self.asemol.get_total_energy()
#     def optimize(self, calcdetails):
#         """
#         :returns: True,
#                   False (if CyclesExceeded)
#         """
#         if calcdetails['gamma'] is not None:
#             ani2x_logger.error('Not implemented in this module. Use turbomole')

#         opt = BFGS(self.asemol,trajectory='job_'+self.job_name+'_ase.traj', logfile='job_'+self.job_name+'_opt.log')
        
#         stat = opt.run(fmax=0.001,steps=1000)
#         if stat == False:
#             ani2x_logger.info('    Optimization failed!')
#             return False
#         write_xyz(self.atoms_list, self.asemol.get_positions(), self.result_xyz_file,
#                     job_name=self.job_name,energy=self.asemol.get_total_energy())
#         ani2x_logger.info('    Optimization done!!')
#         self.optimized_coordinates = self.asemol.get_positions()
#         return True


# def main():
#     pass


from pyar import Molecule
import logging
import os
import subprocess as subp
# import numpy as np
import torchani,torch
from pyar.interface import SF, write_xyz, which
from ase import Atoms
from ase.optimize import BFGS
import ase.calculators.calculator
import ase.units
from pyar.afir import restraints
from pyar.data.atomic_data import chemical_symbols

ani2x_logger = logging.getLogger('pyar.ani2x')

class GamCalcAse(ase.calculators.calculator.Calculator):
    """TorchANI calculator for ASE Modified to incorporate Gamma Values"""

    implemented_properties = ['energy', 'forces']

    def __init__(self, species, model, gamma=None,frags = None,overwrite=False):
        super().__init__()
        self.species_to_tensor = torchani.utils.ChemicalSymbolsToInts(species)
        self.model = model
        # Since ANI is used in inference mode, no gradients on model parameters are required here
        for p in self.model.parameters():
            p.requires_grad_(False)
        self.overwrite = overwrite
        self.gamma = gamma
        self.fragments = frags
        a_parameter = next(self.model.parameters())
        self.device = a_parameter.device
        self.dtype = a_parameter.dtype
        try:
            # We assume that the model has a "periodic_table_index" attribute
            # if it doesn't we set the calculator's attribute to false and we
            # assume that species will be correctly transformed by
            # species_to_tensor
            self.periodic_table_index = model.periodic_table_index
        except AttributeError:
            self.periodic_table_index = False

    def calculate(self, atoms=None, properties=['energy'],
                  system_changes=ase.calculators.calculator.all_changes):
        super().calculate(atoms, properties, system_changes)
        cell = torch.tensor(self.atoms.get_cell(complete=True),
                            dtype=self.dtype, device=self.device)
        pbc = torch.tensor(self.atoms.get_pbc(), dtype=torch.bool,
                           device=self.device)
        pbc_enabled = pbc.any().item()

        if self.periodic_table_index:
            species = torch.tensor(self.atoms.get_atomic_numbers(), dtype=torch.long, device=self.device)
        else:
            species = self.species_to_tensor(self.atoms.get_chemical_symbols()).to(self.device)

        species = species.unsqueeze(0)
        coordinates = torch.tensor([self.atoms.get_positions()],requires_grad=True)
        coordinates = coordinates.to(self.device).to(self.dtype)
        energy = self.model((species, coordinates)).energies
        energy *= ase.units.Hartree
        forces = -torch.autograd.grad(energy.squeeze(), coordinates)[0]
        forces = forces.squeeze(0).to('cpu').numpy()
        if self.gamma:
            species_elements = [chemical_symbols[i] for i in species[-1]]
            rst_en , rst_grad = restraints.isotropic(self.fragments,self.atoms.get_chemical_symbols(),coordinates[-1].detach().numpy(),self.gamma)
            energy += rst_en
            # print(forces,rst_grad)
            forces += rst_grad
            # print(forces)
        self.results['energy'] = energy.item()
        self.results['forces'] = forces

class ANI2X(SF):
    def __init__(self, molecule, qc_params):

        super(ANI2X, self).__init__(molecule)
        self.asemol = Atoms(self.atoms_list,positions=molecule.coordinates)
        if qc_params['gamma']:
            self.asemol.calc = GamCalcAse(self.atoms_list,torchani.models.ANI2x(),qc_params['gamma'],molecule.fragments)
        else:
            self.asemol.calc = GamCalcAse(self.atoms_list,torchani.models.ANI2x(),qc_params['gamma'])
        self.optimized_coordinates = self.asemol.get_positions()
        self.energy = self.asemol.get_total_energy()
    def optimize(self, calcdetails):
        """
        :returns: True,
                  False (if CyclesExceeded)
        """
        if calcdetails['gamma'] is not None:
            ani2x_logger.error('Not implemented in this module. Use turbomole')

        opt = BFGS(self.asemol,trajectory='job_'+self.job_name+'_ase.traj', logfile='job_'+self.job_name+'_opt.log')
        
        stat = opt.run(fmax=0.001,steps=1000)
        if stat == False:
            ani2x_logger.info('    Optimization failed!')
            return False
        write_xyz(self.atoms_list, self.asemol.get_positions(), self.result_xyz_file,
                    job_name=self.job_name,energy=self.asemol.get_total_energy())
        ani2x_logger.info('    Optimization done!!')
        self.optimized_coordinates = self.asemol.get_positions()
        return True


def main():
    pass


if __name__ == "__main__":
    main()
