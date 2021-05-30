import logging
import torchani,torch
from pyar.interface import SF, write_xyz
from ase import Atoms
from ase.optimize import BFGS
import ase.calculators.calculator
from pyar.afir import restraints

ani_logger = logging.getLogger('pyar.anix')



class GamCalcAse(ase.calculators.calculator.Calculator):
    """TorchANI calculator for ASE Modified to incorporate Gamma Values"""

    implemented_properties = ['energy', 'forces']

    def __init__(self, species, model, gamma=None, frags=None, overwrite=False):
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
        if self.periodic_table_index:
            species = torch.tensor(self.atoms.get_atomic_numbers(), dtype=torch.long, device=self.device)
        else:
            species = self.species_to_tensor(self.atoms.get_chemical_symbols()).to(self.device)

        species = species.unsqueeze(0)
        coordinates = torch.tensor([self.atoms.get_positions()], requires_grad=True)
        coordinates = coordinates.to(self.device).to(self.dtype)
        energy = self.model((species, coordinates)).energies
        energy *= ase.units.Hartree
        forces = -torch.autograd.grad(energy.squeeze(), coordinates)[0]
        forces = forces.squeeze(0).to('cpu').numpy()
        if self.gamma:
            rst_en , rst_grad = restraints.isotropic(self.fragments,self.atoms.get_chemical_symbols(),
                                                        coordinates[-1].detach().numpy(),self.gamma)
            energy += rst_en
            forces += rst_grad
        self.results['energy'] = energy.item()
        self.results['forces'] = forces

class ANI(SF):
    def __init__(self, molecule,qc_params):
        super(ANI, self).__init__(molecule)
        self.asemol = Atoms(self.atoms_list,positions=molecule.coordinates)
        if qc_params['software'] == 'ani2x':
            self.animodel = torchani.models.ANI2x()
        if qc_params['software'] == 'ani1x':
            self.animodel = torchani.models.ANI1x()
        if qc_params['software'] == 'ani1cc':
            self.animodel = torchani.models.ANI1ccx()
        if qc_params['gamma']:
            self.asemol.calc = GamCalcAse(self.animodel.species,self.animodel,
                                          qc_params['gamma'],molecule.fragments)
        else:
            self.asemol.calc = self.animodel.ase()
        self.optimized_coordinates = self.asemol.get_positions()
        self.energy = self.asemol.get_total_energy()

    def optimize(self, calcdetails):
        """
        :returns: True,
                  False (if CyclesExceeded)
        """

        opt = BFGS(self.asemol, trajectory='job_' + self.job_name + '_ase.traj',
                   logfile='job_' + self.job_name + '_opt.log')

        stat = opt.run(fmax=0.001, steps=1000)
        if not stat:
            ani_logger.info('    Optimization failed!')
            return False
        write_xyz(self.atoms_list, self.asemol.get_positions(), self.result_xyz_file,
                  job_name=self.job_name, energy=self.asemol.get_total_energy())
        ani_logger.info('    Optimization done!!')
        self.optimized_coordinates = self.asemol.get_positions()
        return True


def main():
    pass


if __name__ == "__main__":
    main()
