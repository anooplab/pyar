import torchani
import ase.optimize
from ase.calculators.calculator import Calculator
import logging

logger = logging.getLogger('ani_interface')

class ANICalculationFailed(Exception):
    pass

class ANI(Calculator):

    def __init__(self, species, model='ANI-1x'):
        self.species = species
        try:
            self.model = torchani.models.__dict__[model]()
        except KeyError:
            msg = f"Model '{model}' not found in torchani models."
            logger.error(msg)
            raise ANICalculationFailed(msg)

    def calculate(self, atoms=None, properties=['energy'], system_changes=[]):
        try:
            positions = atoms.get_positions()
            species = atoms.get_chemical_symbols()
        except Exception as e:
            msg = f"Error getting atoms properties: {e}"
            logger.error(msg)
            raise ANICalculationFailed(msg)

        try:
            energy = self.model((species, positions)).energy
        except Exception as e:
            msg = f"ANI model calculation failed: {e}"
            logger.error(msg)
            raise ANICalculationFailed(msg)
            
        self.results = {'energy': energy}

def optimize(atoms):
    try:
        dyn = ase.optimize.BFGS(atoms, trajectory='optimization.traj')
        dyn.run(fmax=0.001)
    except Exception as e:
        msg = f"Geometry optimization failed: {e}"
        logger.error(msg)
        raise ANICalculationFailed(msg)

    try:
        energy = atoms.get_potential_energy()
    except Exception as e:
        msg = f"Getting potential energy failed: {e}"
        logger.error(msg)
        raise ANICalculationFailed(msg)

    return energy

# singlepoint implementation same as before...

class ANIInterface:

    def __init__(self, xyzfile):
        try:
            self.molecules = ase.io.read(xyzfile)
        except OSError:
            msg = f"Could not read XYZ file {xyzfile}"
            logger.error(msg)
            raise ANICalculationFailed(msg)

    def optimize(self):
        try:
            energy = optimize(self.molecules)
        except ANICalculationFailed as e:
            logger.error("Geometry optimization failed")
            raise e
            
        # write optimized geometry
        return energy
    
        
        

    # singlepoint implementation same as before...

# ani_interface.py

import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=["optimize", "singlepoint"])
    parser.add_argument("xyzfile")
    args = parser.parse_args()

    interface = ANIInterface(args.xyzfile)
    
    if args.command == "optimize":
        energy = interface.optimize()
    elif args.command == "singlepoint":
        energy = interface.singlepoint()
        
    print(energy)