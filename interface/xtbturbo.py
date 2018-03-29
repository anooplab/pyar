"""
turbomole.py - interface to turbomole program

Copyright (C) 2016 by Surajit Nandi, Anoop Ayyappan, and Mark P. Waller
Indian Institute of Technology Kharagpur, India and Westfaelische Wilhelms
Universitaet Muenster, Germany

This file is part of the PyAR project.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
"""
import os
import subprocess as subp
import sys

import interface
import interface.babel
import interface.turbomole
from afir import restraints
from interface import SF
from units import angstrom2bohr, bohr2angstrom

import logging
xtb_turbo_logger = logging.getLogger('pyar.xtbturbo')

class XtbTurbo(SF):

    def __init__(self, molecule, method):

        if interface.which('define') is None:
            print('set Turbomole path')
            sys.exit()
        if interface.which('xtb') is None:
            print('set XTB path')
            sys.exit()

        super(XtbTurbo, self).__init__(molecule)

        self.charge = method['charge']
        self.scf_type = method['scftype']
        self.multiplicity = method['multiplicity']
        self.start_coords = angstrom2bohr(molecule.coordinates)
        self.atoms_in_fragments = molecule.fragments
        self.job_dir = '{}/job_{}'.format(os.getcwd(), self.job_name)
        self.coord_file = 'coord'
        self.energy_file = 'energy'

        self.egrad_program = ['xtb', 'coord', '-grad']
        if self.charge > 0:
            self.egrad_program = self.egrad_program + ['-chrg', str(self.charge)]
        if self.multiplicity != 1:
            self.egrad_program = self.egrad_program + ['-uhf', str(self.multiplicity)]
        self.energy = None
        self.optimized_coordinates = None

    def optimize(self, max_cycles=350, gamma=0.0):

        interface.turbomole.make_coord(self.atoms_list, self.start_coords, self.coord_file)
        interface.turbomole.prepare_control()

        for cycle in range(max_cycles):
            xtb_turbo_logger.debug("Optimization Cycle {}".format(cycle))
            # Calculate energy and gradient
            status, message, energy, gradients = self.calc_engrad
            if status is False:
                xtb_turbo_logger.critical('Energy/Gradient evaluation failed')
                return False

            # Calculate afir gradient if gamma is greater than zero
            if gamma > 0.0:
                augmented_energy, trg, rg = restraints.isotropic(self.atoms_in_fragments,self.atoms_list, interface.turbomole.get_coords(), gamma)
                interface.turbomole.rewrite_turbomole_energy_and_gradient_files(self.number_of_atoms, rg, augmented_energy, trg)

            # Update coordinates and check convergence.
            status = interface.turbomole.update_coord()
            if status is True:
                xtb_turbo_logger.info('converged at {}'.format(cycle))
                self.energy = interface.turbomole.get_energy()
                self.optimized_coordinates = bohr2angstrom(interface.turbomole.get_coords())
                interface.write_xyz(self.atoms_list, self.optimized_coordinates,
                                    self.result_xyz_file,
                                    self.job_name,
                                    energy=self.energy)
                return True
            with open('energy.dat','a') as fe:
                if gamma > 0.0:
                    fe.writelines("{:3d} {:15.8f} {:15.8f}\n".format(cycle, energy, energy+augmented_energy))
                else:
                    fe.writelines("{:3d} {:15.8f}\n".format(cycle, energy))
        else:
            xtb_turbo_logger.info("cycle exceeded")
            status = 'cycle_exceeded'
            self.energy = interface.turbomole.get_energy()
            self.optimized_coordinates = bohr2angstrom(interface.turbomole.get_coords())
            return status

    @property
    def calc_engrad(self):
        with open('job.last', 'a') as fp, open('engrad.out', 'w') as fc:
            subp.check_call(self.egrad_program, stdout=fp, stderr=fc)
        msg = [line for line in open('engrad.out').readlines() if 'ended' in line]
        if os.path.isfile('.sccnotconverged'):
            msg = "SCF Failure. Check files in"+os.getcwd()
            return False, msg, None, None
        if 'abnormally' in msg:
            return False, msg, None, None
        else:
            return True, msg, interface.turbomole.get_energy(), \
                   interface.turbomole.get_gradients(self.number_of_atoms)


def main():
    pass


if __name__ == "__main__":
    from Molecule import Molecule
    my_mol = Molecule.from_xyz(sys.argv[1])
    geometry = XtbTurbo(my_mol, method={})
    geometry.optimize(gamma=100.0)
