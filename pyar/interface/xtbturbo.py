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
import logging
import os
import subprocess as subp
import sys

import pyar.interface.babel
import pyar.interface.turbomole
from pyar import interface
from pyar.afir import restraints
from pyar.data.units import angstrom2bohr, bohr2angstrom
from pyar.interface import SF

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

        self.start_coords = angstrom2bohr(molecule.coordinates)
        self.atoms_in_fragments = molecule.fragments
        self.job_dir = '{}/job_{}'.format(os.getcwd(), self.job_name)
        self.coord_file = 'coord'
        self.energy_file = 'energy'

        self.egrad_program = ['xtb', 'coord', '-grad']
        if self.charge > 0:
            self.egrad_program += ['-chrg', str(self.charge)]
        if self.multiplicity != 1:
            self.egrad_program += ['-uhf', str(self.multiplicity)]
        self.energy = None
        self.optimized_coordinates = None

    def optimize(self):
        # max_cycles = options['opt_cycles']
        # gamma = options['gamma']
        # convergence = options['opt_threshold']
        max_cycles = 250
        gamma = 100.0


        pyar.interface.turbomole.make_coord(self.atoms_list, self.start_coords, self.coord_file)
        pyar.interface.turbomole.prepare_control()

        for cycle in range(max_cycles):
            xtb_turbo_logger.debug("Optimization Cycle {}".format(cycle))
            # Calculate energy and gradient
            status, message, energy, gradients = self.calc_engrad
            if status is False:
                xtb_turbo_logger.critical('Energy/Gradient evaluation failed')
                return 'SCFFailed'

            # Calculate afir gradient if gamma is greater than zero
            afir_energy, afir_gradient = restraints.isotropic(self.atoms_in_fragments, self.atoms_list,
                                                              pyar.interface.turbomole.get_coords(), gamma)
            pyar.interface.turbomole.rewrite_turbomole_energy_and_gradient_files(self.number_of_atoms, afir_energy,
                                                                                 afir_gradient)

            # Update coordinates and check convergence.
            status = pyar.interface.turbomole.update_coord()
            if status is False:
                xtb_turbo_logger.critical('Coordinate update failed in cycle %d' % cycle)
                xtb_turbo_logger.critical('Check the job in %s' % os.getcwd())
                return 'UpdateFailed'

            convergence_status = pyar.interface.turbomole.check_geometry_convergence()
            if convergence_status is True:
                xtb_turbo_logger.info('converged at {}'.format(cycle))
                self.energy = pyar.interface.turbomole.get_energy()
                self.optimized_coordinates = bohr2angstrom(pyar.interface.turbomole.get_coords())
                interface.write_xyz(self.atoms_list, self.optimized_coordinates,
                                    self.result_xyz_file,
                                    self.job_name,
                                    energy=self.energy)
                return True

            with open('energy.dat', 'a') as fe:
                fe.writelines("{:3d} {:15.8f} {:15.8f}\n".format(cycle, energy, energy + afir_energy))
        else:
            xtb_turbo_logger.info("cycle exceeded")
            status = 'cycle_exceeded'
            self.energy = pyar.interface.turbomole.get_energy()
            self.optimized_coordinates = bohr2angstrom(pyar.interface.turbomole.get_coords())
            return status

    @property
    def calc_engrad(self):
        with open('job.last', 'a') as fp, open('engrad.out', 'w') as fc:

            try:
                subp.check_call(self.egrad_program, stderr=fc, stdout=fc)
            except subp.CalledProcessError as e:
                if e.output:
                    with open('xtb.out', 'w') as fb:
                        fb.write(e.output)
                msg = "SCF Failure. Check files in" + os.getcwd()
                xtb_turbo_logger.error(msg)
                return False, msg, None, None

        msg = [line for line in open('engrad.out').readlines() if 'ended' in line]
        if os.path.isfile('.sccnotconverged'):
            msg = "SCF Failure. Check files in" + os.getcwd()
            return False, msg, None, None
        if 'abnormally' in msg:
            return False, msg, None, None
        else:
            return True, msg, pyar.interface.turbomole.get_energy(), \
                   pyar.interface.turbomole.get_gradients()


def main():
    pass


if __name__ == "__main__":
    from pyar.Molecule import Molecule

    my_mol = Molecule.from_xyz(sys.argv[1])
    geometry = XtbTurbo(my_mol, method={})
    options = {'gamma': 100.0}
    geometry.optimize(options)
