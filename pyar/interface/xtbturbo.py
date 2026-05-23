"""Turbomole/AFIR workflow backed by xTB gradients.

This module runs the AFIR/Turbomole optimization loop using xTB for the
energy and gradient evaluations.
"""

import logging
import os
import subprocess as subp
import sys

import numpy as np

import pyar.interface.turbomole
from pyar import interface
from pyar.afir import restraints
from pyar.data.units import angstrom2bohr, bohr2angstrom
from pyar.interface import SF, require_executable
from pyar.interface.xtb_utils import build_xtb_command

xtb_turbo_logger = logging.getLogger('pyar.xtbturbo')


class XtbTurbo(SF):
    """Run the AFIR/Turbomole optimization loop with xTB gradients."""

    def __init__(self, molecule, method):
        """Build the xTB gradient evaluator used by the AFIR loop."""

        self.define_executable = require_executable('define', 'Turbomole')
        self.xtb_executable = require_executable('xtb', 'xTB')

        super(XtbTurbo, self).__init__(molecule)

        self.start_coords = angstrom2bohr(np.asarray(molecule.coordinates, dtype=float))
        self.atoms_in_fragments = molecule.fragments
        self.job_dir = '{}/job_{}'.format(os.getcwd(), self.job_name)
        self.coord_file = 'coord'
        self.energy_file = 'energy'

        self.egrad_program = build_xtb_command(
            self.xtb_executable,
            'coord',
            {
                **method,
                "charge": self.charge,
                "multiplicity": self.multiplicity,
            },
        )
        self.egrad_program.append('-grad')
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
