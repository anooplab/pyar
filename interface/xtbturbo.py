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
import shutil
import subprocess as subp
import sys

import numpy as np

import file_manager
from afir import restraints
import interface.babel
import interface.turbomole


class XtbTurbo(object):
    def __init__(self, molecule, charge=0, multiplicity=1, scftype='rhf'):
        self.job_name = molecule.name
        self.charge = charge
        self.multiplicity = multiplicity
        self.scftype = scftype
        self.atoms_list = molecule.atoms_list
        self.start_xyz_file = 'trial_' + self.job_name + '.xyz'
        self.result_xyz_file = 'result_' + self.job_name + '.xyz'
        self.number_of_atoms = molecule.number_of_atoms
        self.job_dir = '{}/job_{}'.format(os.getcwd(), self.job_name)
        self.start_coords = molecule.coordinates
        self.coord_file = '{}/coord'.format(self.job_dir)
        self.energy_file = '{}/energy'.format(self.job_dir)
        self.egrad_program = ['xtb', 'coord', '-grad']
        if charge > 0:
            self.egrad_program = self.egrad_program + ['-chrg', str(charge)]
        if multiplicity != 1:
            self.egrad_program = self.egrad_program + ['-uhf', str(multiplicity)]
        self.energy = None
        self.optimized_coordinates = None

    def optimize(self, max_cycles=1000, gamma=0.0):

        cwd = os.getcwd()
        job_dir = 'job_' + self.job_name
        file_manager.make_directories(job_dir)

        if gamma > 0.0:
            if os.path.isfile('fragment'):
                shutil.copy('fragment', job_dir)
            else:
                print("file, {}, not found".format('fragment'))
                sys.exit()

        os.chdir(job_dir)

        interface.turbomole.make_coord(self.atoms_list, self.start_coords, self.coord_file)
        interface.turbomole.prepare_control()

        for cycle in range(max_cycles):
            # Calculate energy and gradient
            status, message, energy, gradients = self.calc_engrad
            if status is False:
                print('Energy/Gradient evaluation failed')
                os.chdir(cwd)
                return False

            # Calculate afir gradient if gamma is greater than zero
            if gamma > 0.0:
                re, trg, rg = restraints.isotropic(gamma)
                interface.turbomole.rewrite_turbomole_energy_and_gradient_files(self.number_of_atoms, rg, re, trg)

            # Update coordinates and check convergence.
            status, msg = interface.turbomole.update_coord()
            if status is True:
                print('converged at', cycle)
                self.energy = interface.turbomole.get_energy()
                self.optimized_coordinates = interface.turbomole.get_coords()
                interface.babel.write_xyz(self.atoms_list, self.optimized_coordinates,
                                          self.result_xyz_file,
                                          self.job_name,
                                          energy=self.energy)
                shutil.copy(self.result_xyz_file, cwd)
                os.chdir(cwd)
                return True
            else:
                with open('energy.dat','a') as fe:
                    if gamma > 0.0:
                        fe.writelines("{:3d} {:15.8f} {:15.8f}\n".format(cycle, energy, energy+re))
                    else:
                        fe.writelines("{:3d} {:15.8f}\n".format(cycle, energy))
        else:
            print("cycle exceeded")
            status = 'cycle_exceeded'
        os.chdir(cwd)

    @property
    def calc_engrad(self):
        with open('job.last', 'a') as fp, open('engrad.out', 'w') as fc:
            subp.check_call(self.egrad_program, stdout=fp, stderr=fc)
        msg = [line for line in open('engrad.out').readlines() if 'ended' in line]
        if os.path.isfile('.sccnotconverged'):
            msg = "SCF Failure. Check files in"+os.getcwd()
            return False, msg, None
        if 'abnormally' in msg:
            return False, msg, None
        else:
            return True, msg, interface.turbomole.get_energy(), \
                   interface.turbomole.get_gradients(self.number_of_atoms)


def main():
    from Molecule import Molecule
    mol = Molecule.from_xyz(sys.argv[1])
    geometry = XtbTurbo(mol)
    if len(sys.argv) == 2:
        geometry.optimize()
    elif len(sys.argv) == 3:
        gamma_force = sys.argv[2]
        geometry.optimize(gamma=gamma_force)
    else:
        usage()
        sys.exit()
    return


if __name__ == "__main__":
    from Molecule import Molecule
    mol = Molecule.from_xyz(sys.argv[1])
    geometry = XtbTurbo(mol)
    geometry.optimize(gamma=100.0)
