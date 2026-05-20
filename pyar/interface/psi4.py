"""
psi4.py - interface to PSI4 program

Copyright (C) 2016 by Anoop Ayyappan
Indian Institute of Technology Kharagpur, India

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

import numpy as np

from pyar.interface import SF, require_executable, write_xyz
from pyar.interface.subprocess_utils import run_command

psi4_logger = logging.getLogger("pyar.psi4")


class Psi4(SF):
    """Psi4 optimization wrapper for PyAR."""

    def __init__(self, molecule, method, custom_keyword=None):
        """Prepare a Psi4 job from a PyAR molecule and QC settings."""

        super(Psi4, self).__init__(molecule)
        self.psi4_executable = require_executable("psi4", "Psi4")

        self.start_coords = molecule.coordinates
        self.inp_file = 'trial_' + self.job_name + '.in'
        self.out_file = 'trial_' + self.job_name + '.out'
        self.optimized_coordinates = []
        self.number_of_atoms = len(self.atoms_list)
        self.energy = 0.0
        keyword = "set basis def2-SVP\n"
        if any(x >= 21 for x in molecule.atomic_number):
            keyword += ' def2-ECP'
        if custom_keyword is not None:
            keyword += custom_keyword
        self.prepare_input(keyword=keyword)

    def prepare_input(self, keyword=""):
        """Write the Psi4 input file for the current molecule."""
        coords = self.start_coords
        if self.scftype == 'uks':
            keyword += 'UKS'
        with open(self.inp_file, "w") as f1:
            f1.write(keyword + "\n")
            f1.write("molecule {\n")
            f1.write("{0} {1}\n".format(str(self.charge), str(self.multiplicity)))
            for i in range(self.number_of_atoms):
                f1.write(
                    " " + "%3s  %10.7f  %10.7f %10.7f\n" % (self.atoms_list[i], coords[i][0], coords[i][1], coords[i][2]))
            f1.write("}\n")
            f1.write("optimize(\"B97-D\")")

    def optimize(self, max_cycles=350, gamma=0.0, restart=False, convergence='normal'):
        """
        Optimize the current structure with Psi4.

        Returns:
            bool: ``True`` when Psi4 finishes normally, otherwise ``False``.
        """

        exit_status = run_command([self.psi4_executable, self.inp_file], stdout_path=self.out_file, stderr_path=self.out_file)
        if exit_status == 0:
            with open(self.out_file, "r") as f:
                output = f.read()
            if "  **** Optimization is complete!" in output:
                self.energy = self.get_energy()
                self.optimized_coordinates = self.get_coordinates()
                write_xyz(self.atoms_list,
                          self.optimized_coordinates,
                          self.result_xyz_file, energy=self.energy)
                psi4_logger.info(
                    "psi4 optimization completed job=%s input=%s output=%s result=%s energy=%s",
                    self.job_name, self.inp_file, self.out_file, self.result_xyz_file, self.energy,
                )
                return True
            psi4_logger.error(
                "psi4 optimization did not terminate normally job=%s input=%s output=%s cwd=%s",
                self.job_name, self.inp_file, self.out_file, os.getcwd(),
            )
            return False
        psi4_logger.error(
            "psi4 command returned non-zero job=%s input=%s output=%s cwd=%s",
            self.job_name, self.inp_file, self.out_file, os.getcwd(),
        )
        return False

    def get_energy(self):
        """
        Return the final Psi4 energy in Hartree.
        """
        try:
            energy_in_hartrees = None
            with open(self.out_file, "r") as out:
                lines = out.readlines()
                for line in lines:
                    if "Final energy is" in line:
                        energy_in_hartrees = float((line.strip().split())[-1])
            return energy_in_hartrees
        except IOError:
            psi4_logger.warning("psi4 output file missing job=%s output=%s", self.job_name, self.out_file)

    def get_coordinates(self):
        """
        Return the optimized Cartesian coordinates from the Psi4 log.
        """
        try:
            with open(self.out_file, "r") as out:
                for line in out:
                    if "Final optimized geometry and variables:" in line:
                        break
                coords = []
                for line in out.readlines()[5:5 + self.number_of_atoms + 1]:
                    lc = line.split()
                    if len(lc) == 4:
                        try:
                            _, x, y, z = lc
                            coords.append([float(x), float(y), float(z)])
                        except Exception as e:
                            psi4_logger.exception(
                                "psi4 coordinate parse failed job=%s output=%s line=%s",
                                self.job_name, self.out_file, line.rstrip(),
                            )
            return np.array(coords)
        except IOError:
            psi4_logger.warning("psi4 output file missing job=%s output=%s", self.job_name, self.out_file)


def main():
    """Run the Psi4 workflow from the command line."""
    from pyar.Molecule import Molecule
    import sys
    mol = Molecule.from_xyz(sys.argv[1])
    method = {'charge': 0, 'multiplicity': 1, 'scftype': 'rhf'}
    geometry = Psi4(mol, method)
    geometry.optimize()


if __name__ == "__main__":
    main()
