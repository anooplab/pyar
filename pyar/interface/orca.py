"""
orca.py - interface to ORCA program

Copyright (C) 2016 by Surajit Nandi, Anoop Ayyappan, and Mark P. Waller
Indian Institute of Technology Kharagpur, India and Westfaelische Wilhelms
Universitaet Muenster, Germany

This file is part of the pyar project.

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

from pyar.interface import SF, write_xyz, which
from pyar.interface.subprocess_utils import run_command

orca_logger = logging.getLogger('pyar.orca')


class Orca(SF):
    """ORCA single-point and optimization wrapper for PyAR."""

    def __init__(self, molecule, qc_params):
        """Prepare an ORCA job from a PyAR molecule and QC settings."""

        super(Orca, self).__init__(molecule)

        self.start_coords = molecule.coordinates
        self.inp_file = 'trial_' + self.job_name + '.inp'
        self.out_file = 'trial_' + self.job_name + '.out'
        self.optimized_coordinates = []
        self.energy = 0.0
        # print(custom_keyword)
        keyword = f"! {qc_params['method']} {qc_params['basis']}"

        if any(x >= 21 for x in molecule.atomic_number):
            keyword += ' def2-ECP'
        keyword += ' RI def2/J D3BJ KDIIS'
        if self.scftype == 'uks':
            keyword += ' UKS'
        nprocs = qc_params['nprocs']
        # if custom_keyword is not None:
        #     keyword += custom_keyword
        keyword += f"\n%pal nprocs {nprocs} end\n"
        keyword += f"%scf maxiter {qc_params['scf_cycles']} end\n"
        self.keyword = keyword

    def prepare_input(self):
        """Write the ORCA input file for the current molecule."""
        keyword = self.keyword
        coords = self.start_coords
        with open(self.inp_file, "w") as f1:
            f1.write(keyword + "\n")
            f1.write("*xyz {0} {1}\n".format(str(self.charge), str(self.multiplicity)))
            for i in range(self.number_of_atoms):
                f1.write(
                    " " + "%3s  %10.7f  %10.7f %10.7f\n" % (self.atoms_list[i], coords[i][0], coords[i][1], coords[i][2]))
            f1.write("*")

    def optimize(self):
        """
        Optimize the current structure with ORCA.

        Returns:
            bool: ``True`` when ORCA finishes normally, otherwise ``False``.
        """

        self.keyword = self.keyword + '!Opt'
        self.prepare_input()

        exit_status = run_command([which("orca"), self.inp_file], stdout_path=self.out_file, stderr_path=self.out_file)
        if exit_status == 0:
            with open(self.out_file, "r") as f:
                line = f.readlines()
            if "****ORCA TERMINATED NORMALLY****" in line[-2]:
                self.energy = self.get_energy()
                self.optimized_coordinates = np.loadtxt(self.inp_file[:-4] + ".xyz", dtype=float, skiprows=2,
                                                        usecols=(1, 2, 3))
                write_xyz(self.atoms_list,
                          self.optimized_coordinates,
                          self.result_xyz_file, energy=self.energy)
                orca_logger.info(
                    "orca optimization completed job=%s input=%s output=%s result=%s energy=%s",
                    self.job_name, self.inp_file, self.out_file, self.result_xyz_file, self.energy,
                )
                return True
            orca_logger.error(
                "orca optimization did not terminate normally job=%s input=%s output=%s cwd=%s",
                self.job_name, self.inp_file, self.out_file, os.getcwd(),
            )
            return False
        orca_logger.error(
            "orca command returned non-zero job=%s input=%s output=%s cwd=%s",
            self.job_name, self.inp_file, self.out_file, os.getcwd(),
        )
        return False

    def get_energy(self):
        """
        Return the final single-point energy in Hartree.
        """
        try:
            with open(self.out_file, "r") as out:
                line = out.readlines()
                en_steps = [item for item in line if
                            "FINAL SINGLE POINT ENERGY" in item]
                if en_steps:
                    energy_in_hartrees = float((en_steps[-1].strip().split())[-1])
                else:
                    energy_in_hartrees = 0.0
            return energy_in_hartrees
        except IOError:
            orca_logger.warning("orca output file missing job=%s output=%s", self.job_name, self.out_file)


def main():
    """Module entrypoint retained for parity with the other interfaces."""
    return None


if __name__ == "__main__":
    main()
