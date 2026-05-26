"""
mopac.py - backend for the MOPAC program

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
import logging

import numpy as np

from pyar import backends
from pyar.backends import SF, require_executable
from pyar.backends.subprocess_utils import run_command

mopac_logger = logging.getLogger("pyar.mopac")


class Mopac(SF):
    """MOPAC optimization wrapper for PyAR."""

    def __init__(self, molecule, qc_params):
        """Prepare a MOPAC job from a PyAR molecule and QC settings."""

        require_executable('mopac', 'MOPAC')
        require_executable('obabel', 'OpenBabel')

        super(Mopac, self).__init__(molecule)

        self.inp_file = 'trial_' + self.job_name + '.mop'
        self.arc_file = 'trial_' + self.job_name + '.arc'
        self.start_coords = molecule.coordinates
        self.optimized_coordinates = []
        self.energy = 0.0
        keyword = f'PM7 PRECISE LET DDMIN=0.0 CYCLES=10000 charge={molecule.charge}'
        self.prepare_input(keyword=keyword)

    def prepare_input(self, keyword=""):
        """Create the MOPAC input file for the current structure."""
        keyword_line = '-xkPM7' if not keyword else '-xk' + keyword
        exit_status = run_command([require_executable("obabel", "OpenBabel"), "-ixyz", self.start_xyz_file, "-omop", keyword_line],
                                  stdout_path=self.inp_file, stderr_path='tmp.log')
        if exit_status == 0:
            os.remove('tmp.log')
        return exit_status

    def optimize(self, max_cycles=350, gamma=0.0, restart=False, convergence='normal'):
        """Optimize the current structure with MOPAC.

        Returns:
            bool: ``True`` when MOPAC finishes normally, otherwise ``False``.
        """

        logfile = "trial_{}.log".format(self.job_name)
        exit_status = run_command([require_executable("mopac", "MOPAC"), self.inp_file], stdout_path=logfile, stderr_path=logfile)
        if exit_status == 0:
            if os.path.exists(self.arc_file):
                self.energy = self.get_energy()
                self.optimized_coordinates = self.get_coords()
                backends.write_xyz(self.atoms_list, self.optimized_coordinates, self.result_xyz_file,
                                    self.job_name, energy=self.energy)
                mopac_logger.info(
                    "mopac optimization completed job=%s input=%s arc=%s log=%s result=%s energy=%s",
                    self.job_name, self.inp_file, self.arc_file, logfile, self.result_xyz_file, self.energy,
                )
                return True
            mopac_logger.error(
                "mopac run completed but arc file missing job=%s input=%s arc=%s log=%s cwd=%s",
                self.job_name, self.inp_file, self.arc_file, logfile, os.getcwd(),
            )
            return False
        mopac_logger.error(
            "mopac command returned non-zero job=%s input=%s log=%s cwd=%s",
            self.job_name, self.inp_file, logfile, os.getcwd(),
        )
        return False

    def get_energy(self):
        """Return the MOPAC heat of formation converted to Hartree."""
        en_kcal = 0.0
        en_kj = 0.0
        try:
            with open(self.arc_file, 'r') as arc_out:
                arc_cont = arc_out.readlines()
            for lines in arc_cont:
                if "HEAT OF FORMATION" in lines:
                    line_cont = lines.split('=')
                    en_kcal = float(line_cont[1].split()[0])
                    en_kj = float(line_cont[2].split()[0])
                    break
        except IOError:
            mopac_logger.warning("mopac arc file missing job=%s arc=%s", self.job_name, self.arc_file)
        return en_kcal / 627.51

    def get_coords(self):
        """Return the optimized Cartesian coordinates from the ARC file."""
        number_of_atoms = None
        with open(self.arc_file) as arc_out:
            arc_cont = arc_out.readlines()
        for lines in arc_cont:
            if 'Empirical Formula:' in lines:
                number_of_atoms = int(lines.split()[-2])
        coordinates = arc_cont[-(number_of_atoms + 1):-1]
        coords = []
        atoms_list = []
        for i in coordinates:
            c = i.split()
            try:
                coords.append(np.array([c[1], c[3], c[5]], dtype=str).astype(float))
            except ValueError:
                return
            atoms_list.append(c[0])
        return np.array(coords)


def main():
    """Module entrypoint retained for parity with the other interfaces."""
    return None


if __name__ == "__main__":
    main()
