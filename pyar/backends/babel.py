"""
babel.py - OpenBabel backend helpers for PyAR

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

import os
import logging
import sys

import numpy as np

from pyar.backends import SF, require_executable
from pyar.backends.subprocess_utils import run_command, run_output

babel_logger = logging.getLogger("pyar.babel")


def _obabel_executable():
    """Return the OpenBabel 3 command used for format conversion."""
    return require_executable("obabel", "OpenBabel")


class OBabel(SF):
    """OpenBabel-backed XYZ optimization helper for PyAR."""

    def __init__(self, molecule, forcefield=None):
        """Prepare an OpenBabel workflow helper for a PyAR molecule."""

        _obabel_executable()

        super(OBabel, self).__init__(molecule)

        if not os.path.isfile(self.start_xyz_file):
            molecule.mol_to_xyz(self.start_xyz_file)
        self.optimized_coordinates = []
        self.energy = 0.0

    def optimize(self, max_cycles=350, gamma=0.0, restart=False, convergence='normal'):
        """Optimize the current structure with OpenBabel."""

        exit_status = run_command(
            [require_executable("obminimize", "OpenBabel"), "-ff", "uff", '-n', max_cycles, self.start_xyz_file],
            stdout_path='tmp.xyz',
            stderr_path='tmp.log',
        )
        if exit_status != 0:
            babel_logger.error(
                "obabel optimization failed job=%s input=%s log=%s cwd=%s",
                self.job_name, self.start_xyz_file, 'tmp.log', os.getcwd(),
            )
            return False
        with open('tmp.xyz') as xyzfile, open(self.result_xyz_file, 'w') as result_xyz_file:
            for line in xyzfile:
                if 'WARNING' not in line:
                    result_xyz_file.write(line)
        babel_logger.info(
            "obabel optimization completed job=%s input=%s result=%s tmp_xyz=%s",
            self.job_name, self.start_xyz_file, self.result_xyz_file, 'tmp.xyz',
        )
        self.energy = self.get_energy()
        self.optimized_coordinates = self.get_coords()
        return True

    def get_coords(self):
        """Return the optimized Cartesian coordinates from the XYZ file."""
        return np.loadtxt(self.result_xyz_file, dtype=float, skiprows=2, usecols=(1, 2, 3))

    def get_energy(self):
        """Return the OpenBabel energy in Hartree if the calculation completed."""
        energy_path = self.job_name + '.ene'
        exit_status = run_command([require_executable("obenergy", "OpenBabel"), "-ff", "uff", self.result_xyz_file],
                                  stdout_path=energy_path, stderr_path=energy_path)
        if exit_status == 0:
            with open(self.job_name + '.ene', 'r') as energy_file:
                return float(energy_file.readlines()[-1].split()[3])
        babel_logger.error(
            "obabel energy command failed job=%s input=%s log=%s cwd=%s",
            self.job_name, self.result_xyz_file, energy_path, os.getcwd(),
        )
        return None


def xyz_to_mopac_input(xyzfile, mopac_input_file, keyword=None):
    """Convert an XYZ file to a Mopac input file with OpenBabel."""
    keyword_line = '-xkPM7' if keyword is None else '-xk' + keyword
    run_command([_obabel_executable(), "-ixyz", xyzfile, "-omop", mopac_input_file, keyword_line],
                stdout_path='tmp.log', stderr_path='tmp.log')
    with open('tmp.log') as log_file:
        first_line = log_file.readline().strip()
    babel_logger.info("obabel mopac conversion xyz=%s output=%s note=%s", xyzfile, mopac_input_file, first_line)
    os.remove('tmp.log')


def xyz_to_sdf_file(xyz_input_files, sdf_output_file):
    """Convert one or more XYZ files to a single SDF file."""
    run_command([_obabel_executable(), "-ixyz"] + xyz_input_files + ["-osdf", sdf_output_file],
                stdout_path='tmp.log', stderr_path='tmp.log')
    babel_logger.info("obabel sdf conversion inputs=%s output=%s", xyz_input_files, sdf_output_file)
    os.remove('tmp.log')


def make_inchi_string_from_xyz(xyzfile):
    """Return an InChI string for a molecule stored in an XYZ file."""
    if os.path.isfile(xyzfile):
        inchi = run_output([_obabel_executable(), "-ixyz", str(xyzfile), "-oinchi"], stderr_path='OBabel.log')
        babel_logger.info("obabel inchi conversion input=%s log=%s", xyzfile, 'OBabel.log')
        return inchi.decode("utf-8").strip()
    else:
        raise IOError("file %s does not exists" % xyzfile)


def make_smile_string_from_xyz(xyzfile):
    """Return a SMILES string for a molecule stored in an XYZ file."""
    if os.path.isfile(xyzfile):
        with open('OBabel.log', 'w') as ferr:
            try:
                pre_smile = run_output([_obabel_executable(), "-ixyz", str(xyzfile), "-osmi", "-xn"], stderr_path='OBabel.log')
                smile = pre_smile.decode("utf-8").strip()
            except Exception as e:
                ferr.write(str(e))
                babel_logger.exception("obabel smiles conversion failed input=%s log=%s", xyzfile, 'OBabel.log')
                smile = ''
            else:
                babel_logger.info("obabel smiles conversion input=%s log=%s", xyzfile, 'OBabel.log')
            return smile
    else:
        raise IOError("file %s does not exists" % xyzfile)


def main(input_files):
    """Run the OpenBabel optimization workflow for one or more XYZ files."""
    from pyar.Molecule import Molecule
    for f in input_files:
        mol = Molecule.from_xyz(f)
        g = OBabel(mol)
        g.optimize()


if __name__ == "__main__":
    main(sys.argv[1:])
