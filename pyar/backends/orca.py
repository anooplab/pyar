"""
orca.py - backend for the ORCA program

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
from pathlib import Path

import numpy as np

from pyar.backends import SF, require_executable, write_xyz
from pyar.backends.orca_methods import orca_external_method_block, orca_method_keywords
from pyar.backends.subprocess_utils import run_command
from pyar.io.xyz import XYZParseError, parse_xyz

orca_logger = logging.getLogger('pyar.orca')


class Orca(SF):
    """ORCA single-point and optimization wrapper for PyAR."""

    def __init__(self, molecule, qc_params):
        """Prepare an ORCA job from a PyAR molecule and QC settings."""

        super(Orca, self).__init__(molecule)
        self.orca_executable = require_executable("orca", "ORCA")

        self.start_coords = molecule.coordinates
        self.inp_file = 'trial_' + self.job_name + '.inp'
        self.out_file = 'trial_' + self.job_name + '.out'
        self.optimized_coordinates = None
        self.energy = None
        optimization_keyword = {
            'loose': ' LooseOpt',
            'tight': ' TightOpt',
        }.get(qc_params.get('opt_threshold'), ' Opt').strip()
        keyword, is_xtb = orca_method_keywords(qc_params, optimization_keyword)
        if not is_xtb and not qc_params.get('orca_builtin_method') and any(x >= 21 for x in molecule.atomic_number):
            keyword += ' def2-ECP'
        if not is_xtb and (self.multiplicity != 1 or self.scftype in {'uks', 'uhf'}):
            keyword += ' UKS'
        nprocs = qc_params['nprocs']
        keyword += f"\n%pal nprocs {nprocs} end\n"
        if not is_xtb:
            keyword += f"%scf maxiter {qc_params['scf_cycles']} end\n"
        external_method = orca_external_method_block(qc_params)
        if external_method:
            keyword += external_method + "\n"
        if qc_params.get('opt_cycles') is not None:
            keyword += f"%geom MaxIter {int(qc_params['opt_cycles'])} end\n"
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

        self.prepare_input()

        output_xyz = Path(self.inp_file[:-4] + ".xyz")
        # ORCA writes its optimized geometry to the same basename as the input
        # XYZ. Preserve the input under an unambiguous diagnostic name, then
        # require the calculation to produce a fresh output frame.
        if output_xyz.exists():
            input_copy = output_xyz.with_name(output_xyz.stem + ".input.xyz")
            input_copy.write_bytes(output_xyz.read_bytes())
            output_xyz.unlink()

        exit_status = run_command([self.orca_executable, self.inp_file], stdout_path=self.out_file, stderr_path=self.out_file)
        if exit_status == 0:
            output_text = Path(self.out_file).read_text(encoding="utf-8", errors="replace")
            if "****ORCA TERMINATED NORMALLY****" not in output_text:
                orca_logger.error(
                    "orca did not terminate normally job=%s input=%s output=%s cwd=%s",
                    self.job_name, self.inp_file, self.out_file, os.getcwd(),
                )
                return False
            if "THE OPTIMIZATION HAS CONVERGED" not in output_text:
                orca_logger.error("orca did not report geometry convergence for %s", self.job_name)
                return False

            energy = self.get_energy()
            if energy is None or not np.isfinite(energy):
                orca_logger.error("orca returned no finite final energy for %s", self.job_name)
                return False
            try:
                symbols, coordinates, _name, _title, _xyz_energy = parse_xyz(output_xyz)
            except (OSError, XYZParseError) as exc:
                orca_logger.error("orca returned no usable fresh geometry for %s: %s", self.job_name, exc)
                return False
            coordinates = np.asarray(coordinates, dtype=float).reshape((-1, 3))
            if symbols != list(self.atoms_list) or coordinates.shape != (self.number_of_atoms, 3):
                orca_logger.error("orca geometry atom identities/order do not match input for %s", self.job_name)
                return False
            if not np.isfinite(coordinates).all():
                orca_logger.error("orca returned nonfinite coordinates for %s", self.job_name)
                return False

            self.energy = float(energy)
            self.optimized_coordinates = coordinates
            write_xyz(self.atoms_list, self.optimized_coordinates, self.result_xyz_file, energy=self.energy)
            orca_logger.info(
                "orca optimization completed job=%s input=%s output=%s result=%s energy=%s",
                self.job_name, self.inp_file, self.out_file, self.result_xyz_file, self.energy,
            )
            return True
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
                    energy_in_hartrees = None
                if energy_in_hartrees is not None and not np.isfinite(energy_in_hartrees):
                    return None
            return energy_in_hartrees
        except IOError:
            orca_logger.warning("orca output file missing job=%s output=%s", self.job_name, self.out_file)


def main():
    """Module entrypoint retained for parity with the other interfaces."""
    return None


if __name__ == "__main__":
    main()
