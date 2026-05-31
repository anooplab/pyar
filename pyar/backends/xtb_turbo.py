"""Turbomole/AFIR compatibility backend using xTB for gradients.

The wrapper keeps the AFIR loop separate from the xTB command construction,
which makes the module easier to test and document.
"""

import logging
import os
import subprocess as subp
import sys

import numpy as np

import pyar.backends.turbomole as turbomole
from pyar import backends
from pyar.biases import afir as restraints
from pyar.biases.afir import resolve_gamma
from pyar.data.units import angstrom2bohr, bohr2angstrom
from pyar.backends import SF, require_executable
from pyar.backends.xtb_utils import build_xtb_command

xtb_turbo_logger = logging.getLogger("pyar.xtb_turbo")


class XtbTurbo(SF):
    """Run the AFIR/Turbomole optimization loop with xTB gradients."""

    def __init__(self, molecule, method):
        """Prepare the xTB gradient command, AFIR settings, and Turbomole input state."""
        self.define_command = require_executable("define", "Turbomole")
        self.xtb_executable = require_executable("xtb", "xTB")

        super(XtbTurbo, self).__init__(molecule)

        # Turbomole works in bohr, so normalize the initial coordinates once.
        self.start_coordinates_bohr = angstrom2bohr(np.asarray(molecule.coordinates, dtype=float))
        self.atoms_in_fragments = molecule.fragments
        self.job_directory = "{}/job_{}".format(os.getcwd(), self.job_name)
        self.coordinate_file = "coord"
        self.energy_file = "energy"

        self.gradient_command = build_xtb_command(
            self.xtb_executable,
            "coord",
            {
                **method,
                "charge": self.charge,
                "multiplicity": self.multiplicity,
            },
        )
        self.gradient_command.append("-grad")
        self.gamma = resolve_gamma(method.get("gamma"))

        # Backward-compatible names kept for older call sites.
        self.egrad_program = self.gradient_command
        self.define_executable = self.define_command

        self.energy = None
        self.optimized_coordinates = None

    def optimize(self):
        """Run the AFIR optimization loop until convergence or failure."""
        max_cycles = 250

        turbomole.make_coord(self.atoms_list, self.start_coordinates_bohr, self.coordinate_file)
        turbomole.prepare_control()

        for cycle in range(max_cycles):
            xtb_turbo_logger.debug("optimization cycle %d", cycle)
            status, message, energy, gradients = self.calculate_energy_gradient()
            if status is False:
                xtb_turbo_logger.critical("energy/gradient evaluation failed")
                return "SCFFailed"

            afir_energy, afir_gradient = restraints.isotropic(
                self.atoms_in_fragments,
                self.atoms_list,
                turbomole.get_coords(),
                self.gamma,
            )
            turbomole.rewrite_turbomole_energy_and_gradient_files(
                self.number_of_atoms,
                afir_energy,
                afir_gradient,
            )

            status = turbomole.update_coord()
            if status is False:
                xtb_turbo_logger.critical("coordinate update failed in cycle %d", cycle)
                xtb_turbo_logger.critical("check the job in %s", os.getcwd())
                return "UpdateFailed"

            if turbomole.check_geometry_convergence():
                xtb_turbo_logger.info("converged at cycle %d", cycle)
                self.energy = turbomole.get_energy()
                self.optimized_coordinates = bohr2angstrom(turbomole.get_coords())
                backends.write_xyz(
                    self.atoms_list,
                    self.optimized_coordinates,
                    self.result_xyz_file,
                    job_name=self.job_name,
                    energy=self.energy,
                )
                return True

            with open("energy.dat", "a") as fe:
                fe.writelines("{:3d} {:15.8f} {:15.8f}\n".format(cycle, energy, energy + afir_energy))
        else:
            xtb_turbo_logger.info("cycle limit reached")
            status = "cycle_exceeded"
            self.energy = turbomole.get_energy()
            self.optimized_coordinates = bohr2angstrom(turbomole.get_coords())
            return status

    def calculate_energy_gradient(self):
        """Run one xTB gradient evaluation and return status, message, energy, gradients."""
        with open("job.last", "a") as fp, open("engrad.out", "w") as fc:
            try:
                subp.check_call(self.gradient_command, stderr=fc, stdout=fc)
            except subp.CalledProcessError as exc:
                if exc.output:
                    with open("xtb.out", "w") as fb:
                        fb.write(exc.output)
                msg = "SCF failure. Check files in " + os.getcwd()
                xtb_turbo_logger.error(msg)
                return False, msg, None, None

        with open("engrad.out") as engrad_output:
            message = [line for line in engrad_output.readlines() if "ended" in line]
        if os.path.isfile(".sccnotconverged"):
            msg = "SCF failure. Check files in " + os.getcwd()
            return False, msg, None, None
        if any("abnormally" in line for line in message):
            return False, message, None, None
        return True, message, turbomole.get_energy(), turbomole.get_gradients()

    @property
    def calc_engrad(self):
        """Compatibility alias for older call sites."""
        return self.calculate_energy_gradient()


def main():
    """Module entry point for direct execution."""
    pass


if __name__ == "__main__":
    from pyar.core.molecule import Molecule

    my_mol = Molecule.from_xyz(sys.argv[1])
    geometry = XtbTurbo(my_mol, method={})
    geometry.optimize()
