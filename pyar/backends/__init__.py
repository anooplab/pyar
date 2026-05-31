"""Public backend helpers for PyAR 2.0."""

import os

from pyar.backend_capabilities import (
    BackendCapabilities,
    backend_family,
    backend_supports_geometry_optimization,
    backend_supports_staged_optimization,
    get_backend_capabilities,
    register_backend_capabilities,
    supported_geometry_backends,
    unsupported_qc_options,
)
from .subprocess_utils import run_command, run_output


def which(program):
    """Return the absolute path to an executable if it exists on PATH."""
    def is_exe(exec_path):
        return os.path.isfile(exec_path) and os.access(exec_path, os.X_OK)

    file_path, _ = os.path.split(program)
    if file_path:
        if is_exe(program):
            return program
    else:
        for path in os.environ.get("PATH", "").split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


def require_executable(program, friendly_name=None):
    """Return an executable path or raise a clear FileNotFoundError."""
    executable = which(program)
    if executable is None:
        name = friendly_name or program
        raise FileNotFoundError(f"{name} executable '{program}' was not found on PATH")
    return executable


# from pyar.backends.mlatom_aiqm1 import MlatomAiqm1

class SF(object):
    """Base state for PyAR interface workflows that write XYZ files."""

    def __init__(self, molecule):  # noqa: F811
        """Initialize a workflow helper from a PyAR molecule object."""
        self.job_name = molecule.name
        self.start_xyz_file = 'trial_' + self.job_name + '.xyz'
        self.result_xyz_file = 'result_' + self.job_name + '.xyz'
        self.atoms_list = molecule.atoms_list
        self.number_of_atoms = molecule.number_of_atoms
        self.title = 'opt {}'.format(molecule.title)
        self.charge = molecule.charge
        self.multiplicity = molecule.multiplicity
        self.scftype = molecule.scftype
        if not os.path.isfile(self.start_xyz_file):
            write_xyz(self.atoms_list,
                      molecule.coordinates, self.start_xyz_file,
                      job_name=self.job_name)


def write_xyz(atoms_list, coordinates, filename, job_name='no_name', energy=0.0):
    """Write a simple XYZ file with an optional energy label."""
    with open(filename, 'w') as fp:
        fp.write("%3d\n" % len(coordinates))
        fp.write(job_name + ':' + str(energy) + '\n')
        for a, c in zip(atoms_list, coordinates):
            fp.write("{:<2}{:12.5f}{:12.5f}{:12.5f}\n".format(a, c[0], c[1], c[2]))


__all__ = [
    "BackendCapabilities",
    "SF",
    "backend_family",
    "backend_supports_geometry_optimization",
    "backend_supports_staged_optimization",
    "get_backend_capabilities",
    "register_backend_capabilities",
    "require_executable",
    "supported_geometry_backends",
    "unsupported_qc_options",
    "which",
    "write_xyz",
    "run_command",
    "run_output",
]
