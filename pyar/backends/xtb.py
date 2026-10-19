"""Pure xTB geometry optimization backend."""
import logging
import os
from pathlib import Path
import shutil
import subprocess as subp
import sys

import numpy as np

from pyar.backends import SF, require_executable, write_xyz
from pyar.backends.xtb_utils import (
    build_xtb_command, canonical_xtb_model, check_xtb_output, xtb_model_arguments,
)
from pyar.io.xyz import XYZParseError, parse_xyz

xtb_logger = logging.getLogger('pyar.xtb')


def xtb_supports_gxtb(xtb_executable):
    try:
        help_text = subp.check_output(
            [xtb_executable, '--help'],
            stderr=subp.STDOUT,
            text=True,
        )
        return '--gxtb' in help_text
    except Exception as e:
        xtb_logger.warning(f'Could not check xTB --gxtb support: {e}')
        return False


class Xtb(SF):
    """Run a standalone xTB optimization for a single molecule."""

    def __init__(self, molecule, method):
        """Build the xTB command for a molecule and QC parameter set."""
        self.xtb_executable = require_executable('xtb', 'xTB')

        super(Xtb, self).__init__(molecule)

        self.cmd = build_xtb_command(
            self.xtb_executable,
            self.start_xyz_file,
            method,
            opt_threshold=method['opt_threshold'],
        )

        xtb_model = canonical_xtb_model(method.get('xtb_model'))
        if xtb_model == 'gxtb' and not xtb_supports_gxtb(self.xtb_executable):
            raise RuntimeError(
                f"xTB executable {self.xtb_executable!r} does not advertise --gxtb support; "
                "refusing to run a different Hamiltonian as g-xTB"
            )
        self.cmd.extend(xtb_model_arguments(xtb_model))
        xtb_logger.info('Using xTB model %s', xtb_model)

        self.trajectory_xyz_file = 'traj_' + self.job_name + '.xyz'

    def optimize(self, max_cycles=350, gamma=None, restart=False, convergence='normal'):
        """Execute xTB and return a success flag or failure status string."""
        if gamma is not None:
            xtb_logger.error('not implemented in this module. Use xtb_turbo')

        # xTB uses fixed marker and result filenames. Archive them before each
        # launch so artifacts from an earlier attempt cannot authorize this one.
        stale_names = (
            '.xtboptok', '.sccnotconverged', 'NOT_CONVERGED', 'xtbopt.xyz',
            'energy', 'xtb.out', 'xtbopt.log',
        )
        stale_paths = [Path(name) for name in stale_names if Path(name).exists()]
        if stale_paths:
            archive_root = Path('stale_attempts')
            archive_root.mkdir(exist_ok=True)
            index = 0
            while (archive_root / f'run_{index:03d}').exists():
                index += 1
            archive = archive_root / f'run_{index:03d}'
            archive.mkdir()
            for stale_path in stale_paths:
                shutil.move(str(stale_path), str(archive / stale_path.name))

        with open('xtb.out', 'w') as output_file_pointer:
            try:
                out = subp.check_call(
                    self.cmd,
                    stdout=output_file_pointer,
                    stderr=output_file_pointer,
                )
            except Exception as e:
                from pyar.recovery.observations import record_process
                record_process(self.cmd, getattr(e, "returncode", None), exception=e)
                xtb_logger.info('    Optimization failed')
                xtb_logger.error(f'      {e}')
                return False
        from pyar.recovery.observations import record_process
        record_process(self.cmd, out)

        try:
            with open('xtb.out') as output_file_pointer:
                output_text = output_file_pointer.read()
                check_xtb_output(output_text)
        except RuntimeError as exc:
            xtb_logger.error('%s', exc)
            return False

        lowered_output = output_text.lower()
        if any(marker in lowered_output for marker in (
            'scc failed to converge', 'scc convergence failed', 'scc not converged',
        )):
            xtb_logger.error('Current xTB run reported SCC nonconvergence')
            return 'SCFFailed'

        if os.path.isfile('.xtboptok'):
            try:
                symbols, coordinates, _name, _title, _xyz_energy = parse_xyz('xtbopt.xyz')
                coordinates = np.asarray(coordinates, dtype=float).reshape((-1, 3))
                energy = self.energy
                if (symbols != list(self.atoms_list)
                        or coordinates.shape != (self.number_of_atoms, 3)
                        or not np.isfinite(coordinates).all()
                        or energy is None or not np.isfinite(float(energy))):
                    raise ValueError('xTB result has invalid symbols, coordinates, or energy')
            except (OSError, ValueError, TypeError, XYZParseError) as exc:
                xtb_logger.error('Invalid current xTB optimization result: %s', exc)
                return False
            write_xyz(self.atoms_list, coordinates, self.result_xyz_file,
                      job_name=self.job_name, energy=float(energy))

            if os.path.isfile('xtbopt.log'):
                os.rename('xtbopt.log', self.trajectory_xyz_file)

            os.remove('.xtboptok')
            return True

        elif os.path.isfile('.sccnotconverged') or os.path.isfile('NOT_CONVERGED'):
            xtb_logger.info(
                '      SCF Convergence failure in {} run in {}'.format(
                    self.start_xyz_file,
                    os.getcwd(),
                )
            )
            return 'SCFFailed'

        else:
            xtb_logger.info(
                '      Something went wrong with {} run in {}'.format(
                    self.start_xyz_file,
                    os.getcwd(),
                )
            )
            return False

    @property
    def optimized_coordinates(self):
        """Read the optimized xTB coordinates from ``xtbopt.xyz``."""
        return np.loadtxt(
            'xtbopt.xyz',
            dtype=float,
            skiprows=2,
            usecols=(1, 2, 3),
        )

    @property
    def energy(self):
        """Read the final total energy from xTB output files."""
        if os.path.exists('energy'):
            with open('energy') as fp:
                return float(fp.readlines()[-2].split()[1])

        elif os.path.exists('xtb.out'):
            energy = None
            with open('xtb.out') as fp:
                for line in fp.readlines():
                    if 'total E' in line:
                        energy = float(line.split()[-1])
                    if 'TOTAL ENERGY' in line:
                        energy = float(line.split()[3])
            return energy

        else:
            return None


def main():
    pass


if __name__ == '__main__':
    main()
