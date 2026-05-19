"""Shell helpers for launching MLatom and the bundled MLatomF binary."""

import sys
import os
import time
from mlatom.MLatom import run
from shutil import which


def _resolve_mlatomf_bin():
    """Return the path to the MLatomF executable used by the shell helper."""
    local_bin = os.path.abspath(os.path.join(os.path.dirname(__file__), 'MLatomF'))
    if os.path.isfile(local_bin) and os.access(local_bin, os.X_OK):
        return local_bin

    installed_bin = which('MLatomF')
    if installed_bin:
        return installed_bin

    raise FileNotFoundError(
        'Unable to locate the MLatomF executable. '
        'Expected it next to pyar.mlatom or on PATH.'
    )

def mlatom():
    """Run the main MLatom Python entrypoint."""
    run()
    time.sleep(1)

def MLatomF():
    """Invoke the MLatomF executable with the current command-line args."""
    os.system(f'{_resolve_mlatomf_bin()} {" ".join(sys.argv[1:])}')
    time.sleep(1)

if __name__ == '__main__':
    run()
