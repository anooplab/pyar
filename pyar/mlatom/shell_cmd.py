import sys
import os
import time
from mlatom.MLatom import run
from shutil import which


def _resolve_mlatomf_bin():
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
    run()
    time.sleep(1)

def MLatomF():
    os.system(f'{_resolve_mlatomf_bin()} {" ".join(sys.argv[1:])}')
    time.sleep(1)

if __name__ == '__main__':
    run()
