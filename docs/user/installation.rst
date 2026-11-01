Install PyAR
============

Choose an environment
---------------------

PyAR requires Python 3.10 or newer. Use an isolated Python environment. Python 3.10–3.13 is the classified,
CI-tested range. Python 3.11 is a convenient starting point. Linux is the primary
command-line environment; on Windows use a Linux environment such as WSL.
Commands below use a POSIX shell. Installation on your own macOS/WSL environment
must still be checked with ``pyar doctor`` and the first tutorial.

These pages describe PyAR 1.6. This development version is installed from the
repository until an approved package release is available. The distribution is
called ``pyar-chem``; the command and Python package are called ``pyar``.

.. code-block:: bash

   python3.11 -m venv .venv
   source .venv/bin/activate
   python -m pip install --upgrade pip
   python -m pip install "pyar-chem[selection,xtb,conformer] @ git+https://github.com/anooplab/pyar.git@master"
   pyar --version
   pyar --help

For a reproducible environment replace ``master`` with the commit you intend
to use and retain the resolved package versions. Avoid installing into your
system Python. A plain ``pip install pyar-chem`` installs the currently published
release, which may differ from the version described here.

Install a calculation backend
------------------------------

The Python extra and the external executable are separate requirements. The
``xtb`` extra does not install the ``xtb`` executable. Install xTB using its
supported platform instructions, then put the executable on ``PATH``:

.. code-block:: bash

   xtb --version
   pyar install status
   pyar backends
   pyar doctor

Read :doc:`/reference/backends` for backend-specific prerequisites. Start with
GFN2-xTB for the small examples here. An ORCA installation, for example, must
satisfy ORCA's own distribution/license requirements and provide its executable.
You do not need every backend to use PyAR.

Check the selected task
-----------------------

Download the water input from :doc:`first_calculation` and validate without
starting a calculation:

.. code-block:: bash

   pyar energy water.xyz --backend xtb --check

A successful check confirms the requested route is ready. It is not a numerical
calculation. Missing packages/programs should be corrected in this same activated
environment. Check the executable named in the diagnostic instead of installing
unrelated extras. ``pyar install status --json`` records discovered capabilities.

Removing an old installation
-----------------------------

Use a fresh environment for 1.6. Stale launcher files from unrelated installations
can remain on your shell's PATH; check which environment supplies ``pyar``.
The supported 1.6 package installs a single executable. See
:doc:`/reference/upgrade` before using older calculation directories.

Optional Python capabilities
----------------------------

* ``aggregate``: install the declared `aggregate` Python capability set; external executable/model prerequisites still apply.
* ``grow``: install the declared `grow` Python capability set; external executable/model prerequisites still apply.
* ``conformer``: install the declared `conformer` Python capability set; external executable/model prerequisites still apply.
* ``reaction``: install the declared `reaction` Python capability set; external executable/model prerequisites still apply.
* ``microsolvate``: install the declared `microsolvate` Python capability set; external executable/model prerequisites still apply.
* ``neb``: install the declared `neb` Python capability set; external executable/model prerequisites still apply.
* ``geometric``: install the declared `geometric` Python capability set; external executable/model prerequisites still apply.
* ``aimnet2``: install the declared `aimnet2` Python capability set; external executable/model prerequisites still apply.
* ``sella``: install the declared `sella` Python capability set; external executable/model prerequisites still apply.
* ``geodesic``: install the declared `geodesic` Python capability set; external executable/model prerequisites still apply.
* ``structure-comparison``: install the declared `structure-comparison` Python capability set; external executable/model prerequisites still apply.
* ``identity``: install the declared `identity` Python capability set; external executable/model prerequisites still apply.
* ``analysis``: install the declared `analysis` Python capability set; external executable/model prerequisites still apply.
* ``ml``: install the declared `ml` Python capability set; external executable/model prerequisites still apply.
* ``xtb``: install the declared `xtb` Python capability set; external executable/model prerequisites still apply.
* ``selection``: install the declared `selection` Python capability set; external executable/model prerequisites still apply.
* ``all``: install the declared `all` Python capability set; external executable/model prerequisites still apply.

Source-development dependency groups
------------------------------------

* ``test``: source-development dependencies, installed with pip dependency-group support.
* ``docs``: source-development dependencies, installed with pip dependency-group support.
* ``release``: source-development dependencies, installed with pip dependency-group support.
* ``quality``: source-development dependencies, installed with pip dependency-group support.
* ``dev``: source-development dependencies, installed with pip dependency-group support.
