Installation
============

This page explains two different things:

* installing the Python package ``pyar``
* installing optional external programs such as xTB, ORCA, Gaussian, Psi4,
  MOPAC, MLatom, and OpenBabel

PyAR can be installed with Python, but several computational backends are
separate programs. They do not automatically come with PyAR.

Install PyAR
------------

From a local checkout:

.. code-block:: bash

   git clone https://github.com/anooplab/pyar.git
   cd pyar
   python -m pip install .

For development:

.. code-block:: bash

   python -m pip install -e .

Check that the command-line entry point is available:

.. code-block:: bash

   pyar-cli --help

What PyAR installs automatically
--------------------------------

A normal pip-based PyAR install installs the Python package metadata and Python
dependencies declared by PyAR. Python dependencies include libraries used by
PyAR itself, such as clustering, data handling, and optimisation helpers.

``geomeTRIC`` is a Python dependency rather than an external electronic-
structure program. It is used as the internal-coordinate optimiser for AFIR
reaction calculations with registered energy/force providers such as ``xtb``,
``aimnet_2``, ``orca``, and ``gaussian``.

OpenBabel has two parts:

* the Python binding may be installed as a Python dependency
* command-line executables such as ``obabel``, ``obminimize``, and
  ``obenergy`` must still be available on your system if you use the OpenBabel
  backend route

External programs not installed by PyAR
---------------------------------------

The following programs are separate from PyAR. Install only the ones you need.

.. list-table::
   :header-rows: 1
   :widths: 18 28 54

   * - Program
     - Licence/access
     - Where to get installation instructions
   * - xTB
     - Free/open source
     - Homepage: https://github.com/grimme-lab/xtb ; documentation and setup:
       https://xtb-docs.readthedocs.io/en/latest/setup.html
   * - ORCA
     - Free for academic use after registration; see ORCA terms
     - Homepage/downloads: https://www.faccts.de/orca/
   * - Gaussian
     - Commercial software
     - Product site: https://gaussian.com/ ; obtain and install through your
       institution, licence administrator, or Gaussian vendor instructions
   * - Psi4
     - Free/open source
     - Homepage: https://psicode.org/ ; install page:
       https://psicode.org/installs/latest/
   * - MOPAC
     - Check current OpenMOPAC/MOPAC distribution terms
     - Homepage: https://github.com/openmopac/mopac
   * - MLatom
     - Free academic software; check upstream terms
     - Homepage: https://mlatom.com/ ; installation manual:
       https://mlatom.com/manual/installation.html
   * - OpenBabel
     - Free/open source
     - Homepage: https://openbabel.org/ ; installation instructions:
       https://openbabel.org/docs/Installation/install.html
   * - Turbomole
     - Commercial software
     - Homepage: https://www.turbomole.org/

Backend-specific notes
----------------------

xTB
~~~

xTB is the easiest backend for most new users. It is fast and is the best first
backend to test PyAR.

After installing xTB, check:

.. code-block:: bash

   xtb --version

See :doc:`xtb` for PyAR-specific xTB notes.

ORCA
~~~~

ORCA is not installed by PyAR. Download it from the official ORCA/FACCTs site
after registration, then follow the platform-specific installation instructions
provided there.

Check that your shell can find ORCA before using it from PyAR:

.. code-block:: bash

   orca --version

ORCA support in the geomeTRIC/AFIR route is wired through an energy-gradient
provider, but you should validate a small ORCA test case on your own machine
before relying on it for production calculations.

Gaussian
~~~~~~~~

Gaussian is commercial software. PyAR does not install it. Your institution or
licence administrator must provide access and installation instructions.

Check your local Gaussian command name. On many systems it is ``g16``:

.. code-block:: bash

   g16 < input.com > output.log

Gaussian support in the geomeTRIC/AFIR route is wired through an energy-
gradient provider, but local validation is required because Gaussian setup and
output handling are installation-specific.

Psi4
~~~~

Psi4 is free and open-source. Use the official Psi4 install page for the
current recommended conda or package-manager command:

https://psicode.org/installs/latest/

Check:

.. code-block:: bash

   psi4 --version

At present, ``psi4`` is not registered as an energy-gradient provider for the
geomeTRIC-backed AFIR reaction route.

MLatom
~~~~~~

MLatom is a separate Python package and program ecosystem used by some PyAR
routes. In many environments it can be installed with pip:

.. code-block:: bash

   python -m pip install mlatom

Use the upstream installation manual for current instructions:

https://mlatom.com/manual/installation.html

OpenBabel
~~~~~~~~~

OpenBabel command-line tools must be installed separately if you use the
OpenBabel route. Check:

.. code-block:: bash

   obabel -V
   obminimize --help

Use the upstream installation instructions:

https://openbabel.org/docs/Installation/install.html

Turbomole
~~~~~~~~~

Turbomole is commercial software and is used only for legacy or compatibility
routes in PyAR. New AFIR/geomeTRIC development should not depend on Turbomole
coordinate updates.

Troubleshooting
---------------

If PyAR says an executable is missing:

1. Check whether the program is installed.
2. Check whether the executable is on your ``PATH``.
3. Run the executable directly, for example ``xtb --version`` or
   ``psi4 --version``.
4. Then rerun the PyAR command.

A good first PyAR check is:

.. code-block:: bash

   pyar-cli --help

Then try the example in :doc:`first_run`.
