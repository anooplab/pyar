External Program Requirements
=============================

PyAR can be installed with Python, but many chemistry backends are separate
programs or libraries that must be installed independently.

The fastest path is usually:

* ``python -m pip install "pyar-chem[selection]"`` for clustering and
  descriptor helpers such as ``ase``, ``dscribe``, ``hdbscan``, ``pandas``,
  ``scikit-learn``, and ``MDAnalysis``
* ``python -m pip install "pyar-chem[ml]"`` for ML helpers such as
  ``torchani``, ``mlatom``, ``pyh5md``, and ``h5py``
* ``python -m pip install "pyar-chem[aimnet2]"`` when you need the AIMNet2
  runtime stack
* ``python -m pip install "pyar-chem[xtb]"`` for the geomeTRIC channel used by
  xTB-backed AFIR work
* ``python -m pip install "pyar-chem[openbabel]"`` for the Open Babel Python
  bindings

If you only want a single Python package, install it directly with
``python -m pip install <package-name>``. For example:

* ``python -m pip install ase``
* ``python -m pip install geometric``
* ``python -m pip install dscribe``
* ``python -m pip install hdbscan``
* ``python -m pip install pandas``
* ``python -m pip install scikit-learn``

Use the upstream site or vendor installer for tools that are not distributed
that way, especially executables.

Official project sites
----------------------

.. list-table::
   :header-rows: 1
   :widths: 18 32 50

   * - Program
     - Official site
     - Notes
   * - ASE
     - https://wiki.fysik.dtu.dk/ase/install.html
     - Install with ``python -m pip install ase``.
   * - geomeTRIC
     - https://github.com/leeping/geomeTRIC
     - Install with ``python -m pip install geometric``.
   * - xTB
     - https://xtb-docs.readthedocs.io/en/latest/setup.html
     - Use the upstream setup guide for binaries, environment variables, and
       local installation details.
   * - ORCA
     - https://orca-manual.mpi-muelheim.mpg.de/contents/quickstartguide/installation.html
     - Download from the ORCA forum or the official installer package.
   * - Gaussian
     - https://gaussian.com/
     - Commercial software; install through the official vendor channels.
   * - Psi4
     - https://psicode.org/psi4manual/master/index.html
     - Open-source quantum chemistry package with official installation docs.
   * - DScribe
     - https://singroup.github.io/dscribe/latest/
     - Install with ``python -m pip install dscribe`` or
       ``python -m pip install "pyar-chem[selection]"``.
   * - MOPAC
     - https://openmopac.net/download/installer/
     - Use the official installer or the upstream installation manual.
   * - MLatom
     - https://mlatom.com/docs/installation.html
     - Install with ``python -m pip install mlatom`` or
       ``python -m pip install "pyar-chem[ml]"``.
   * - Open Babel
     - https://openbabel.org/docs/Installation/install.html
     - Install the Python bindings with ``python -m pip install openbabel`` or
       ``python -m pip install "pyar-chem[openbabel]"``; the command-line
       tools still come from the upstream installer.
   * - TURBOMOLE
     - https://www.turbomole.org/
     - Commercial software; use the official documentation and installer.
   * - DFT-D4
     - https://github.com/dftd4/dftd4
     - Separate dispersion-correction executable used by some workflows.

What PyAR installs
------------------

PyAR installs the Python package and small support files. It does not install
the external executables above, and the ``pyar-chem`` wheel does not bundle
large AIMNet2 or MLatom model/vendor assets. Optional extras install Python
dependencies only.

See also
--------

* :doc:`installation`
* :doc:`backends`
