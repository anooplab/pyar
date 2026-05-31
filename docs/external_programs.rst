External Program Requirements
=============================

PyAR can be installed with Python, but many chemistry backends are separate
programs or libraries that must be installed independently.

Official project sites
----------------------

.. list-table::
   :header-rows: 1
   :widths: 18 32 50

   * - Program
     - Official site
     - Notes
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
   * - MOPAC
     - https://openmopac.net/download/installer/
     - Use the official installer or the upstream installation manual.
   * - MLatom
     - https://mlatom.com/docs/installation.html
     - Python package and workflow ecosystem with upstream install guidance.
   * - Open Babel
     - https://openbabel.org/docs/Installation/install.html
     - Command-line tools and optional Python bindings are installed
       separately from PyAR.
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
