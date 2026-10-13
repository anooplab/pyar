External program requirements
=============================

PyAR installs Python requirements through workflow and capability extras.
It discovers external computational programs on PATH before calculation.
None of the programs below is installed by PyAR's pip extras.

.. list-table:: Installation and discovery
   :header-rows: 1
   :widths: 15 25 60

   * - Program
     - Executable / verification
     - Official installation source
   * - xTB
     - ``xtb --version``
     - https://xtb-docs.readthedocs.io/en/latest/setup.html;
       ``mamba install -c conda-forge xtb`` is an available installation route.
   * - ORCA
     - ``command -v orca``; small calculation per vendor tutorial
     - https://www.faccts.de/docs/orca/6.1/tutorials/first_steps/install.html
   * - Gaussian 16
     - ``command -v g16``; vendor sample calculation
     - https://gaussian.com/g16/g16src_install.pdf (licensed distribution).
   * - Turbomole
     - ``command -v define ridft rdgrad statpt jobex``
     - https://www.turbomole.org/turbomole/turbomole-documentation/
   * - MOPAC
     - ``command -v mopac``; vendor sample calculation
     - https://openmopac.net/
   * - Psi4
     - ``psi4 --version``
     - https://psi4.github.io/psi4docs/master/conda.html
   * - OpenBabel
     - ``obabel -V``; ``command -v obminimize obenergy``
     - https://openbabel.org/docs/Installation/install.html;
       ``mamba install -c conda-forge openbabel``.

On Windows use ``where PROGRAM`` in place of ``command -v PROGRAM``.
OpenBabel is used for XYZ-to-chemical-identity conversion, force-field
optimization and MOPAC input generation. AIMNet2's explicit-charge XYZ route
requires only its Python extra; other file formats or charge perception in
the standalone helper may require OpenBabel Python bindings.

MLatom/AIQM1 routes may also require Gaussian and model resources according
to https://mlatom.com/docs/installation.html. DFT-D4 is a backend method
setting where supported; PyAR currently has no standalone DFT-D4 adapter.

``pyar doctor`` reports external discovery without launching programs.
``--check`` validates the selected route, including Python requirements and
external model assets. Use the vendor's small calculation to validate a
computational program beyond executable discovery.

See :doc:`installation` for workflow installation examples.
