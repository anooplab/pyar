Backend capabilities and installation guides
============================================

A backend supplies energies/gradients or native optimization; an optimizer uses
those values. Capability is workflow-specific. An installed package does not
prove an executable/model works for every calculation.

.. code-block:: bash

   pyar backends
   pyar install status
   pyar doctor

Use ``--check`` with the actual selected command to validate its route. The
standalone energy/gradient/frequency commands and the geometry-optimization
workflows have different supported routes. Do not infer a universal backend
matrix from an import succeeding.

.. toctree::
   :maxdepth: 1

   backends/xtb
   backends/orca_backend
   backends/aimnet2_backend
   backends/gaussian_backend
   backends/mopac_backend
   backends/obabel_backend
   backends/psi4_backend

For Turbomole and optional ML programs, consult the installed capability and
missing-dependency diagnostic before running. Keep proprietary executable/model
installation separate from PyAR packaging. Use the vendor's current license,
platform and execution requirements. Optional backends must remain optional for
root help and tasks that do not use them.
