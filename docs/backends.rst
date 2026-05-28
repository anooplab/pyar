Backend Guides
==============

PyAR can use different computational backends for different chemistry tasks.
A backend may be a semiempirical program, a quantum-chemistry program, a
machine-learning potential, or a hybrid route.

For most users, the practical question is simple: choose a backend that is
installed, appropriate for your chemistry, and supported by the PyAR task you
want to run.

Current backend families
------------------------

.. list-table::
   :header-rows: 1
   :widths: 20 24 56

   * - Backend
     - Family
     - Typical role in PyAR
   * - ``xtb``
     - semiempirical
     - Fast optimisation, aggregation, solvation, and AFIR/geomeTRIC reaction
       searches.
   * - ``aimnet_2``
     - machine-learning potential
     - Fast molecular energy/force evaluation for optimisation and
       AFIR/geomeTRIC reaction searches.
   * - ``orca``
     - quantum chemistry
     - Higher-level DFT-style backend. The geomeTRIC route is wired, but local
       executable setup and validation are required.
   * - ``gaussian``
     - quantum chemistry
     - Higher-level DFT-style backend. The geomeTRIC route is wired, but local
       executable setup and validation are required.
   * - ``mopac``
     - semiempirical
     - Legacy semiempirical optimisation route.
   * - ``obabel``
     - molecular mechanics / OpenBabel route
     - Lightweight structure preparation or optimisation where charge and
       multiplicity handling are not central.
   * - ``psi4``
     - quantum chemistry
     - Available as a backend family, but not currently on the registered
       AFIR/geomeTRIC energy-gradient route.
   * - ``turbomole`` / ``xtb_turbo``
     - quantum chemistry / compatibility route
     - Legacy or compatibility workflows. New AFIR/geomeTRIC development should
       not depend on Turbomole coordinate updates.
   * - ``xtb-aimnet2`` / ``xtb-aiqm1``
     - hybrid
     - Staged routes that combine fast pre-optimisation with ML or AIQM-style
       refinement where available.

Which backend should I start with?
----------------------------------

For a new user:

* Start with ``xtb`` when you want a robust and fast first calculation.
* Try ``aimnet_2`` when you want a machine-learning-potential route and your
  molecule is within the chemistry where the model is expected to behave well.
* Use ``orca`` or ``gaussian`` only when the executable is installed and you
  have validated a small test case on your machine.
* Use legacy or hybrid routes only when you know why you need them.

AFIR/geomeTRIC support
----------------------

The current registered energy-gradient providers for geomeTRIC-backed AFIR
reaction optimisation are:

.. code-block:: text

   xtb
   aimnet_2
   orca
   gaussian

This means PyAR can, in principle, build the objective
``backend energy + AFIR bias`` and give it to geomeTRIC/TRIC. In practice,
backend installation, executable paths, charge/multiplicity support, and local
validation still matter.

Backend pages
-------------

.. toctree::
   :maxdepth: 1

   xtb
   aimnet2_backend
   orca_backend
   gaussian_backend
   mopac_backend
   obabel_backend
   psi4_backend
