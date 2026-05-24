xTB Guide
=========

This page summarizes how PyAR uses xTB at runtime.

Wrappers
--------

PyAR uses xTB in four places:

* ``pyar.interface.xtb`` for standalone xTB optimization
* ``pyar.interface.xtb_aiqm1`` for xTB followed by AIQM1 refinement
* ``pyar.interface.xtb_aimnet2`` for xTB followed by AIMNet2 refinement
* ``pyar.interface.xtb_turbo`` for the current Turbomole/AFIR workflow that
  calls xTB for gradients

The ``xtb_turbo`` path is retained only as current behaviour during the
redesign. It requires Turbomole for coordinate updates and convergence
checking even though xTB supplies the physical gradient. The target reaction
optimization design replaces it with an internal-coordinate optimizer over a
composed xTB plus AFIR energy/gradient objective. See
:doc:`reaction_optimization`.

Shared behavior
---------------

* ``--nprocs`` is passed through to xTB as ``--parallel``.
* Charge, multiplicity, and spin handling are normalized in
  ``pyar.interface.xtb_utils``.
* The wrappers fail fast if ``xtb`` is not available on ``PATH``.

See the full API documentation in :doc:`api`.
