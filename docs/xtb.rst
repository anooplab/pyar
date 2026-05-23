xTB Guide
=========

This page summarizes how PyAR uses xTB at runtime.

Wrappers
--------

PyAR uses xTB in four places:

* ``pyar.interface.xtb`` for standalone xTB optimization
* ``pyar.interface.xtb_aiqm1`` for xTB followed by AIQM1 refinement
* ``pyar.interface.xtb_aimnet2`` for xTB followed by AIMNet2 refinement
* ``pyar.interface.xtbturbo`` for Turbomole/AFIR workflows that call xTB for gradients

Shared behavior
---------------

* ``--nprocs`` is passed through to xTB as ``--parallel``.
* Charge, multiplicity, and spin handling are normalized in
  ``pyar.interface.xtb_utils``.
* The wrappers fail fast if ``xtb`` is not available on ``PATH``.

See the full API documentation in :doc:`api`.
