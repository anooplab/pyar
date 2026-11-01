Build explicit solvent shells
=============================

Microsolvation places solvent molecules relative to the original solute surface.
It has a shell/coverage objective that differs from unconstrained cluster-energy
minimization. Choose the objective and confinement model for your question.

.. code-block:: bash

   pyar microsolvate solute.xyz water.xyz --count 2 --orientations 4 --check
   pyar microsolvate solute.xyz water.xyz --count 2 --orientations 4
   pyar microsolvate solute.xyz water.xyz --count 2 --orientations 4 --backend xtb --check

Use an appropriate solute XYZ and the water example. Geometry-only mode checks
placement/shell behavior; add the selected backend to optimize. Inspect coverage
history, solute atom indices, solvent fragments, objective and selected outputs
under ``microsolvation/results``. Inspect actual structures for bridging,
coordination and unwanted dissociation rather than relying on coverage alone.

Search and resampling budgets, confinement and solvent count are scientific
settings. Keep them with the ensemble. Current restarts preserve committed
proposal/RNG state; do not manually reseed or remove a pending attempt to make a
run appear completed. A shell ensemble is not an implicit-solvent free energy
or a bulk-solvation calculation.


See :doc:`/reference/cli` for current options and :doc:`results` for interpretation.
