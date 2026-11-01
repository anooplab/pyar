Energy, gradients, optimization and frequencies
===============================================

Use fixed-geometry commands to measure the supplied structure. ``energy`` writes
an electronic energy; ``engrad`` adds Cartesian gradients; ``frequency`` uses
gradient finite differences and qualifies a minimum. The frequency workflow
never performs a hidden optimization.

.. code-block:: bash

   pyar energy water.xyz --backend xtb --output energy_water --check
   pyar engrad water.xyz --backend xtb --output gradient_water
   pyar optimize water.xyz --backend xtb
   pyar frequency job_water/result_water.xyz --backend xtb --output frequency_water

Download :download:`water.xyz </examples/water.xyz>`. Standalone output roots
contain ``request.json``, ``summary.json``/``summary.csv``, and isolated per-input
jobs. For the optimized input the name is ``result_water``, so its frequency job
is ``job_result_water``. Read ``result.json`` for units, qualification and evidence
paths. Finite-difference frequencies scale with the number of Cartesian degrees
of freedom and require many backend evaluations; use the printed budget before
applying this to a large molecule.

Charge/spin, solvent and method settings apply to every input in a batch. Use
separate commands for different electronic states. Request-level validation
checks all inputs before calculations begin. Keep the backend evidence and
avoid comparing absolute energies from different methods or compositions.


See :doc:`/reference/cli` for current options and :doc:`results` for interpretation.
