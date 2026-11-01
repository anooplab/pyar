Add ligands or fragments to a fixed seed
========================================

Growth repeatedly adds one species while retaining the seed's role. Use this
for a stepwise ligand/fragment question; use ``aggregate`` for an unprivileged
final composition and ``microsolvate`` for original-solute-centred shell sampling.

.. code-block:: bash

   pyar grow seed.xyz ligand.xyz --count 2 --backend xtb --orientations 4 --check
   pyar grow seed.xyz ligand.xyz --count 2 --backend xtb --orientations 4

Provide your own seed/ligand XYZs with valid states. Without a backend, supported
growth routes produce geometry-only candidates; those do not have electronic
energies or optimization evidence. Do not present them as optimized structures.

Inspect each step's selected pool and ``grow/results``. Confirm atom counts,
fragment identities and intended connectivity. An early selection decision
changes later candidates, so preserve per-step budgets and selection settings.
Restart uses supported committed records and verified artifact references;
editing a selected XYZ invalidates its authority. See :doc:`restart`.


See :doc:`/reference/cli` for current options and :doc:`results` for interpretation.
