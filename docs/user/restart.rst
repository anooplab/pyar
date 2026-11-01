Restart interrupted runs safely
===============================

Use the same supported PyAR series, scientific request, input ordering and output
directory. The workflow checks hashes, immutable commits and request identity
before reusing accepted evidence. A derived summary or copied final XYZ alone
does not authorize restart.

.. code-block:: bash

   pyar inspect grow
   pyar grow seed.xyz ligand.xyz --count 2 --backend xtb --orientations 4 --check

Repeat the original command only after checking its supported restart behavior.
Some tasks intentionally require a fresh directory; conformer search and
standalone fixed-geometry batches have different lifecycle rules from growth or
reaction discovery. Check the task's diagnostic rather than forcing a generic
resume policy onto it.

A stale-lock recovery option is explicit permission to recover a demonstrably
stale writer; it must not break an active process. Preserve interrupted-attempt
records and RNG/proposal state. Replacing files, deleting locks blindly or editing
status fields can corrupt provenance or change the scientific ensemble.

Older run generations are unsupported in 1.6. Use their original version and
matching dependencies for reproduction. Starting a new 1.6 calculation from an
old XYZ is a new request and must explicitly specify the desired electronic state
and settings. See :doc:`/reference/upgrade`.
