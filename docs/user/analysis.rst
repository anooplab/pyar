Analyze existing results
========================

Choose the right evidence
-------------------------

Work from a copy or an explicitly selected result manifest; analysis should not
modify the calculation's immutable records. First confirm matching composition,
charge/spin and energy model. Geometry similarity and chemical identity answer
different questions: a small RMSD is not necessarily the same connectivity, and
a matching graph can still contain several conformations.

.. code-block:: bash

   pyar energies selected_a.xyz selected_b.xyz
   pyar compare selected_a.xyz selected_b.xyz
   pyar identify selected_a.xyz
   pyar clustering selected_a.xyz selected_b.xyz --help
   pyar trace reaction --help

Replace the example filenames with selected XYZs from your run. ``energies``
reads stored XYZ energy metadata and reports relative energies in kcal/mol.
It does not calculate missing energies or establish cross-method comparability.
``compare`` reports structural/energetic differences under its chosen semantics;
``identify`` reports perceived identity with documented limitations.

Selection and diversity
------------------------

Use ``deduplicate`` to remove structural duplicates and ``select`` or ``clustering``
to construct a representative subset. Inspect the feature/comparator, thresholds,
ordering and diagnostics; changing them changes the selected ensemble. Preserve
source names/paths so every selected structure can be traced to its calculation.
See :doc:`/reference/tools` for the specialist pool similarity and descriptor tools.

Interpret a search
-------------------

Count attempted, successful, rejected and unqualified candidates. Compare
representatives across pathways, orientations and seeds. A larger sample can
change a ranking without a software defect. Inspect suspicious close contacts,
fragment rearrangements and spin/model choices. Report the sampling protocol and
selection policy with scientific conclusions, not just the lowest number.
