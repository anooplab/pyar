Scan a bond distance
====================

A relaxed bond scan constrains one chosen coordinate while relaxing the others.
It probes that coordinate; it does not guarantee a minimum-energy reaction path.
The scan evidence and the freely relaxed product are different artifacts.

.. code-block:: bash

   pyar scan-bond water.xyz --atoms 0 1 --scan-end 1.2 --scan-points 3 --backend xtb --check
   pyar scan-bond water.xyz --atoms 0 1 --scan-end 1.2 --scan-points 3 --backend xtb

For one molecule, indices are zero-based within that molecule. For two fragments,
the two indices are local to fragment A and fragment B. Check atom ordering and
the scan direction in the resolved request before running. The water example
stretches an O–H coordinate; it is an illustrative probe, not a claim about a
qualified reaction mechanism.

Inspect ``scan_bond/summary.json`` and per-orientation profile, trajectory and
result evidence. Keep original constrained artifacts unchanged. Product grouping
uses the freely relaxed candidate and the workflow's identity policy. Energy
maxima can seed further investigation, but a sampled maximum is not itself a
validated transition state. ``--through`` and unique-product characterization
can continue through the staged pathway checks; inspect every required gate.


See :doc:`/reference/cli` for current options and :doc:`results` for interpretation.
