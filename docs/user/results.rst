Read and interpret results
==========================

Find the scientific result
---------------------------

Read the task summary first, then the per-job result and raw backend evidence.
``pyar inspect RUN_DIRECTORY`` helps locate recorded status. Search workflows
publish selected geometries under their documented results directory; execution
scratch files are not interchangeable with selected, accepted structures.

.. list-table::
   :header-rows: 1

   * - Quantity/evidence
     - Interpretation
   * - Electronic energy (hartree)
     - Compare only matching composition, model, electronic state and conventions
   * - Coordinates (angstrom)
     - Preserve atom order, fragments and constraints when comparing geometries
   * - Gradient (hartree/bohr in standalone results)
     - Inspect both component units and the workflow's convergence criterion
   * - Native optimizer convergence
     - A numerical geometry criterion; it does not independently prove a minimum
   * - ``minimum_confirmed``
     - The supplied geometry passed the implemented minimum qualification
   * - First-order saddle gate
     - Requires the implemented frequency/Hessian acceptance, not just TS optimization
   * - Reaction connection gate
     - Requires both IRC directions, endpoint minima and endpoint identity checks

Compare ensembles
-------------------

.. code-block:: bash

   pyar energies results/*.xyz
   pyar compare structure_a.xyz structure_b.xyz
   pyar identify structure_a.xyz

Use paths that actually exist in your workflow; see :doc:`analysis` for selecting
them. An energy ranking is not a thermal population. Differences in atom count,
charge, method, basis, solvent or bias objective can invalidate a naive table.
A biased AFIR objective must not be reported as an unbiased product energy.
Electronic energy differences alone omit zero-point, thermal, entropic and
solvation contributions unless explicitly calculated and consistently combined.

Keep uncertainty visible
--------------------------

Retain failed and unqualified jobs when interpreting a search. A small successful
subset can hide missing basins, failed electronic states or disconnected
pathways. Check counts and budgets in the run record. Use structural diversity,
chemical identity and relevant vibrational/path evidence alongside energies.
Do not overwrite original calculations to make a report easier to read.
