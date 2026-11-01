Choose a workflow
=================

.. list-table::
   :header-rows: 1

   * - Your question
     - Start with
     - What to inspect
   * - What is the energy at this geometry?
     - ``pyar energy`` / ``pyar engrad``
     - Energy, gradients, electronic state and units
   * - Is this a local minimum?
     - ``pyar optimize`` followed by ``pyar frequency``
     - Accepted geometry, stationarity and vibrational qualification
   * - Which conformations are accessible to sampling?
     - ``pyar conformer``
     - Selected conformers, diversity and refinement method
   * - How do these fragments assemble?
     - ``pyar aggregate``
     - Composition, selected ensemble and comparable energies
   * - What changes when I add another ligand?
     - ``pyar grow``
     - Stepwise seed composition and selected pools
   * - Where can explicit solvent molecules sit?
     - ``pyar microsolvate``
     - Shell coverage, solvent count and optimization objective
   * - Which products might be discovered under a bias?
     - ``pyar react``
     - Bias/release evidence and product identity
   * - How does a chosen bond coordinate evolve?
     - ``pyar scan-bond``
     - Constrained profile and freely relaxed products
   * - Do endpoints connect through a validated saddle?
     - ``pyar neb``, ``pyar ts``, ``pyar irc``
     - Each stage's gates, Hessian and both endpoint assignments

Use :doc:`input_and_settings` to choose charge, spin and backend before sampling.
A search samples candidates; absence from its results does not establish that
a structure or pathway cannot exist. Increase sampling only with a clear budget
and compare like-for-like chemistry.
