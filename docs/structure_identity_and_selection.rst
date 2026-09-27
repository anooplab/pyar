Structure identity and selection
================================

PyAR treats three related questions as separate services:

.. code-block:: text

   Chemical Identity
           ↓
   Geometrical Equivalence
           ↓
   Diversity Selection

Chemical identity asks whether structures describe the same chemical species
or topology. Geometrical equivalence asks whether chemically compatible
structures occupy the same spatial arrangement, conformer, or potential-energy
surface basin. Diversity selection is a resource-allocation decision: given
distinct candidates, which should be retained for expensive calculations.

No single descriptor is required to answer all three questions. In particular,
clustering is not a definition of chemical identity. The comparison interfaces
in ``pyar.structure_comparison`` allow identity providers and geometry
comparators to evolve independently of diversity selection and workflow
orchestration.

This architectural foundation preserves current algorithms and defaults. The
OpenBabel/InChI product identity decision, Coulomb-fingerprint prefilter,
permutation-aware Kabsch RMSD, adaptive duplicate threshold, and current
clustering/selection behavior remain in use. Alternative geometrical metrics
and selection algorithms require separate scientific evaluation.
