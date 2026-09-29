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

Product identity remains separate from geometry comparison. The
``CoulombEigenvalueRMSDComparator`` applies the sorted Coulomb eigenvalue
prefilter followed by permutation-aware Kabsch RMSD. Candidate-pool
deduplication now uses an element-labeled inferred graph followed by
permutation-aware Kabsch RMSD. It removes a candidate only when graph comparison
completes and RMSD is below the adaptive threshold. If graph connectivity
matches but graph mapping reaches its cap, PyAR runs iRMSD in both argument
orders in an isolated process. It removes the candidate only if both calls
finish without diagnostics and both distances are below the threshold; it uses
the larger distance. Graph mismatches, backend warnings, errors, and remaining
incomplete comparisons retain both candidates. This follows “in doubt, keep.”

Optional RMSD strategies
------------------------

``IRMSDComparator`` is an optional element-permutation-invariant geometry
metric. Install it with ``pip install 'pyar-chem[structure-comparison]'``.
It does not establish molecular identity or check connectivity; equal-formula
constitutional isomers can receive a finite distance. Keep identity checks
separate when interpreting its result.

``GraphRMSDComparator`` infers an element-labeled graph from the XYZ geometry
using covalent radii and a configurable ``bond_scale`` (default ``1.15``).
It computes a proper-rotation Kabsch RMSD only when those inferred graphs are
isomorphic. This prevents an RMSD-only match across different inferred
connectivities, but inferred bonds remain a geometric heuristic and do not
encode bond orders. It is the default comparator used by deduplication.
``IRMSDComparator`` is also used as a gated secondary check by the deduplication
policy. A graph match is required first, and native output is captured so an
internal topology fallback cannot silently authorize deletion.
