Aggregate
=========

Modern interface
----------------

.. code-block:: bash

   pyar aggregate A.xyz B.xyz
   pyar aggregate water.xyz --size 6 --backend xtb
   pyar aggregate water.xyz ammonia.xyz --size 4 2 --backend xtb
   pyar aggregate --formula CH4 --check

Multiple inputs default to one of each. A single input requires ``--size``;
explicit sizes contain one positive count per input. Inputs describe building
block types for a final composition, with existing alternative build pathways;
no component is a permanently privileged seed.

Backend refinement is optional. Geometry-only search uses the existing bounded
structural selection without artificial energies. Orientations default to eight,
and the survivor budget remains eight. Selection/connectivity/pathway overrides
retain their established policies. Formula input specifies composition, not bonds.

``--check`` validates inputs, electronic states, backend requirements when
selected, and the existing ``aggregates/state.json`` restart request without
creating output or modifying state. Legacy aggregate syntax remains supported.

Use aggregation when you want PyAR to build and screen low-energy structures
from fragments or a formula. This is the main workflow for clusters,
noncovalent complexes, and other build-up problems.

Disconnected covalent graphs are expected for molecular aggregates and other
noncovalent complexes. Use ``--connectivity-policy auto`` for the default
chemistry-aware choice, or override with ``off``, ``prefer``, or ``strict``
when you need a different selection rule.

Basic commands
--------------

.. code-block:: bash

   pyar-cli aggregate C H -as 1 4 -N 8
   pyar-cli -a C H -as 1 4 -N 8
   pyar-cli --aggregate --formula C5H4 -N 8

What it does
------------

* generates trial geometries
* evaluates and ranks candidates
* removes near-duplicates
* persists restart state in ``aggregates/state.json``

Useful outputs
--------------

* ``aggregates/state.json`` for restart and provenance
* ``selected/`` for the chosen structures
* the energy table for quick inspection of relative energies

See also
--------

* :doc:`quickstart`
* :doc:`installation`
* :doc:`workflows`

Aggregate versus grow
---------------------

Aggregate searches the requested final composition without a permanently
privileged seed and may explore alternative build pathways. For a specified
seed with one repeatedly added species and retained intermediate stages, use
:doc:`grow`. Both use the same single-addition service.

Output layout for fresh runs
----------------------------

Fresh runs use an explicit v2 marker under the unchanged ``aggregates/`` root.
Original components are snapshotted in ``inputs/``. Stable pathway metadata in
``pathways/pathway_NNN/pathway.json`` records component identities and the exact
addition sequence. Each ``steps/step_NNN/`` retains a committed selected pool;
its native ``ag_*`` calculator/cache tree stays inside the step's execution
attempt. Verified additions can be reused after interruption. Existing v1 runs
and established ``pyar.log`` imports retain their original implementation and
paths, and are never converted. Completed aggregate runs remain non-resumable.

Final cross-path selection consumes explicit committed source references.
Legacy logical names and ordering are retained for tie-breaking and duplicate
selection. Native final-clustering outputs remain internal under
``selection/attempts/attempt_000/execution/``; flat final copies and a manifest
are in ``results/``. ``final_selection.json`` records selection provenance.
The root summary and manifest are derived reports.

The existing selector groups surviving structures by composition, including
smaller compositions retained during growth. This behavior is preserved, not
silently filtered to the requested composition. Manifest entries state their
actual composition and whether it matches the requested target. These are final
cross-path selections for each composition, distinct from unselected candidates.

New v2 execution uses a local POSIX single-writer guard. The modern CLI accepts
``--recover-stale-lock`` after a process crash; it cannot break an active lock.
Read-only preflight validates existing records without locking or mutating the
calculation output.
