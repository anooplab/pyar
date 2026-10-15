Microsolvate
============

``microsolvate`` builds an explicit solvent-shell ensemble around one fixed
solute. During construction, placement targets are sampled from the accessible
surface of the original solute. Previously added solvent molecules block
occupied regions and participate in clash checks and backend energy
calculations, but they do not become new first-shell target surfaces.

This is different from generic ``grow``:

* ``grow`` places each addend against the complete current structure; adding
  water to existing water can therefore be a legitimate generic growth step.
* ``microsolvate`` keeps targeting the original solute surface while building
  the first shell, reducing one-sided solvent-island growth.

For example:

.. code-block:: console

   pyar microsolvate methane.xyz water.xyz --count 8
   pyar microsolvate methane.xyz water.xyz --count 8 --backend xtb
   pyar microsolvate solute.xyz water.xyz --count 3 --site 7 --backend xtb

``--site`` accepts one or more 0-based atom indices in the original solute and
restricts targets to accessible surface samples associated with those atoms.
Without it, the entire discretized solute surface is eligible.

The surface is generated deterministically using Fibonacci points on
atom-centred van der Waals spheres expanded by a spherical probe (1.4 Å by
default). This is a sampled placement surface, not an exact analytical SASA.
Placement checks clashes against the full current cluster. After backend
optimization, a shell-validity post-filter rejects solvent centres that have
moved too far from the targeted surface. PyAR does not claim to apply an
optimization wall or permanent site tether.

Geometry-only mode generates and selects bounded, unoptimized candidates. It
does not assign artificial energies. Backend mode uses energies only after
shell-validity filtering; coverage and structural diversity remain part of
survivor selection. Coverage is recorded at each solvent count. A finite
gas-phase cluster is not by itself a bulk-solution model, and no implicit
continuum model is inferred from the solvent XYZ.

This workflow is motivated by solute-centred accessible-surface and
first-shell sampling concepts used in methods such as systematic microsolvation,
CREST QCG, ORCA SOLVATOR, and Fibonacci-surface approaches. PyAR does not claim
to reproduce any of these methods exactly.

The legacy ``pyar-cli solvate`` interface remains available through its
compatibility algorithm and existing ``solvation/state.json`` restart format.
It emits a deprecation warning. Existing legacy state is not interpreted as a
new microsolvation run.

Output layout and continuation
------------------------------

Fresh runs retain the public ``microsolvation/`` root, with explicit v2
``layout.json`` metadata, immutable ``inputs/{solute,solvent}.xyz``, and the
original ``surface/initial_surface.json``. Each ``steps/step_NNN/`` contains
selected snapshots, coverage and an immutable ``selection.json`` commit. Native
candidates and calculator job/cache files stay together under
``attempts/attempt_000/execution/candidates/``. Final copies and their manifest
are in ``results/``; the root summary and manifest are rebuildable reports.

Committed records retain the original-solute prefix, solvent-fragment identity,
shell coverage, scientific request and upstream byte hashes. Restart validates
all these records and adopts a committed step whose progress update was
interrupted without repeating that step. Completed result bytes and stable
artifact IDs are retained on rerun. An altered original surface or committed
survivor blocks reuse. The existing shell post-filter, site restrictions and
coverage/energy/diversity selection rules are unchanged.

Existing v1 microsolvation runs continue through their original ``step_NNN/``
and ``final/`` paths and are never converted. Legacy ``pyar-cli solvate`` stays
on its own v1 solvation implementation. Modern ``pyar solvate`` remains an
alias for this microsolvation command. New v2 execution has a local POSIX
single-writer guard; after a process crash use ``--recover-stale-lock``
explicitly. It cannot break an active writer. ``--check`` and ``--dry-run``
only validate and do not acquire a writer lock or mutate calculation output.
