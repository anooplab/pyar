Microsolvate
============

``microsolvate`` builds explicit solvent ensembles around one fixed solute.
Its default ``shell`` objective samples targets from the accessible surface of
the original solute. The opt-in ``cluster`` objective also proposes new
molecules around already placed solvent, allowing finite water networks to
grow beyond a uniform solute shell. In both modes the full current cluster is
used for clash checks and backend calculations.

This is different from generic ``grow``:

* ``grow`` places each addend against the complete current structure; adding
  water to existing water can therefore be a legitimate generic growth step.
* ``microsolvate --objective shell`` keeps targeting the original solute
  surface while building the first shell, reducing solvent-island growth.
* ``microsolvate --objective cluster`` mixes surface and solvent-centred
  proposals, then selects by fixed-composition backend energy.

For example:

.. code-block:: console

   pyar microsolvate methane.xyz water.xyz --count 8
   pyar microsolvate methane.xyz water.xyz --count 8 --backend xtb
   pyar microsolvate solute.xyz water.xyz --count 3 --site 7 --backend xtb
   pyar microsolvate solute.xyz water.xyz --count 3 --backend xtb \
       --site-group donor=4 --site-group acceptor=8,9 --contact-distance 3.4

``--site`` accepts one or more 0-based atom indices in the original solute and
restricts targets to accessible surface samples associated with those atoms.
Without it, the entire discretized solute surface is eligible.

``--site`` continues to control initial placement targeting. For separate
post-placement or post-optimization evidence, repeat ``--site-group
NAME=ATOM[,ATOM...]`` and supply an explicit ``--contact-distance`` in
Angstrom. Atom indices refer to the original 0-based solute prefix. Each
reported pair is an atom-distance contact only; PyAR does not infer donor,
acceptor, bond, or metal-coordination chemistry. A solvent molecule may
occupy more than one overlapping group, and each group counts that molecule
once. The report also labels molecules as first shell when they contact the
original solute, second shell when they contact any first-shell solvent, or
outer otherwise, using the same explicit cutoff. These are operational
geometric labels, not thermodynamic shell assignments. Without the opt-in
groups, existing coverage and ``--site`` behavior is unchanged.
With contacts enabled, the report also compares the original-solute prefix
with its starting geometry using a covalent-radii distance heuristic. It flags
a changed proximity graph or reports ``unknown`` if atom identity mapping is
unavailable; this is not a bond-order or reaction analysis.

The shell surface is generated deterministically using Fibonacci points on
atom-centred van der Waals spheres expanded by a spherical probe (1.4 Å by
default). This is a sampled placement surface, not an exact analytical SASA.
Placement checks clashes against the full current cluster. After backend
optimization, a shell-validity post-filter rejects solvent centres that have
moved too far from the targeted surface. PyAR does not claim to apply an
optimization wall or permanent site tether.

The default ``--objective shell --placement-score geometric`` preserves the
coverage-first behavior. Geometry-only mode generates bounded, unoptimized
candidates and does not assign artificial energies. With a backend, candidates
are optimized and shell-invalid structures are rejected before survivor
selection.

The explicit ``--placement-score energy`` option requires a backend with a
registered energy-gradient provider. It computes a full-system single-point
energy for each generated candidate at the current solvent count, retains a
bounded union of low-energy and high-coverage candidates, and optimizes those
survivors. Energies are only compared within one fixed composition and one
backend/electronic-state setup. The prescreen values and survivor decisions are
saved in each step's ``energy_prescreen.json``. Final ranking uses the
post-optimization energy for energy-aware selection.

``--objective cluster --placement-score energy --backend xtb`` explicitly
searches finite clusters rather than requiring every solvent molecule to stay
in the original solute's first shell. The first solvent is placed on the
solute surface; later placements can grow around existing solvent molecules.
After optimization, the shell post-filter is disabled, coverage remains a
reported diagnostic, and candidates are ranked by unbiased full-system energy
with a 0.005 Hartree survivor frontier and structural diversity. Geometry-only
cluster ranking is rejected because it cannot identify low-energy clusters.
This objective is intended for finite-cluster minima, not uniform shell
coverage or bulk hydration. It does not add an optimization wall or a physical
restraint.

Both new options are opt-in and fingerprinted into fresh run requests. Existing
default request serialization remains unchanged for v1/v2 continuation.
Coverage is recorded at each solvent count. A finite gas-phase cluster is not
by itself a bulk-solution model, and no implicit continuum model is inferred
from the solvent XYZ.

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
