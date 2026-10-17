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

An additional experimental, opt-in outer-sphere wall is available only for
the ``shell`` objective with an energy-gradient backend and geomeTRIC:

.. code-block:: console

   pyar microsolvate solute.xyz water.xyz --count 3 --backend xtb \
       --geometry-optimizer geometric --wall outer-sphere \
       --wall-radius 8.0 --wall-force-constant 0.5

The flat-bottom harmonic penalty acts on solvent-fragment centers of mass
outside a sphere centered on the instantaneous solute center of mass. Its
force translates each solvent fragment without distorting its internal
geometry. The radius is in
Angstrom and the force constant is in hartree/Angstrom squared. Both values
are required. The wall is absent by default. A wall-guided optimization is
followed by a separate unbiased geomeTRIC relaxation; shell validity and final
selection use that relaxed geometry and its physical backend energy. The
search objective and wall-energy contribution are retained in candidate
diagnostics. A changed wall request cannot resume an existing run.

This is an experimental search bias, not a solvation free-energy model. A
small methane-water xTB pilot verified execution and energy separation, but
did not demonstrate improved minimum recovery; do not infer a recommended
radius or strength from that pilot. The wall is rejected for ``cluster`` mode
and is not passed to native geometry optimizers.

Use ``--wall auto`` to apply the wall conditionally during a shell search. PyAR
first performs an unbiased optimization and checks shell coverage. If that
geometry is shell-valid, no wall optimization is run. If it escaped, the
candidate is retried with a derived wall and then relaxed again without the
wall. Survivor ranking still follows the shell objective: surface coverage,
then final physical energy as a tie-breaker, then structural diversity. The
confinement trigger does not change selection policy. The radius is the
farthest sampled accessible solute-surface point from the mass-weighted solute
center plus the solvent fragment's mass-weighted COM-to-atom extent and a
0.5 Å margin. The harmonic constant is calibrated so a 0.5 Å overrun costs
3 kBT at a 298.15 K reference, yielding about 0.02266 Eh/Å². This is a soft
search-stiffness convention, not a temperature-dependent sampling model or a
claim about solvation thermodynamics. The policy and per-candidate decision
are recorded in request and candidate diagnostics. Auto mode requires an
energy-gradient backend and geomeTRIC for shell searches; for the finite
``cluster`` objective it is recorded as disabled. Explicit wall values remain
available for controlled experiments.

The geometry-derived envelope follows the general sizing principle used by
solvent-placement methods such as `ORCA SOLVATOR
<https://www.faccts.de/docs/orca/6.1/manual/contents/structurereactivity/solvator.html>`_
and `CREST QCG <https://crest-lab.github.io/crest-docs/page/documentation/qcg.html>`_.
The trigger, margin and 3 kBT calibration above are PyAR's policy choices and
have not been validated as universally optimal.
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

Whole-cluster basin exploration is available as a separate opt-in v2 stage::

   pyar microsolvate methane.xyz water.xyz --count 4 --objective cluster \
     --placement-score energy --backend xtb --geometry-optimizer geometric \
     --resampling-moves 40 --resampling-evaluation-budget 2000 \
     --resampling-time-limit 7200 --resampling-seed 17

Moves rigidly translate and rotate one solvent fragment or a small group of
fragments, and occasionally exchange the local rigid geometries at two solvent
sites; each proposal then reoptimizes the complete fixed-composition cluster. The
original solute atom order and solvent-fragment partition are checked before a
candidate can enter the survivor pool. A van der Waals overlap screen rejects
severe interfragment clashes before optimization. Proposals, PCG64 state,
physical backend evaluation reservations, and survivor snapshots are persisted
under ``resampling/``. The evaluation cap is strict; the time limit is checked
before each backend energy/gradient call, so one already-running call can finish
after the deadline. Use a fresh v2 run: MS4 does not change legacy continuation.
These stochastic trials are bounded basin exploration and do not establish
global-minimum recovery.

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
