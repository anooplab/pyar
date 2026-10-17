Bond Scan
=========

``scan-bond`` performs a relaxed surface scan along one distance between
two molecular fragments. Atom indices are zero-based and local to each input
fragment. The scan endpoint defaults to 0.8 times the sum of the selected
atoms' PyAR covalent radii, and the default spacing is 0.10 Angstrom.

.. code-block:: bash

   pyar-cli scan-bond A.xyz B.xyz --atoms 0 1 --software orca -N 4

ORCA remains the default for ``--software`` and uses its native relaxed scan.
Other registered energy/gradient backends (currently standalone xTB, Gaussian,
and AIMNet2) use constrained geomeTRIC optimization at each scan distance,
starting each point from the previous optimized geometry. Unsupported backends
are rejected; there is no silent replacement of the requested calculator.

Select an ORCA method with ``--method``. DFT methods require
an explicit basis set, while ORCA's built-in xTB method keywords do not take
one:

.. code-block:: bash

   pyar-cli scan-bond A.xyz B.xyz --atoms 0 1 -N 4 --method BP86 --basis def2-SVP
   pyar-cli scan-bond A.xyz B.xyz --atoms 0 1 -N 4 --method GFN2-xTB

ORCA recognizes aliases such as ``XTB2`` for ``GFN2-xTB`` and ``GFN-xTB``
for GFN1-xTB. PyAR omits DFT-only keywords (basis, RI, D3BJ, and KDIIS) for
these methods. The selected physical model is used throughout the requested calculation.
Gaussian requires ``--method`` and ``--basis``. AIMNet2 uses its existing
PyAR model. Standalone xTB selects its Hamiltonian using ``--xtb-model``:

.. code-block:: bash

   pyar-cli scan-bond A.xyz B.xyz --atoms 0 1 -N 4 --software xtb

``gfn2`` is the standalone xTB default. Both models use the normal
``xtb`` executable; no equivalence with a different implementation is claimed.

Set fragment charges, multiplicities, and SCF types with ``--charge``,
``--multiplicity``, and ``--scftype`` (one value applies to both fragments or
provide one value per fragment). The combined A+B spin normally follows PyAR's
historical multiplicity heuristic. For open-shell fragments this coupling is
ambiguous, so select the total spin explicitly with ``--total-multiplicity``;
the value is checked against the allowed angular-momentum couplings, recorded
in the request and orientation provenance, and carried through free relaxation
and reaction-path stages:

.. code-block:: bash

   pyar-cli scan-bond radical_a.xyz radical_b.xyz --atoms 0 0 --software orca \
     --method BP86 --basis def2-SVP --multiplicity 2 --scftype uhf \
     --total-multiplicity 1 --through all

Omitting ``--total-multiplicity`` preserves the legacy heuristic and emits a
warning for open-shell fragment coupling. ORCA DFT and built-in xTB gradient
inputs use method-specific keywords; DFT basis/RI/dispersion/SCF settings are
not applied to built-in xTB routes. ORCA external g-xTB remains native-scan-only
because its current wrapper route does not provide PyAR with a supported
energy/gradient interface for continuation.

ORCA's g-xTB support uses its external-method wrapper rather than a built-in
method keyword. Supply the executable ``oet_gxtb`` wrapper path with
``--gxtb-wrapper`` (the wrapper and g-xTB parameter files must be installed):

.. code-block:: bash

   pyar-cli scan-bond A.xyz B.xyz --atoms 0 1 -N 4 --method g-xTB \
     --gxtb-wrapper /path/to/oet_gxtb

The ORCA 6.1 tutorial describes this g-xTB interface as preliminary and
Linux-only; it uses numerical gradients, so scans can be substantially slower
than native GFN-xTB methods. This external ORCA wrapper is supported for native
scans only. To continue a g-xTB scan into the reaction-path cycle, select
``--software xtb --xtb-model gxtb``. ORCA built-in GFN-xTB methods support
continuation through the ORCA energy/gradient provider. ORCA xTB execution
requires its external xTB installation, for example the documented
``XTBEXE`` environment setting.

Use ``--scan-end`` to set an explicit endpoint and either ``--scan-step`` or
``--scan-points`` to control the scan grid. ``--opt-threshold`` selects the
scan optimizer's convergence preset: native ORCA for ORCA scans and geomeTRIC
for generic scans. Those presets are not numerically identical across
optimizers. The CLI defaults to scan only.
Explicitly request continuation using ``--through``. The Python API preserves
its historical scan-plus-free-relaxation behavior when ``through`` is omitted;
pass ``through="scan"`` for a scan without free relaxation.

For one reaction-path characterization per unique relaxed product, opt in with
``--characterize-unique-products`` and select the last path stage with
``--through``. This mode completes the scan and free product relaxation first,
groups equivalent products, then seeds one representative path from that
orientation's complete scan trajectory. Duplicate orientations remain recorded
as alternative routes. For example:

.. code-block:: bash

   pyar scan-bond A.xyz B.xyz --atoms 0 0 --backend xtb --orientations 5 \
     --through all --characterize-unique-products --output scan_bond

Without this flag, ``--through`` keeps its existing per-orientation behavior.
The opt-in requires a fresh or v2 output directory and a stage beyond ``scan``.
Changing path-stage settings reuses verified scan and relaxation work while
creating a distinct path attempt. A scan-derived relaxed product identity does
not itself confirm a reaction; only validated endpoint/IRC evidence sets the
``reaction_connected`` product status.

Results are written to ``scan_bond/`` by default. The directory contains the
request and summary JSON/CSV files, per-orientation starting structures, ORCA
input/output, the complete scan trajectory, and the final scan frame.
Continued runs also contain standardized reaction-path artifacts and summaries
in each orientation's ``reaction_path/`` directory. Re-running the identical request in the same
output directory reuses completed work; a different request is rejected so
results cannot be silently mixed.

The target contact after relaxation and the product identity are separate
diagnostics. Contact presence uses a distance threshold based on the selected
atoms' covalent radii (1.3 times their sum; the threshold is recorded in each
orientation result); product identity compares canonical molecular identity
before the scan with the freely relaxed structure. Neither alone establishes
a reaction mechanism.

Each successful scan also writes ``scan/scan_profile.csv`` and
``scan/scan_profile.json`` from ORCA's ``scan.relaxscanact.dat`` actual-energy
table or the generic scan's final point energies. The energy source is recorded. ``ts_candidates/`` contains the highest-energy scan geometry and its
adjacent frames when available. These are scan-maximum candidates only: a
relaxed distance scan does not locate or confirm a transition state, and its
maximum may occur at an endpoint. Inspect the profile and optimize/validate
any proposed transition structure independently.

The default endpoint factor (0.8 times the covalent-radius sum) is a practical
compression heuristic, not a universal chemical threshold. Use ``--scan-end``
to specify a system-appropriate endpoint explicitly.

An ORCA endpoint-factor check over five small association examples is recorded
in ``benchmarks/scan_bond/VALIDATION.md``. It is a coarse implementation
validation, not a universal calibration of the default factor.


Continuing the scan
-------------------

Use one cumulative stopping option rather than combinations of independent
flags. Every mode includes the scan and the prerequisites of its final stage.

.. list-table:: Cumulative modes
   :header-rows: 1
   :widths: 25 75

   * - ``--through``
     - Calculation
   * - ``scan`` (default)
     - Constrained relaxed scan only.
   * - ``neb``
     - Scan, free relaxation of reactant/product references, then climbing-image NEB.
   * - ``ts``
     - Previous stages and TS optimization from the final NEB candidate.
   * - ``frequency``
     - Previous stages and independent TS frequency/stationarity validation.
   * - ``irc``
     - Previous stages and both IRC directions.
   * - ``endpoints``
     - Previous stages and free optimization of both IRC endpoints.
   * - ``endpoint-frequency``
     - Previous stages, both endpoint frequency checks, and reaction connection validation.
   * - ``all``
     - Alias for ``endpoint-frequency``: the full validated cycle.

NEB includes climbing-image support; there is no separate "NEBTS" stage.
The highest-energy interior scan frame is a waypoint initializer. If a scan
has only two points, their coordinate midpoint supplies the waypoint. A scan
maximum at an endpoint is recorded honestly and never promoted to a confirmed
TS. The final valid geomeTRIC climbing image is preferred for TS refinement;
the highest-energy interior NEB image is the fallback. Scan waypoints and NEB
candidates require independent TS frequency and reaction connection validation.

.. code-block:: bash

   pyar-cli scan-bond A.xyz B.xyz --atoms 0 1 -N 4 \
     --software xtb --xtb-model gfn2 --through neb --output reaction_scan

   # Extend the same scan to NEB or to the full cycle without rescanning:
   pyar-cli scan-bond A.xyz B.xyz --atoms 0 1 -N 4 \
     --software xtb --xtb-model gfn2 --through all --output reaction_scan

``--ts-optimizer geometric`` remains the default; ``sella`` is optional.
``--ts-fmax`` sets the existing shared TS force threshold. Sella internal
coordinates are explicitly selected with ``--sella-internal-coordinates``.
The NEB interpolation, climbing threshold and convergence options retain the
``pyar-neb`` meanings; ``--neb-max-cycles`` controls NEB cycles.

Each orientation stops at the first exception or failed scientific gate.
The CLI reports the overall result and exits nonzero if any orientation fails.
A stage recorded as complete may still fail its scientific gate; both are
reported. Optimized IRC endpoints without frequency checks are not certified
minima and cannot confirm the reaction connection. Neither a scan maximum,
a climbing image nor optimizer convergence establishes a transition state.
The final connection requires a stationary first-order saddle, converged IRCs,
two verified distinct endpoint minima, and agreement with the intended endpoints.

Changing ``--through`` does not change scan identity. Completed reaction
stages are reused only after validating physical settings, stage parameters,
input/artifact hashes and upstream summary hashes. Extending to endpoint
frequencies reuses the optimized endpoint geometries without optimizing them
again or rewriting their upstream artifacts. Changing the Hamiltonian requires
a new scan output directory; changing downstream convergence settings reruns
the affected reaction stages and invalidates dependent stages.

Run layout v2
-------------

Fresh scans use explicit ``layout.json`` and immutable copied inputs. Each
``orientation_NNN`` retains ``start.xyz`` and ``result.json``; native scan,
profile, candidate and free-relaxation files live together under
``attempts/attempt_000/execution``. A successful scan publishes
``scan_commit.json`` with request identity and hashed relative artifact refs.
Changing ``--through`` does not change that identity. Corrupt committed files
are rejected rather than silently rescanned. Free relaxation has a separate
``relax_commit.json``; mutable reports are not restart authority.

Reaction-path continuation copies its inputs into ``reaction_path/inputs``.
Native cumulative stage files share ``reaction_path/attempts/attempt_NNN/execution``;
per-stage commits under ``reaction_path/stages`` index their provenance and
scientific gates. Changing continuation settings selects another attempt and
preserves prior commits. Native stage names and scientific reuse checks remain
unchanged. Existing v1 runs continue with their original paths.

``results/manifest.json`` indexes scan-maximum candidates in their native
locations. These candidates are explicitly **not validated transition states**.

When free relaxation has completed, layout v2 also writes
``products/index.json``. It groups only orientations whose relaxed molecular
identity was successfully established and differs from the separated
reactants. The product record points to the hashed scan and relaxation
artifacts for every route, and marks one deterministic representative while
retaining alternative orientations. An identity failure remains an individual
unknown entry; an unchanged reactant is listed as a non-product. Target-bond
status is recorded per route, so an off-target product remains visible.

If a failed orientation succeeds on restart, its route can join an existing
product without changing that product's ID. Earlier ``product.json`` records
remain intact; additional membership snapshots live under the product's
``revisions`` directory, and the index names and hashes the current snapshot.
If the selected representative changes, its characterization uses a separate
directory beneath ``characterization/representatives/orientation_NNN`` so the
previous representative's inputs and path evidence remain intact.

Representative ranking uses a complete validated scan trace, an available
internal scan maximum, the constrained scan peak relative to its first point,
and finally orientation number. That peak is only a route-selection heuristic;
it is not an activation barrier or transition-state energy. Registry creation
does not start NEB or other reaction-path stages. Existing v1 output and the
default ``--through`` behavior are unchanged.
