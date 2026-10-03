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

   pyar-cli scan-bond A.xyz B.xyz --atoms 0 1 -N 4 --software xtb --xtb-model gfn2

``gxtb`` remains the standalone xTB default. Both models use the normal
``xtb`` executable; no equivalence with a different implementation is claimed.

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
