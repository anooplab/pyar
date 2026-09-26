Bond Scan
=========

``scan-bond`` is an ORCA-only relaxed surface scan along one distance between
two molecular fragments. Atom indices are zero-based and local to each input
fragment. The scan endpoint defaults to 0.8 times the sum of the selected
atoms' PyAR covalent radii, and the default spacing is 0.10 Angstrom.

.. code-block:: bash

   pyar-cli scan-bond A.xyz B.xyz --atoms 0 1 --software orca -N 4

Use ``--scan-end`` to set an explicit endpoint and either ``--scan-step`` or
``--scan-points`` to control the scan grid. Each successful scan is followed
by an unconstrained ORCA optimization from its final constrained geometry.

Results are written to ``scan_bond/`` by default. The directory contains the
request and summary JSON/CSV files, per-orientation starting structures, ORCA
input/output, the complete scan trajectory, the final scan frame, and the
unconstrained relaxed structure. Re-running the identical request in the same
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
table. ``ts_candidates/`` contains the highest-energy scan geometry and its
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
