Bond Scan
=========

``scan-bond`` is an ORCA-only relaxed surface scan along one distance between
two molecular fragments. Atom indices are zero-based and local to each input
fragment. The scan endpoint defaults to 0.8 times the sum of the selected
atoms' PyAR covalent radii, and the default spacing is 0.10 Angstrom.

.. code-block:: bash

   pyar-cli scan-bond A.xyz B.xyz --atoms 0 1 --software orca -N 4

ORCA is currently the only scan-bond backend and is the default for
``--software``. Select an ORCA method with ``--method``. DFT methods require
an explicit basis set, while ORCA's built-in xTB method keywords do not take
one:

.. code-block:: bash

   pyar-cli scan-bond A.xyz B.xyz --atoms 0 1 -N 4 --method BP86 --basis def2-SVP
   pyar-cli scan-bond A.xyz B.xyz --atoms 0 1 -N 4 --method GFN2-xTB

ORCA recognizes aliases such as ``XTB2`` for ``GFN2-xTB`` and ``GFN-xTB``
for GFN1-xTB. PyAR omits DFT-only keywords (basis, RI, D3BJ, and KDIIS) for
these methods. The same method is used for scan and final relaxation. Other
scan-bond backends are not available yet.

ORCA's g-xTB support uses its external-method wrapper rather than a built-in
method keyword. Supply the executable ``oet_gxtb`` wrapper path with
``--gxtb-wrapper`` (the wrapper and g-xTB parameter files must be installed):

.. code-block:: bash

   pyar-cli scan-bond A.xyz B.xyz --atoms 0 1 -N 4 --method g-xTB \
     --gxtb-wrapper /path/to/oet_gxtb

The ORCA 6.1 tutorial describes this g-xTB interface as preliminary and
Linux-only; it uses numerical gradients, so scans can be substantially slower
than native GFN-xTB methods.

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
