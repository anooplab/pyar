Bond Scan
=========

``scan-bond`` is an ORCA-only relaxed surface scan along one distance between
two molecular fragments. Atom indices are zero-based and local to each input
fragment. The scan endpoint defaults to the sum of the selected atoms'
PyAR covalent radii, and the default spacing is 0.10 Angstrom.

.. code-block:: bash

   pyar-cli scan-bond A.xyz B.xyz --atoms 0 1 --software orca -N 4

Use ``--scan-end`` to set an explicit endpoint and either ``--scan-step`` or
``--scan-points`` to control the scan grid. Each successful scan is followed
by an unconstrained ORCA optimization from its final constrained geometry.

Results are written to ``scan_bond/`` by default. The directory contains the
request and summary JSON/CSV files, per-orientation starting structures, ORCA
input/output, the complete scan trajectory, the final scan frame, and the
unconstrained relaxed structure. A surviving target contact is reported as a
diagnostic; it does not establish a transition state or reaction mechanism.
