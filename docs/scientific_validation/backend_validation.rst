Backend validation matrix
=========================

The machine-readable matrix is
``benchmarks/scientific_validation/backend_matrix.json``. A test requires its
backend entries and prerequisite metadata to match
``pyar/backend_capabilities.py``. The capability registry and
missing-dependency checks establish routing and diagnostics, not numerical
chemistry. The matrix records evidence and limitations for xTB, ORCA, Gaussian,
Turbomole, Psi4, MOPAC, OpenBabel, AIMNet2, MLatom, ANI and the hybrid routes.
ASE, geomeTRIC and Sella are Python optimizer/integration components rather
than entries in the electronic-structure backend registry.

For a real-backend result, record executable/library/model version, command and
settings, charge, multiplicity, ordered symbols, input hash, output/logs,
energy and gradient units, optimizer convergence evidence and candidate SHA.
Mocks and provider fixtures remain interface tests. Proprietary software is
not required in public CI; unavailable backends are ``NOT_TESTED``.
