Scientific validation troubleshooting
======================================

Read the machine-readable JSON and human-readable Markdown report together
with the case manifest and raw backend logs. The validator fails when the
candidate method, backend, electronic state, input hash, backend version,
metric, tolerance or log reference is missing or inconsistent.

.. list-table:: Failure classes
   :header-rows: 1

   * - Class
     - First evidence to inspect
   * - Scientific regression
     - Expected/observed metrics, geometry, identity and pathway validation
   * - Numerical tolerance failure
     - Units, absolute/relative error, precision, threads and frozen tolerance
   * - Backend incompatibility
     - Executable/version, generated input, raw stdout/stderr and parser output
   * - Infrastructure failure
     - Runner status, resource limits, filesystem and retained stage logs
   * - Missing backend / NOT_TESTED
     - Dependency/executable/model availability; never treat as PASS
   * - Unsupported workflow
     - Capability registry and documented implementation status
   * - Nondeterministic result
     - Seed, scheduling, hardware, thread count and repeated raw results
   * - Reference inconsistency
     - Provenance, immutable hashes, source version and mapping

Repair the supported cause and rerun the same frozen case. Do not widen a
tolerance or update a reference solely to turn a failing result green.
