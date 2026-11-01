Diagnose failed or uncertain calculations
=========================================

Start with the first failure
-----------------------------

Read the command's error, ``pyar inspect`` output and the native backend log in
the named job directory. Distinguish missing dependencies, invalid requests,
backend execution failure, optimizer nonconvergence and failed scientific gates.
A process exit code alone cannot establish convergence or chemical validity.

* Missing executable: activate the intended environment and check ``PATH`` with
  ``pyar install status``. Python extras and external executables are separate.
* Invalid charge/spin: inspect atom counts, fragment charges and intended state;
  do not alter multiplicity merely to silence a diagnostic.
* Nonconverged optimization: inspect the backend's actual iterations and final
  force evidence. Retain the failed attempt before changing a budget or method.
* Imaginary modes: inspect mode vectors, geometry and thresholds in the recorded
  protocol. A TS has a different acceptance target from a minimum.
* Stale/tampered restart: use a fresh directory or recover only through the
  documented workflow mechanism. Do not edit hashes or status flags.
* Unsupported older format: use the version/environment recorded with that
  calculation. PyAR 1.6 does not migrate older runs automatically.

Before an expensive retry
----------------------------

Use ``--check`` or ``--dry-run`` on the supported calculation command. Preserve
input state and methodology; record any intentional change as a new calculation.
Trial coordinates are regenerated for each optimization attempt, even when a
job name is reused. A valid cache still requires matching provenance and input.

When reporting a bug, include the PyAR/backend versions, exact command, a minimal
nonprivate input, request/result files and the relevant raw error. Redact private
paths or credentials. Attach the earliest failing evidence rather than only the
last progress message. See :doc:`/developer/contributing` for a focused reproducer.
