Using PyAR 1.6 after an older version
=====================================

PyAR 1.6 requires Python 3.11 or newer; Python 3.10 is no longer supported.
Use a fresh environment with a supported interpreter.

PyAR 1.6 is an intentionally breaking cleanup. The single executable is ``pyar``.
Standalone commands such as ``pyar-react``, ``pyar-cli`` and ``pyar-neb`` and
compatibility-only Python import paths are removed. Use the corresponding
``pyar react``, ``pyar neb`` and other commands in :doc:`cli`.

Old run/config/checkpoint generations are not migrated. Use the environment and
version recorded by the original work to reproduce it. Keep that environment
separate from 1.6. A historical XYZ can be used as a new input, but its comment
line does not replace explicit method, charge/spin, mapping and provenance.
Do not relabel old JSON/TOML with the new schema to bypass validation.

Useful specialist utilities now live under subcommands; see :doc:`tools`.
Current data integrity, restart and scientific acceptance checks remain required.
PyAR 1.6 release/package publication is separate from the repository development
version described by these pages.
