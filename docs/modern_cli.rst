Task-oriented command interface
===============================

The ``pyar`` executable is the task-oriented entry point. Legacy commands such
as ``pyar-cli`` and the ``pyar-*`` executables remain available during the
migration.

Environment and run inspection
------------------------------

.. code-block:: console

   pyar info
   pyar backends
   pyar backends xtb
   pyar backends --json
   pyar doctor
   pyar doctor --all-backends
   pyar inspect reaction/reaction/state.json

``info`` reports the installed package and Python environment. ``backends``
lists registered capabilities without importing every backend. ``doctor``
checks the core installation; ``--all-backends`` also reports optional
executables and modules. ``inspect`` reads a PyAR run record or workflow state.

Run configurations and profiles
--------------------------------

Save a validated task request as an editable, versioned TOML configuration:

.. code-block:: console

   pyar react A.xyz B.xyz --backend xtb --bias-max 100 \
       --write-config reaction.toml --check
   pyar run reaction.toml --check
   pyar run reaction.toml

Configurations store resolved command settings in a ``[settings]`` table,
including inferred electronic states, alongside request provenance. Edit that
table or override an option with ``pyar run reaction.toml --bias-max 150``.
Replay reuses the normal scientific resolver and preflight; saved ``--check``
requests do not turn future executions into validation-only runs. Relative
input/output paths resolve against the current working directory.
``pyar config validate reaction.toml`` performs the same read-only preflight
as ``pyar run reaction.toml --check`` for calculation workflows.

Profiles store methodology settings and are applied before explicit
command-line options, so CLI values win:

.. code-block:: console

   pyar profile create careful-react --from reaction/pyar-run.toml
   pyar profile list
   pyar react C.xyz D.xyz --backend xtb --bias-max 100 --profile careful-react

Profiles live in ``$XDG_CONFIG_HOME/pyar/profiles`` (or
``~/.config/pyar/profiles``). They intentionally omit input and output paths.

Validation and dry runs
-----------------------

``--check`` asks whether the exact request can run and performs the command's
read-only validation and preflight. ``--dry-run`` performs that same validation
and, only after it succeeds, prints a shared execution-plan summary: command
options with their sources, inputs, planned stages, and known output locations.
Plans include backend-resolved settings and electronic states. Neither mode
starts chemistry or creates workflow output. Explicit ``--write-config`` writes
only the requested configuration after validation. For example:

.. code-block:: console

   pyar optimize molecule.xyz --backend xtb --dry-run

``--write-config``, ``--profile`` and ``-v`` are shared modern CLI controls.
The command-specific resolver and preflight remain authoritative for scientific
validation; the common plan does not replace workflow logic. Completed runs
and failures after the workflow starts write a ``pyar-run.toml`` record with
resolved scientific settings, electronic states, dependency checks, the workflow
result, and CLI/config/profile/default provenance. Later invocations retain
previous records; ``inspect`` selects the newest record in a run directory.
Validation-only requests do not create run records. ``--json`` on calculation
commands emits a single structured result; ordinary progress goes to stderr.
``-q`` suppresses progress and ``-v``/``-vv`` enable diagnostic logging without
changing scientific settings.

Reaction-path tasks
-------------------

The modern path-stage commands use the existing validated ``pyar.neb`` engine:

.. code-block:: console

   pyar neb reactant.xyz product.xyz --ts-guess guess.xyz --backend xtb --check
   pyar ts guess.xyz --backend xtb --check
   pyar irc neb_run --backend xtb --check

ORCA and Gaussian requests require explicit method/basis settings when needed.
``neb`` relaxes the endpoints and runs NEB, stopping before TS optimization.
``ts`` optimizes the candidate and checks its frequency signature. ``irc``
consumes a frequency-validated TS stage directory, verifying saved physical
settings and artifact hashes before starting. It runs both IRC directions;
endpoint optimization/validation remains a separate engine/scan continuation
stage. Normal numerical and input/output-alias checks happen during preflight.
The scan-bond ``--through`` workflow continues to use the same stage engine.

Compatibility spelling
----------------------

``pyar solvate`` is a compatibility alias for the canonical modern
``pyar microsolvate`` command. The legacy ``pyar-cli solvate`` syntax and old
restart-state format remain on their compatibility path.
