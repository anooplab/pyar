PyAR 2.0 API Migration
======================

Use the documented package entrypoints for new Python code. Compatibility
modules and the legacy ``pyar-cli`` command remain available during the 1.x
transition; this guide does not announce their removal. New calculations should
use the task-oriented ``pyar COMMAND`` interface.

Supported Public Imports
------------------------

* ``from pyar import Molecule, __version__``
* ``from pyar.workflows import aggregate, grow, conformer_search, react``
* ``from pyar.workflows import solvate, microsolvate, run_scan_bond``
* ``from pyar.workflows import WorkflowResult, AggregateResult, GrowResult``
* ``from pyar.sampling import generate_trial_vectors``
* ``from pyar.selection import choose_geometries``
* ``from pyar.backends import get_backend_capabilities``
* ``from pyar.biases import afir`` or ``from pyar.biases.afir import isotropic``

Feature packages own request models and restart-contract payloads for their
domains. Use workflow functions to run calculations; use request models when
validating or inspecting workflow configuration:

* ``pyar.aggregation.AggregateRequest``
* ``pyar.solvation.SolvationRequest``
* ``pyar.conformer.ConformerRequest``

Import Path Changes
-------------------

.. list-table::
   :header-rows: 1
   :widths: 35 40 25

   * - Old import path
     - New import path
     - Notes
   * - ``pyar.Molecule``
     - ``pyar`` or ``pyar.core.molecule``
     - Use ``from pyar import Molecule`` for application code.
   * - ``pyar.molecule_geometry``
     - ``pyar.core.geometry``
     - Internal geometry helpers moved under ``core``.
   * - ``pyar.molecule_io``
     - ``pyar.io.xyz``
     - XYZ parsing and writing live under ``io``.
   * - ``pyar.orientation_sampling``
     - ``pyar.sampling``
     - Public sampling export is ``generate_trial_vectors``.
   * - ``pyar.trial_generation``
     - ``pyar.sampling.trial_generator``
     - Internal trial-generation helpers live in the sampling package.
   * - ``pyar.reactor``
     - ``pyar.workflows.reaction``
     - Use ``from pyar.workflows import react`` for application code.
   * - ``pyar.aggregator``
     - ``pyar.workflows.aggregate``
     - Use ``from pyar.workflows import aggregate`` for application code.
   * - ``pyar.interface``
     - ``pyar.backends``
     - Backend adapters live under ``backends``; legacy aliases may remain for compatibility.
   * - ``pyar.afir.restraints``
     - ``pyar.biases.afir``
     - AFIR bias functions are implemented under ``biases``.
   * - ``pyar.afir.utils``
     - ``pyar.biases.afir``
     - ``resolve_gamma`` moved with AFIR bias helpers.
   * - ``pyar.data_analysis.clustering``
     - ``pyar.selection``
     - Use ``from pyar.selection import choose_geometries``.

Console Scripts
---------------

``pyar`` is the task-oriented command-line entry point for new workflows.
Examples include ``pyar grow``, ``pyar microsolvate``, ``pyar react`` and
``pyar scan-bond``. ``pyar solvate`` is a compatibility alias for modern
microsolvation. The legacy ``pyar-cli solvate`` command retains its historical
algorithm and restart state. ``pyar-cli`` and auxiliary commands remain
available, including
``pyar-reaction-trace``, ``pyar-benchmark-orientations``,
``pyar-benchmark-clustering``, ``pyar-descriptor``, ``pyar-mlopt``, and
``pyar-aimnet2-ase-opt``.

Use ``pyar install status`` to inspect installed workflow dependencies and
``pyar react A.xyz B.xyz --backend xtb --bias-max 100 --through all`` to request
the full validated reaction path. ``--check`` and ``--dry-run`` validate the
requested stages without starting calculations.

Release boundary
----------------

The ``v1.5.1`` tag does not include the run-layout-v2 changes on the current
``master`` branch. At the 2026-10-16 release check, ``master`` still declared
version ``1.5.1`` and had no corresponding release tag. Check the
installed package version and release notes before assuming that a tagged
release contains current ``master`` behavior. Existing recognized v1 runs
continue on their legacy path; they are not silently converted to v2.
