PyAR 2.0 API Migration
======================

PyAR 2.0 narrows the supported Python API to documented package entrypoints.
Legacy compatibility modules from the 1.x transition remain as thin shims
while the migration is completed. Command-line workflows continue to use
``pyar-cli``.

Supported Public Imports
------------------------

* ``from pyar import Molecule, __version__``
* ``from pyar.workflows import aggregate, react, solvate``
* ``from pyar.sampling import generate_trial_vectors``
* ``from pyar.selection import choose_geometries``
* ``from pyar.backends import get_backend_capabilities``
* ``from pyar.biases import afir`` or ``from pyar.biases.afir import isotropic``

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

``pyar-cli`` remains the primary command-line entry point. Auxiliary commands
that point to maintained script modules remain available, including
``pyar-reaction-trace``, ``pyar-benchmark-orientations``,
``pyar-benchmark-clustering``, ``pyar-descriptor``, ``pyar-mlopt``, and
``pyar-aimnet2-ase-opt``.
