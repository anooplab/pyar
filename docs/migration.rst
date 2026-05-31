Migration Guide
===============

PyAR keeps a narrow public API under the canonical package paths listed below.
Legacy import paths remain available as compatibility shims during the
cleanup period, but new code should use the canonical modules.

Canonical imports
-----------------

* ``from pyar import Molecule``
* ``from pyar.workflows import aggregate, react, solvate``
* ``from pyar.sampling import generate_trial_vectors``
* ``from pyar.selection import choose_geometries``
* ``from pyar.backends import ...`` for backend adapters and capability helpers
* ``from pyar.biases.afir import isotropic``

Common path changes
-------------------

.. list-table::
   :header-rows: 1
   :widths: 36 40 24

   * - Legacy path
     - Canonical path
     - Notes
   * - ``pyar.Molecule``
     - ``pyar`` or ``pyar.core.molecule``
     - Use ``from pyar import Molecule`` for application code.
   * - ``pyar.molecule_geometry``
     - ``pyar.core.geometry``
     - Geometry helpers now live under ``core``.
   * - ``pyar.molecule_io``
     - ``pyar.io.xyz``
     - XYZ parsing and writing live under ``io``.
   * - ``pyar.orientation_sampling``
     - ``pyar.sampling``
     - The sampling package owns the public entrypoints.
   * - ``pyar.trial_generation``
     - ``pyar.sampling.trial_generator``
     - Trial generation is implemented in the sampling package.
   * - ``pyar.aggregator``
     - ``pyar.workflows.aggregate``
     - Use ``pyar.workflows.aggregate`` for new code.
   * - ``pyar.reactor``
     - ``pyar.workflows.reaction``
     - Use ``pyar.workflows.reaction`` for new code.
   * - ``pyar.interface``
     - ``pyar.backends``
     - Backend adapters live under ``backends``.
   * - ``pyar.afir.restraints``
     - ``pyar.biases.afir``
     - AFIR bias functions live under ``biases``.
   * - ``pyar.data_analysis.clustering``
     - ``pyar.selection``
     - Selection services are now canonical under ``pyar.selection``.

See also
--------

* :doc:`api`
* :doc:`architecture`
