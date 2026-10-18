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

Request-model packages
----------------------

Workflow functions remain the supported entry points for running calculations.
Feature packages own typed request models and restart-contract payloads:

* ``pyar.aggregation.AggregateRequest``
* ``pyar.solvation.SolvationRequest``
* ``pyar.conformer.ConformerRequest``

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

Packaging and installation
--------------------------

Install the workflow you need; see :doc:`installation` for examples.
Existing CLI commands remain available. The packaging migration is:

.. list-table:: Extras and dependency groups
   :header-rows: 1

   * - Previous installation
     - Current installation
   * - Base package with Autograd
     - Base package without the unused Autograd requirement
   * - ``selection``
     - Retained alias for ``aggregate`` plus ``analysis``; prefer workflow extras
   * - ``conformer``
     - Retained RDKit workflow extra
   * - ``xtb``
     - Retained Python optimizer alias for ``geometric``; xTB remains external
   * - ``aimnet2``
     - ASE/Torch extra; configure the external model with ``PYAR_AIMNET2_MODEL``
   * - ``ml``
     - MLatom/TorchANI extra, including actual eager import requirements
   * - ``identity``, ``structure-comparison``, ``geodesic``, ``sella``
     - Retained composable Python capabilities
   * - ``openbabel``
     - Removed fragile source-build extra; use upstream or conda-forge installation
   * - ``test`` and ``docs``
     - Checkout dependency groups: ``pip install --group test`` / ``--group docs``
   * - ``all``
     - End-user Python capabilities only; developer tools use dependency groups

Workflow extras and dependency groups change as capabilities evolve. The
installation guide and ``pyproject.toml`` are authoritative for the current
lists and compositions. Dependency groups require pip 25.1+ and a checkout
containing ``pyproject.toml``.

Task-oriented command interface
--------------------------------

Use ``pyar COMMAND`` for new workflows, including ``aggregate``, ``grow``,
``react``, ``scan-bond``, ``microsolvate``, ``conformer`` and ``optimize``.
``pyar install status`` reports installed workflow capabilities. The legacy
``pyar-cli`` and ``pyar-*`` entry points remain available during the 1.x
transition; this guide does not announce their removal. ``pyar solvate`` is a
modern microsolvation alias, whereas ``pyar-cli solvate`` retains its historic
algorithm and restart state. Use ``--check`` for read-only validation and
``--dry-run`` to include the requested execution stages.

This checkout declares package version |package_version|. Development-source
behavior can differ from an installed distribution. Check ``pyar --version``
and the matching release notes before relying on a feature from an unreleased
branch. Recognized v1 runs remain on their legacy continuation path and are not
silently converted.
