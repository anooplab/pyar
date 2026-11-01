AIMNet2 Backend
===============

AIMNet2 is a machine-learning-potential backend available in PyAR. It is useful
when you want fast molecular energy and force evaluations without calling an
external quantum-chemistry executable.

When to use AIMNet2
-------------------

Use ``aimnet_2`` when you want:

* fast optimisation of molecular structures
* a machine-learning-potential route for aggregation, solvation, or reaction
  exploration
* AFIR/geomeTRIC reaction searches without relying on Turbomole coordinate
  updates

Example commands
----------------

Molecular cluster search:

.. code-block:: bash

   pyar microsolvate water.xyz water.xyz --backend aimnet_2 --count 10 --orientations 16 --charge 0 0 --multiplicity 1 1

Reaction search:

.. code-block:: bash

   pyar react A.xyz B.xyz --orientations 8 --bias-max 1000 --backend aimnet_2

Practical notes
---------------

* AIMNet2 model assets are not bundled into the ``pyar-chem`` wheel. If you
  rely on them outside a source checkout, provide the model files separately
  from the upstream project or another model-package source.
* AIMNet2 is a third-party project from the Isayev Lab and is MIT licensed
  upstream. See ``THIRD_PARTY_LICENSES/AIMNet2-LICENSE`` and
  ``THIRD_PARTY_LICENSES/AIMNet2-PROVENANCE.md``.
* As with all ML potentials, use AIMNet2 aggressively for exploration but
  validate important structures and energetics with an appropriate higher-level
  method.

AFIR/geomeTRIC status
---------------------

``aimnet_2`` is registered as an energy-gradient provider for the
geomeTRIC-backed AFIR reaction route.
