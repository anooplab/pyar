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

Model setup
-----------

Install the Python runtime and select your external model:

.. code-block:: bash

   python -m pip install "pyar-chem[aimnet2]"
   export PYAR_AIMNET2_MODEL="/absolute/path/to/aimnet2_wb97m-d3_0.jpt"
   pyar install status
   pyar optimize water.xyz --backend aimnet_2 --check

Use the PyAR version/environment selected in :doc:`/user/installation`; the
plain package command above selects the published version. Export the model
path before the process starts. Relative paths are resolved at startup so
workflow directory changes do not change the selected model. A successful
``--check`` establishes asset/runtime availability, not numerical qualification.
Retain the model source, version and checksum with calculation provenance.
The standalone ``pyar aimnet2-ase-opt`` utility takes an explicit model argument;
see :doc:`/reference/tools`.

Practical notes
---------------

* AIMNet2 models are not bundled in the source checkout, wheel or sdist.
  Obtain the compatible ``aimnet2_wb97m-d3_0.jpt`` TorchScript model from the
  `upstream AIMNet2 project <https://github.com/isayevlab/AIMNet2>`_ and set
  ``PYAR_AIMNET2_MODEL`` before starting PyAR. Installing the Python extra does
  not install model weights. Do not substitute a different model format or
  Hamiltonian without establishing compatibility and scientific suitability.
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
