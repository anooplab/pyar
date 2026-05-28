Psi4 Backend
============

Psi4 is listed as a quantum-chemistry backend family in PyAR. It can be useful
for users who have a local Psi4 installation and want to connect PyAR-generated
structures to Psi4 calculations.

When to use Psi4
----------------

Use ``psi4`` when:

* Psi4 is installed and configured on your system
* you specifically want to test a Psi4-based route
* you are comfortable validating backend behaviour on small examples before
  larger PyAR runs

Example command
---------------

.. code-block:: bash

   pyar-optimiser structure.xyz -c 0 -m 1 --software psi4

Practical notes
---------------

* Psi4 must be installed separately.
* The current capability registry does not mark ``psi4`` as a registered
  energy-gradient provider for the geomeTRIC AFIR route.
* Use a small validation case before relying on Psi4 outputs in a larger PyAR
  workflow.

AFIR/geomeTRIC status
---------------------

``psi4`` is not currently registered as an energy-gradient provider for the
geomeTRIC-backed AFIR reaction route.
