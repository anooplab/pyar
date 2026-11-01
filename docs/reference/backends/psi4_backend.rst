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

   pyar optimize structure.xyz -c 0 -m 1 --backend psi4 --method B97-D --basis def2-SVP

Practical notes
---------------

* Psi4 must be installed separately.
* ``--method`` and ``--basis`` select the Psi4 method and basis; they default to
  ``B97-D`` and ``def2-SVP``. Psi4 attaches def2 effective core potentials itself.
* Open-shell molecules (multiplicity above 1, or ``--scftype uhf``) run with
  ``set reference uhf``. Psi4 derives the HF/DFT functional from the requested
  method; this reference also supports open-shell correlated methods such as CCSD.
* The current capability registry does not mark ``psi4`` as a registered
  energy-gradient provider for the geomeTRIC AFIR route.
* Use a small validation case before relying on Psi4 outputs in a larger PyAR
  workflow.

AFIR/geomeTRIC status
---------------------

``psi4`` is not currently registered as an energy-gradient provider for the
geomeTRIC-backed AFIR reaction route.
