First run
=========

A first calculation needs an installed PyAR workflow extra and a calculation
backend. This example uses the microsolvation workflow with xTB. The xTB
executable is installed separately; see :doc:`external_programs`.

Install the required Python capability and verify the command-line interface:

.. code-block:: bash

   python -m pip install "pyar-chem[microsolvate]"
   pyar --version
   pyar doctor

Create ``water.xyz`` with the following contents:

.. code-block:: text

   3
   water
   O  0.000000  0.000000  0.000000
   H  0.758602  0.000000  0.504284
   H -0.758602  0.000000  0.504284

First check the request. It validates inputs and backend availability without
starting the calculation:

.. code-block:: bash

   pyar microsolvate water.xyz water.xyz --count 1 --backend xtb --check

Then run it by removing ``--check``:

.. code-block:: bash

   pyar microsolvate water.xyz water.xyz --count 1 --backend xtb

PyAR reports the output location and run result. Use ``pyar inspect`` with that
run directory to inspect its recorded state. The exact output files depend on
the workflow and version; selected structures are XYZ files and can be opened
in a molecular viewer.

If preflight reports a missing xTB executable, install it separately and ensure
``xtb`` is on ``PATH``. A successful preflight checks requirements and input
validity; it is not evidence that a calculation completed or that a scientific
result is valid.

Next, read :doc:`aggregate`, :doc:`react`, :doc:`bond_scan`, or :doc:`conformer`
for the task you want to run.
