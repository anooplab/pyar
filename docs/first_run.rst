First Successful Run
====================

This page gives a small first calculation for a new user who wants to check
that PyAR works before learning all command-line options.

The example builds a small water dimer search. It is simple, fast, and easy to
inspect visually.

1. Create a working folder
--------------------------

.. code-block:: bash

   mkdir pyar_first_run
   cd pyar_first_run

2. Create the input file
------------------------

Create a file named ``water.xyz`` with a text editor. Put the following lines
in the file:

.. code-block:: text

   3
   water
   O  0.000000  0.000000  0.000000
   H  0.758602  0.000000  0.504284
   H -0.758602  0.000000  0.504284

3. Run a small job
------------------

Run a small solvation-style growth job with xTB:

.. code-block:: bash

   pyar-cli solvate water.xyz water.xyz --software xtb -ss 1 -N 4 -c 0 0 -m 1 1

If xTB is not installed, PyAR should report that the executable is missing.
See :doc:`installation`, :doc:`external_programs`, and :doc:`xtb` for backend
setup.

4. Inspect the output
---------------------

After the run, look for folders such as:

.. code-block:: text

   solvation/
     state.json
     state/
       geometries/
   aggregate_002/

For most users, the important files are the selected XYZ structures from the
final growth or aggregate directory. Open those XYZ files in a molecular
viewer that supports XYZ files.

5. Make an energy table
-----------------------

If you have a folder of selected structures, print their relative energies with:

.. code-block:: bash

   pyar-energy-table selected/*.xyz

If your selected structures are in another directory, use that path instead.
For example:

.. code-block:: bash

   pyar-energy-table aggregate_002/selected/*.xyz

6. What success looks like
--------------------------

A successful first run should give you:

* no missing-executable error
* one or more output XYZ files
* an energy table if the selected XYZ files contain energy information
* a ``state.json`` file recording restart and provenance information

7. What the common options mean
-------------------------------

.. list-table::
   :header-rows: 1
   :widths: 18 82

   * - Option
     - Meaning
   * - ``--software xtb``
     - Use xTB as the computational backend.
   * - ``-N 4``
     - Generate four trial structures or orientations at each placement step.
   * - ``-ss 1``
     - Add one solvent or growth unit around the starting structure.
   * - ``-c 0 0``
     - Use charge 0 for the two input fragments.
   * - ``-m 1 1``
     - Use singlet multiplicity for the two input fragments.

Next steps
----------

After this first run, try the task pages:

* :doc:`aggregate` for clusters and noncovalent complexes
* :doc:`solvate` for adding molecules or ligands around a core
* :doc:`react` for AFIR-style reaction searches
