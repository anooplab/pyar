Quickstart
==========

PyAR is a command-line tool for aggregation, reaction, solvation, and bond-scan
workflows in molecular systems.

Install
-------

.. code-block:: bash

   python -m pip install .

For development:

.. code-block:: bash

   python -m pip install -e .

Check the CLI
-------------

.. code-block:: bash

   pyar-cli --help

Useful workflows
----------------

Aggregate two one-atom fragments:

.. code-block:: bash

   pyar-cli -a C H -as 1 4 -N 8

Generate a starting aggregate from a formula:

.. code-block:: bash

   pyar-cli --aggregate --formula C5H4 -N 8

Run a reaction search:

.. code-block:: bash

   pyar-cli -r A.xyz B.xyz -N 8 -gmin 100 -gmax 1000 --software xtb

Scan a bond:

.. code-block:: bash

   pyar-cli --scan-bond 1 2 A.xyz B.xyz -N 8

Run a solvation search:

.. code-block:: bash

   pyar-cli -s solute.xyz solvent.xyz --software xtb -ss 10 -N 16

Notes
-----

* In aggregate mode, bare element symbols such as ``C`` and ``H`` are accepted.
* If ``--software`` is omitted for aggregate mode, PyAR logs that it is
  generating trial geometries only and skips quantum-chemistry optimization.
* ``--formula`` is only valid together with ``--aggregate``.
* Reaction runs with ``--software xtb`` or ``--software aimnet_2`` use
  geomeTRIC/TRIC for the AFIR-biased optimization, then use ``gamma=0.0`` for
  product relaxation.
