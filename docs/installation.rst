Installation
============

PyAR is distributed on PyPI as ``pyar-chem``. The import namespace remains
``pyar`` and the command-line entry point remains ``pyar-cli``.

Install
-------

From PyPI:

.. code-block:: bash

   python -m pip install pyar-chem

From a local checkout:

.. code-block:: bash

   python -m pip install -e .

Optional Python extras
----------------------

* ``selection`` for clustering and selection helpers
* ``ml`` for MLatom, TorchANI, and related ML interfaces
* ``aimnet2`` for AIMNet2 runtime support
* ``openbabel`` for the Open Babel Python binding
* ``test`` for the test and build toolchain
* ``docs`` for Sphinx documentation builds
* ``all`` for the combined optional dependency set

Check the install
-----------------

.. code-block:: bash

   pyar-cli --help

External programs
-----------------

Many workflows also need external executables such as xTB, ORCA, Gaussian,
Psi4, MOPAC, TURBOMOLE, Open Babel, MLatom, or DFT-D4. See
:doc:`external_programs` for the official upstream sites.

See also
--------

* :doc:`quickstart`
* :doc:`external_programs`
