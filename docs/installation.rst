Installation
============

Install PyAR from a local checkout:

.. code-block:: bash

   python -m pip install .

This installs the package metadata, Python dependencies, and the bundled
runtime assets used by the MLatom and AIMNet2 interfaces. The packaging
metadata now carries the dependency set, so users do not need to install
``hdbscan``, ``DBCV``, or ``geomeTRIC`` manually for a normal pip-based
install. The OpenBabel Python binding is also installed as part of the package
requirements.

For development, an editable install works the same way:

.. code-block:: bash

   python -m pip install -e .

For documentation builds on Read the Docs, the project uses a minimal Sphinx
configuration and installs the package from the repository itself.

External program requirements are still separate from Python dependencies.
PyAR will report a clear error and exit if an executable backend is missing.
The backends that must be installed on the system are:

* ORCA
* Gaussian
* MOPAC
* xTB
* Turbomole
* OpenBabel executables such as ``obabel``, ``obminimize``, and ``obenergy``

``geomeTRIC`` is a Python dependency rather than an external electronic
structure program. It is used as the internal-coordinate optimizer for AFIR
reaction calculations with the ``xtb``, ``aimnet_2``, ``orca``, and
``gaussian`` energy/force providers.
