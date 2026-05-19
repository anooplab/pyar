Installation
============

Install PyAR from a local checkout:

.. code-block:: bash

   python -m pip install .

This installs the package metadata, Python dependencies, and the bundled
runtime assets used by the MLatom and AIMNet2 interfaces. The packaging
metadata now carries the dependency set, so users do not need to install
``hdbscan`` or ``DBCV`` manually for a normal pip-based install.

For development, an editable install works the same way:

.. code-block:: bash

   python -m pip install -e .

For documentation builds on Read the Docs, the project uses a minimal Sphinx
configuration and installs the package from the repository itself.
