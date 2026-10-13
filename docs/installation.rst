Installation
============

What do you want to do?
-----------------------

PyAR is distributed as ``pyar-chem`` and imported as ``pyar``. Use Python
3.10 or newer. Install the workflow you need:

.. code-block:: bash

   # Aggregation with an externally installed backend such as xTB
   python -m pip install "pyar-chem[aggregate]"

   # Conformer generation and AIMNet2 refinement
   python -m pip install "pyar-chem[conformer,aimnet2]"

   # Reaction exploration with an externally installed backend such as ORCA
   python -m pip install "pyar-chem[reaction]"

   # Explicit solvent shells or reaction paths
   python -m pip install "pyar-chem[microsolvate]"
   python -m pip install "pyar-chem[neb,sella]"

The base ``python -m pip install pyar-chem`` installs only shared molecule,
geometry, I/O and CLI dependencies: NumPy, SciPy, NetworkX, and tomli on Python
3.10. Workflow extras install their Python requirements automatically.
``grow`` shares aggregation dependencies; ``solvate`` is a compatibility alias
for ``microsolvate``. Backend selection remains a runtime option.

Check your installation
-----------------------

.. code-block:: bash

   pyar --help
   pyar-cli --help
   pyar doctor
   pyar install status
   pyar install status --json
   pyar aggregate molecule.xyz --size 2 --backend xtb --check
   pyar react A.xyz B.xyz --backend orca --method BP86 --basis def2-SVP --bias-max 100 --check

Normal execution automatically checks requirements. ``--check`` validates
without calculations. Missing optional components in ``doctor`` are
informational; missing core requirements make it fail.

``pyar install status`` lists package versions, available workflow Python
stacks, external executables on PATH and the AIMNet2 model file. It shows
the interpreter and environment being inspected. This command is read-only;
it does not install packages or launch external programs. Status reports
metadata and discovery; workflow ``--check`` validates the selected route.

Python capabilities
-------------------

Extras compose normally: ``[conformer,aimnet2]`` or ``[neb,sella,geodesic]``.
``geometric`` installs ASE and geomeTRIC; ``sella`` selects Sella support;
``geodesic`` supplies geodesic interpolation; ``structure-comparison`` supplies
iRMSD; ``identity`` supplies RDKit; ``analysis`` supplies legacy descriptor,
CSV and plotting utilities; ``ml`` supplies MLatom/TorchANI interfaces.
The retained TorchANI 2.0 interface uses ``pkg_resources`` at runtime, so
``ml`` also installs a compatible ``setuptools>=69,<81`` version.
``all`` contains end-user Python capabilities only and is not the recommended
installation path.

AIMNet2 model configuration
---------------------------

``aimnet2`` installs ASE and Torch. Large model files remain external. Set
``PYAR_AIMNET2_MODEL`` before starting PyAR to the legacy
``aimnet2_wb97m-d3_0.jpt`` file from https://github.com/isayevlab/AIMNet2.
The upstream ``aimnet`` distribution currently requires Python 3.11+ and uses
a different model/calculator interface; PyAR retains its working adapter to
preserve Python 3.10 and scientific behavior. See :doc:`aimnet2_backend`.

External backends
-----------------

PyAR does not install xTB, ORCA, Gaussian, Turbomole, MOPAC, Psi4, or OpenBabel
executables through pip. Install the selected program separately and make its
executables available on PATH. See :doc:`external_programs` for installation
sources and verification commands.

Development and migration
-------------------------

In a checkout, use pip 25.1+ for standardized dependency groups:

.. code-block:: bash

   python -m pip install --upgrade pip
   python -m pip install -e ".[aggregate,reaction,neb]" --group dev
   python -m pytest
   python -m build

The old ``test`` and ``docs`` extras move to ``--group test`` and
``--group docs``; release tools use ``--group release``. ``selection`` remains
a compatibility alias for aggregate plus analysis. ``xtb`` remains an alias
for Python optimizer support (``geometric``), and does not install xTB.
The unreliable ``openbabel`` source-build extra is removed: use the official
installer or ``mamba install -c conda-forge openbabel`` instead.
