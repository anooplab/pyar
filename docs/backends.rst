Backends
========

PyAR uses backend adapters to evaluate energies, gradients, and in some cases
native optimisations. The canonical backend modules live under
``pyar.backends``.

Installable pieces
------------------

Some backend dependencies are Python packages that can be installed directly
with ``pip`` or through PyAR's optional extras:

* ``python -m pip install "pyar-chem[selection]"`` for ``ase``, ``dscribe``,
  ``hdbscan``, ``pandas``, ``scikit-learn``, and ``MDAnalysis``
* ``python -m pip install "pyar-chem[ml]"`` for ``torchani``, ``mlatom``,
  ``pyh5md``, and ``h5py``
* ``python -m pip install "pyar-chem[aimnet2]"`` for the AIMNet2 runtime
* ``python -m pip install "pyar-chem[xtb]"`` for ``geometric``
* ``python -m pip install "pyar-chem[openbabel]"`` for the Open Babel Python
  bindings

If you only want one package, install it directly. For example:

* ``python -m pip install ase``
* ``python -m pip install geometric``
* ``python -m pip install dscribe``
* ``python -m pip install mlatom``

Current backend families
------------------------

.. list-table::
   :header-rows: 1
   :widths: 18 24 58

   * - Backend
     - Family
     - Typical role
   * - ``xtb``
     - semiempirical
     - Fast optimisation, aggregation, solvation, and AFIR/geomeTRIC
       reaction searches.
   * - ``aimnet_2``
     - machine-learning potential
     - Fast energy/force evaluation for optimisation and AFIR/geomeTRIC
       reaction searches.
   * - ``orca``
     - quantum chemistry
     - Higher-level DFT-style backend that requires a local ORCA installation.
   * - ``gaussian``
     - quantum chemistry
     - Higher-level DFT-style backend that requires a local Gaussian installation.
   * - ``mopac``
     - semiempirical
     - Legacy semiempirical optimisation route.
   * - ``obabel``
     - Open Babel route
     - Lightweight structure preparation or optimisation.
   * - ``psi4``
     - quantum chemistry
     - Available as a backend family, but not currently on the registered AFIR route.
   * - ``xtb_aimnet2`` / ``xtb_aiqm1``
     - hybrid
     - Staged routes that combine fast pre-optimisation with refinement.

Choose a backend
-----------------

* Start with ``xtb`` for a fast first calculation.
* Use ``aimnet_2`` when you want an ML potential and the chemistry is in-scope.
* Use ``orca`` or ``gaussian`` only when the executable is installed and
  validated locally.
* Use legacy or hybrid routes only when you know why you need them.

See also
--------

* :doc:`external_programs`
* :doc:`architecture`
* :doc:`api`
