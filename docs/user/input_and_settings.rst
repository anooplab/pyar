Inputs, charge, spin and computational settings
===============================================

XYZ files and atom identity
----------------------------

Coordinates are in angstrom and atom indices in PyAR constraints are zero-based.
Keep ordered elements consistent between pathway endpoints. The same formula
with a different atom mapping is not automatically a valid NEB endpoint pair.
XYZ comments are not a reliable electronic-state specification. Supply charge
and multiplicity explicitly when they matter.

Charge and multiplicity
-------------------------

Charge is an integer total charge; multiplicity is 2S+1, not the number of
unpaired electrons. Omitted spin is inferred from electron parity on supported
request paths. Parity rejects impossible combinations but cannot choose the
chemically correct state among several plausible spins. For multiple fragments,
check how individual states are combined and supply a total multiplicity where
the command supports it. Native optimizers receive the state of the molecule
being optimized, including a merged aggregate.

.. code-block:: bash

   pyar aggregate sodium.xyz water.xyz --size 1 1 --charge 1 0 --multiplicity 1 1 --backend xtb --orientations 1 --check

The sodium input is downloadable with :doc:`aggregate`. This example is a
closed-shell +1 aggregate. A geometry file alone does not carry all this state.

Backend and optimizer
-----------------------

Choose ``--backend`` from the capabilities in ``pyar backends``. Backend means
the energy/gradient provider; a geometry optimizer is the algorithm that uses
those values. Native and geomeTRIC routes have different capabilities. GFN2-xTB
and g-xTB are different Hamiltonians. Requesting g-xTB must not silently run
GFN2-xTB. A method/basis meaningful to ORCA is not automatically meaningful to
an xTB route. Backend-specific settings are explained in :doc:`/reference/backends`.

For DFT calculations explicitly record method, basis, dispersion, solvent,
charge, multiplicity, numerical controls and program version. An extra dispersion
term must be justified for your method. Model files and machine-learning backend
versions are part of scientific provenance.
