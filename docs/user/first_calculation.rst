Your first calculation: water
=============================

Prepare a working directory
----------------------------

Activate the environment from :doc:`installation`. Download
:download:`water.xyz </examples/water.xyz>` into a new working directory.
XYZ atom coordinates are in angstrom; the file contains an atom count, comment
line, and one element/coordinate row per atom. This is an illustrative geometry,
not a prequalified minimum. Water is neutral and singlet in this example.

Calculate an energy
--------------------

.. code-block:: bash

   pyar energy water.xyz --backend xtb --output first_energy --check
   pyar energy water.xyz --backend xtb --output first_energy
   pyar inspect first_energy

The first command checks prerequisites and creates no calculation output. The
second performs a fixed-geometry single-point calculation. Read
``first_energy/summary.csv`` and ``first_energy/job_water/result.json``. The
``energy_hartree`` field is an electronic energy, not an activation free energy.
The qualification ``minimum_not_assessed`` is expected: an energy calculation
has not tested the curvature or stationarity of the geometry.

Optimize the geometry
-----------------------

.. code-block:: bash

   pyar optimize water.xyz --backend xtb --check
   pyar optimize water.xyz --backend xtb

Inspect ``job_water/xtb.out`` for the native optimization and
``job_water/result_water.xyz`` for the accepted geometry. Input and accepted
coordinates are distinct artifacts. Compare the O–H distances or open the XYZ
in a molecular viewer. A successful optimization must have finite, compatible
coordinates and energy; inspect diagnostics when it fails.

Check whether the result is a minimum
--------------------------------------

.. code-block:: bash

   pyar frequency job_water/result_water.xyz --backend xtb --output first_frequency
   pyar inspect first_frequency

This frequency workflow evaluates the supplied geometry; it does not optimize
it. Read ``first_frequency/job_result_water/result.json`` and its frequency
artifacts. ``minimum_confirmed`` requires the implemented stationarity and
vibrational checks. A completed calculation with an unconfirmed minimum is a
scientific warning, not a certificate of success. Significant imaginary modes
call for examining the geometry, convergence and mode displacement.

Observed teaching example
---------------------------

The input above was run with native GFN2-xTB during preparation of this guide.
The single-point energy was -5.065646329768 hartree; the optimized frequency
calculation reported -5.070544406221 hartree and ``minimum_confirmed``. These
are observations from one environment, not universal reference values or new
acceptance tolerances. Record your own backend version and settings; compare
results only under a consistent model and electronic state.

Keep the evidence
-------------------

Retain the input, ``pyar-run.toml``, request/result files and native backend
output. Record the PyAR and xTB versions in your notebook. See :doc:`results`
for units and interpretation, :doc:`analysis` to compare structures and energies,
and :doc:`troubleshooting` if a task fails. Never conclude that an aggregate,
reaction mechanism or transition state is validated from this water example.
