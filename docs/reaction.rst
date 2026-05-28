Reaction Search
===============

Use reaction search when you want PyAR to explore possible bond formation,
bond rearrangement, or close-contact reaction candidates between two input
reactants. The current reaction route uses AFIR-style biased optimisation to
push fragments together and then checks whether unbiased relaxation gives a
new product.

Typical chemistry questions
---------------------------

Reaction search is useful when you want to ask questions such as:

* Can two reactants form a new bonded product from any of many starting
  orientations?
* Which initial orientations lead back to the starting materials, and which
  lead to new structures?
* Where along an AFIR path does the physical backend energy rise most strongly?
* Which geometry should be inspected first as a candidate for later NEB,
  string, dimer, or transition-state work?

Basic command
-------------

Run a reaction search with xTB:

.. code-block:: bash

   pyar-cli react A.xyz B.xyz -N 8 -gmin 100 -gmax 1000 --software xtb
   pyar-cli -r A.xyz B.xyz -N 8 -gmin 100 -gmax 1000 --software xtb

``A.xyz`` and ``B.xyz`` are Cartesian coordinate files for the two reactants.

Backends
--------

The geomeTRIC/TRIC reaction route is used for registered energy-gradient
providers. At present, this route is wired for ``xtb``, ``aimnet_2``,
``orca``, and ``gaussian``. In practice, ``xtb`` and ``aimnet_2`` are the
easier immediately usable options. ``orca`` and ``gaussian`` require installed
executables and should be validated on the target installation.

How the reaction workflow works
-------------------------------

At a high level, PyAR does the following:

1. Generate many starting orientations of the two reactants.
2. Optimise each orientation with a backend energy/force plus an AFIR bias.
3. Detect bonded candidates.
4. Relax bonded candidates again with ``gamma=0.0``.
5. Compare the relaxed structure against the starting reactants using molecular
   identity checks.
6. Save new products and analyse the AFIR trace.

The AFIR-biased optimisation is a product-generation step. The structures in
``candidate_ts/`` are not confirmed transition states.

Reaction outputs
----------------

A reaction run stores restart state and product information under
``reaction/``:

.. code-block:: text

   reaction/
     state.json
     state/
       geometries/
     gamma_.../
       orientation_.../
     products/

Useful files to inspect in a successful orientation job are:

.. code-block:: text

   reaction_trace/
     trace.jsonl
     steps/step_*.xyz
   path_summary.csv
   candidate_ts/
     highest_backend_energy.xyz
     highest_total_energy.xyz
     pre_product_geometry.xyz
     max_bond_change.xyz
     metadata.json
   trace_plots/
     reaction_profile.png

Energy terms in the trace
-------------------------

The trace distinguishes physical energy from the biased optimisation objective:

* ``backend_energy_hartree`` is the backend electronic, ML, or xTB energy
  without AFIR.
* ``afir_energy_hartree`` is the artificial AFIR contribution.
* ``total_energy_hartree`` is the optimisation objective followed by
  geomeTRIC.
* ``backend_relative_kcalmol`` reports the backend-energy change relative to
  the first recorded trace frame.

The file ``candidate_ts/highest_backend_energy.xyz`` is usually the first
structure to inspect for later path-refinement work because it is selected from
the physical backend energy, not the AFIR-biased total energy.

Trace analysis
--------------

After a successful geomeTRIC-backed AFIR reaction job, analyse one orientation
job directory with:

.. code-block:: bash

   pyar-reaction-trace reaction/gamma_0100/orientation_xxxxx --plot

or, from inside the orientation directory:

.. code-block:: bash

   pyar-reaction-trace . --plot

This writes ``path_summary.csv``, candidate geometries in ``candidate_ts/``,
and optional plots in ``trace_plots/``.

Important caution
-----------------

Files in ``candidate_ts/`` are candidate geometries for later path refinement.
They are not confirmed transition states. Confirmation requires a separate
NEB, string, dimer, or TS optimisation and an IRC or downhill connectivity
check.
