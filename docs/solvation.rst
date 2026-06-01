Solvation and Growth Around a Core
==================================

Use solvation when you want PyAR to add one or more copies of a fragment around
a central structure. Microsolvation is the classic use case, but the same idea
also applies to ligand addition, coordination growth, and building local
cluster environments around a molecular or metal-containing core.

Disconnected covalent graphs are expected in these noncovalent growth
problems, so the workflow does not drop candidates just because they are not
single covalent components.

Typical chemistry questions
---------------------------

Solvation is useful when you want to ask questions such as:

* How do solvent molecules arrange around a solute?
* What are plausible first-shell or second-shell microsolvation structures?
* How can ligands attach around a metal centre or molecular core?
* Which solvent, ligand, or additive arrangements should be refined with a
  higher-level backend?

Basic command
-------------

Run a solvation search with xTB:

.. code-block:: bash

   pyar-cli solvate solute.xyz solvent.xyz --software xtb -ss 10 -N 16
   pyar-cli -s solute.xyz solvent.xyz --software xtb -ss 10 -N 16

The first XYZ file is the core or solute. The second XYZ file is the fragment
that will be added around the core.

How the solvation workflow works
--------------------------------

At a high level, PyAR does the following:

1. Start from the solute or central core.
2. Generate trial orientations of the added fragment.
3. Optimise and select diverse low-energy structures.
4. Use selected structures as seeds for the next growth cycle.
5. Repeat until the requested number of added fragments is reached.

This means the command can be used for more than water solvation. It can also
model ligand addition, ion coordination, and local growth around molecular
clusters.

Outputs and restart state
-------------------------

Solvation restart state is stored as readable JSON:

.. code-block:: text

   solvation/
     state.json
     state/
       geometries/
   aggregate_002/
   aggregate_003/

``state.json`` records the input seed, added fragment, calculation settings,
next cycle, completed cycles, and current seeds. Re-running an interrupted
solvation with the same request resumes from the last completed cycle and
reuses the stored seed geometries.

Useful files to inspect:

* ``solvation/state.json`` for restart and cycle progress
* ``solvation/state/geometries/`` for saved seed structures
* selected structures from the final cycle

Common follow-up steps
----------------------

After a solvation or ligand-growth run, common follow-up steps are:

* inspect selected structures visually
* generate an energy table for the selected structures
* refine a smaller selected set with a higher-level backend
* compare different solvents, ligands, or ion-pair environments
