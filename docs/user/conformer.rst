Search molecular conformers
===========================

RDKit generates molecular conformers and applies the selected force field and
diversity controls. Optional backend refinement is a separate scientific model;
force-field rankings and electronic energies must be interpreted accordingly.

.. code-block:: bash

   pyar conformer CCO --num-conformers 20 --top-n 3 --check
   pyar conformer CCO --num-conformers 20 --top-n 3
   pyar conformer CCO --num-conformers 20 --top-n 3 --backend xtb --check

Install the ``conformer`` extra. The input may be SMILES or a supported molecular
file. XYZ input needs the documented OpenBabel conversion route; it does not
encode all bonding/stereochemical information as reliably as a molecular format.
Check inferred protonation, formal charge and stereochemistry before generation.

Every stereo element must be specified. If the input leaves a stereocentre or
stereogenic double bond undefined (``CC(O)CC`` rather than ``C[C@H](O)CC``), the
run stops before generation and names the undefined elements, because each
structure would otherwise receive a random configuration and different
stereoisomers would be mixed into one conformer set. Run each stereoisomer as
its own search.

Each record in ``summary.csv`` and the run state carries a ``stereo_signature``
(canonical isomeric SMILES) perceived from its 3D coordinates. Generated
structures whose configuration differs from the input are discarded, and a
backend refinement that changes the configuration is excluded with status
``stereo_changed``. If input stereo perception fails, both ``--check`` and
execution stop. Generated structures with unverified configuration are discarded;
backend refinements with failed perception are excluded as ``stereo_unverified``.
An unavailable signature does not establish agreement with the input.

Inspect ``conformers/results`` and its manifest for selected structures, with
state and generation/refinement records in the run directory. Compare the
selected torsions and distinct basins rather than judging diversity solely by
count. Keep seed/thread/settings provenance; different sampling budgets can
change the ensemble. A selected conformer is not automatically a vibrationally
verified minimum. Run a consistent frequency protocol when that matters.


See :doc:`/reference/cli` for current options and :doc:`results` for interpretation.
