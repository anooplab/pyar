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

Stereochemistry is tracked from the 3D coordinates. Each record in
``summary.csv`` and the run state carries a ``stereo_signature`` (canonical
isomeric SMILES), and the state's ``stereochemistry`` entry lists any
unspecified stereo elements in the input. When the input specifies every stereo
element, generated structures with a different configuration are discarded,
and a backend refinement that changes the configuration is excluded with
status ``stereo_changed``. When stereo elements are unspecified, PyAR warns and
keeps the stereoisomers that ETKDG produces, labelled by their signatures; give
explicit stereochemistry (for example ``C[C@H](O)CC``) to search one isomer.

Inspect ``conformers/results`` and its manifest for selected structures, with
state and generation/refinement records in the run directory. Compare the
selected torsions and distinct basins rather than judging diversity solely by
count. Keep seed/thread/settings provenance; different sampling budgets can
change the ensemble. A selected conformer is not automatically a vibrationally
verified minimum. Run a consistent frequency protocol when that matters.


See :doc:`/reference/cli` for current options and :doc:`results` for interpretation.
