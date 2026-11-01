PyAR: molecular structures and reaction exploration
===================================================

Find candidate structures, calculate their properties, and examine the evidence
behind a chemical interpretation. PyAR connects structure sampling and selection
to electronic-structure programs and machine-learning potentials through one
command: ``pyar``.

Start with :doc:`user/installation`, then follow :doc:`user/first_calculation`.
The first tutorial calculates, optimizes and checks a water molecule with xTB.
It explains which output to read and what each calculation establishes.

.. toctree::
   :maxdepth: 2
   :caption: Learn PyAR

   user/installation
   user/first_calculation
   user/workflow_chooser
   user/input_and_settings
   user/results
   user/troubleshooting

.. toctree::
   :maxdepth: 2
   :caption: Chemistry workflows

   user/singlepoint
   user/conformer
   user/aggregate
   user/grow
   user/microsolvate
   user/react
   user/scan_bond
   user/reaction_paths
   user/analysis
   user/configuration
   user/restart

.. toctree::
   :maxdepth: 2
   :caption: Reference

   reference/cli
   reference/tools
   reference/backends
   reference/formats
   reference/citations
   reference/upgrade

.. toctree::
   :maxdepth: 2
   :caption: For developers

   developer/index

PyAR generates and screens hypotheses about molecular structures. A low energy,
a successful optimizer, a minimum, a transition state, and a connected reaction
path establish different things. Each workflow guide describes the relevant
acceptance checks and the evidence still needed for your research question.
