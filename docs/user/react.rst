Explore candidate reaction products
===================================

Reaction exploration uses an explicit bias and subsequent release/selection
logic to search for candidate products. Bias choices guide a search; they do not
establish a barrier or a physically sampled time evolution.

.. code-block:: bash

   pyar react reactant_a.xyz reactant_b.xyz --backend xtb --bias-max 1000 --orientations 4 --check
   pyar react reactant_a.xyz reactant_b.xyz --backend xtb --bias-max 1000 --orientations 4

Supply chemically appropriate reactants, electronic states and bias ceiling.
The value above demonstrates syntax, not a universally suitable research value.
Use the workflow help and methodology to choose your own settings. A preflight
check does not establish that the reaction will be found.

Inspect the reaction run's discovery/product records and ``pyar trace`` output.
Separate physical energy from biased objective; inspect unbiased release and
identity decisions before calling a candidate a product. Preserve failed jobs
and budget accounting. A discovered rearrangement is not a verified mechanism.
Use consistent endpoint optimization, a validated saddle and two-direction IRC
when making a connection claim; see :doc:`reaction_paths`.


See :doc:`/reference/cli` for current options and :doc:`results` for interpretation.
