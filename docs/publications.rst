Selected Publications
=====================

PyAR has been used in a range of chemistry problems beyond pure algorithm
development. This page highlights a few representative papers so a student or
researcher can see what kinds of questions the workflows have already been
used to study. The list is selective rather than exhaustive.

Aggregation And Clusters
------------------------

PyAR's aggregation and cluster-building workflows are useful whenever the
problem is to grow a low-energy structure from smaller fragments.

* [Nandi2017]_ Tabu-guided build-up for molecular aggregates and binary
  reactions.
* [Khatun2019]_ Cluster optimization and global search for nanoclusters.
* [Sherpa2026]_ Noncovalent-cluster application of the build-up workflow.

Reaction Discovery And Prebiotic Chemistry
------------------------------------------

PyAR's AFIR-based reaction workflow is suitable for exploring products,
reaction pathways, and chemically plausible intermediates.

* [Nandi2018]_ Automated reaction search for prebiotic HCN tetramerization.
* [Panda2024]_ Prebiotic reaction channels in an extraterrestrial atmosphere.

Catalyst Formation And Sequential Ligand Addition
--------------------------------------------------

The same reaction-search machinery can also follow stepwise ligand addition
or catalyst assembly.

* [Roy2022]_ Active catalyst formation from dinuclear palladium acetate by
  sequential ligand addition.

Exploration Of Chemical Space
-----------------------------

PyAR can also be used as a general build-up engine for exploring clusters and
other parts of chemical space.

* [Giri2025]_ Chemical-space exploration of noncovalent molecular clusters.

Next Steps
----------

If you want to use PyAR for a similar chemistry problem:

* start with :doc:`quickstart` to choose the right command
* use :doc:`aggregate`, :doc:`react`, and :doc:`solvate` for the main tasks
* use :doc:`workflows` to see what each workflow writes to disk
* use :doc:`usage` for helper commands and analysis tools

References
----------

.. [Nandi2017] Surajit Nandi, S. R. McAnanama-Brereton, M. P. Waller, and
   Anakuthil Anoop, *A tabu-search based strategy for modeling molecular
   aggregates and binary reactions*, *Computational and Theoretical Chemistry*
   1111, 69-81 (2017). https://doi.org/10.1016/j.comptc.2017.03.040
.. [Khatun2019] Maya Khatun, Rajat Shubhro Majumdar, and Anakuthil Anoop,
   *A Global Optimizer for Nanoclusters*, *Frontiers in Chemistry* 7:644
   (2019). https://doi.org/10.3389/fchem.2019.00644
.. [Sherpa2026] Lazumla T. Sherpa, Maya Khatun, Sunanda Panda, and
   Anakuthil Anoop, *Cooperative many-body interactions and spectroscopic
   signatures in (HCN)_{n} and (HNC)_{n} clusters up to n = 15*, *Physical
   Chemistry Chemical Physics* 28, 2227-2238 (2026).
   https://doi.org/10.1039/D5CP03273C
.. [Nandi2018] Surajit Nandi, Debankur Bhattacharyya, and Anakuthil Anoop,
   *Prebiotic Chemistry of HCN Tetramerization by Automated Reaction Search*,
   *Chemistry - A European Journal* 24(19), 4885-4894 (2018).
   https://doi.org/10.1002/chem.201705492
.. [Panda2024] Sunanda Panda and Anakuthil Anoop, *Potential Prebiotic
   Pathways in Extraterrestrial Atmosphere: A Computational Exploration of
   HCN and NH_{3} Reactions*, *ACS Earth and Space Chemistry* (2024).
   https://doi.org/10.1021/acsearthspacechem.3c00321
.. [Roy2022] Saikat Roy and Anakuthil Anoop, *Insights into the Active
   Catalyst Formation from Dinuclear Palladium Acetate in Pd-Catalyzed
   Coupling Reactions: A DFT Study*, *The Journal of Physical Chemistry A*
   126(46), 8562-8576 (2022). https://doi.org/10.1021/acs.jpca.2c03762
.. [Giri2025] Sandip Giri and Anakuthil Anoop, *Exploring the Chemical Space
   of Noncovalent Molecular Clusters Using Automated Cluster Building
   Algorithm and Neural Network Potential*, *Journal of Computational
   Chemistry* 46(32), e70287 (2025). https://doi.org/10.1002/jcc.70287
