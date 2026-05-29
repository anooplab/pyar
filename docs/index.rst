PyAR Documentation
==================

PyAR is a chemistry-focused structure-search package. It helps generate,
optimise, select, and inspect plausible molecular structures before committing
to higher-level electronic-structure calculations.

PyAR is most useful when you want to:

* explore low-energy cluster geometries for a set of fragments
* search for products or close-contact reaction candidates between two
  reactants
* place solvent molecules or ligands around a central structure
* scan a bond distance between two fragments as a simple reaction-coordinate
  probe

For most chemistry users, start with the task that matches your problem:

* :doc:`first_run` for a complete first calculation
* :doc:`aggregation` for clusters, aggregates, and noncovalent complexes
* :doc:`reaction` for AFIR-style product and close-contact searches
* :doc:`solvation` for microsolvation, ligand addition, and growth around a core
* :doc:`bond_scan` for a simple distance-coordinate probe

.. rubric:: Cite this work

If you use PyAR in a project, start with :doc:`publications` and cite the
paper that best matches your chemistry problem. The sections below show the
most common mappings.

Why chemists use PyAR
---------------------

PyAR is useful when you want an automated way to generate and screen plausible
structures before committing to higher-level electronic-structure work. In
practice that means:

* getting candidate low-energy geometries for clusters and complexes
* finding plausible reaction products or prebiotic intermediates
* exploring ligand addition or catalyst formation without hand-building every
  guess structure
* testing solvation and coordination growth around a central core
* following a simple bond coordinate when you want a quick structural probe

If you are a student or researcher, start with :doc:`quickstart`, then move to
one of the chemistry task pages above. If you want examples of published uses,
see :doc:`publications`.

Which paper matches my problem?
-------------------------------

.. list-table::
   :header-rows: 1
   :widths: 30 22 48

   * - Chemistry problem
     - Command/task
     - Example publications
   * - Noncovalent cluster growth and aggregation
     - ``aggregate``
     - [Nandi2017]_, [Khatun2019]_, [Sherpa2026]_
   * - Prebiotic reaction discovery and bond rearrangement
     - ``react``
     - [Nandi2018]_, [Panda2024]_
   * - Microsolvation or ligand addition around a central core
     - ``solvate``
     - Use the same build-up logic for solvent shells, coordination
       complexes, and organometallic assembly
   * - Catalyst formation and sequential ligand addition
     - ``react``
     - [Roy2022]_
   * - Chemical-space exploration of noncovalent clusters
     - ``aggregate``
     - [Giri2025]_

Start here
----------

* :doc:`first_run`
* :doc:`quickstart`
* :doc:`aggregation`
* :doc:`reaction`
* :doc:`solvation`
* :doc:`bond_scan`
* :doc:`usage`
* :doc:`publications`
* :doc:`installation`

Tools and references
--------------------

* :doc:`energy_table`
* :doc:`molecule`
* :doc:`orientation_sampling`
* :doc:`backends`

Backend guides
--------------

* :doc:`backends`
* :doc:`xtb`
* :doc:`aimnet2_backend`
* :doc:`orca_backend`
* :doc:`gaussian_backend`
* :doc:`mopac_backend`
* :doc:`obabel_backend`
* :doc:`psi4_backend`

Developer and API details
-------------------------

These pages are useful when you are modifying PyAR or trying to understand its
internal design:

* :doc:`workflows`
* :doc:`reaction_optimization`
* :doc:`api`
* :doc:`reference`
* :doc:`generated_api`
* :doc:`architecture`

How to cite PyAR
----------------

There is no single citation that covers every PyAR use case. In general:

* cite the paper that best matches what you used
* if you used more than one task, cite the relevant papers
* use :doc:`publications` as the short list of chemistry-facing examples

For a general citation, the two original papers that introduced the main build
up and reaction-search ideas are [Nandi2017]_ and [Khatun2019]_.

.. note::

   If your work used AFIR-style reaction search, cluster growth, solvation, or
   another specific application, also cite the corresponding paper from
   :doc:`publications`.

The fastest way to verify a local install is:

.. code-block:: bash

   pyar-cli --help

.. toctree::
   :maxdepth: 2
   :caption: Getting started

   first_run
   quickstart
   installation
   usage
   publications

.. toctree::
   :maxdepth: 2
   :caption: Chemistry tasks

   aggregation
   reaction
   solvation
   bond_scan

.. toctree::
   :maxdepth: 2
   :caption: Tools and reference

   energy_table
   molecule
   orientation_sampling

.. toctree::
   :maxdepth: 2
   :caption: Backend guides

   backends
   xtb
   aimnet2_backend
   orca_backend
   gaussian_backend
   mopac_backend
   obabel_backend
   psi4_backend

.. toctree::
   :maxdepth: 2
   :caption: Developer documentation

   workflows
   reaction_optimization
   api
   reference
   generated_api
   architecture
