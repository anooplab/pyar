PyAR Documentation
==================

PyAR is a command-line package for chemistry-focused structure search:
aggregation of noncovalent complexes and clusters, AFIR-style reaction
searches, microsolvation, and bond scans.

It is most useful when you want to:

* explore low-energy cluster geometries for a set of fragments
* search for products or close-contact reaction candidates between two
  reactants
* place solvent molecules around a solute and follow growth cycles
* scan a bond distance between two fragments as a simple reaction-coordinate
  probe

.. rubric:: Cite this work

If you use PyAR in a project, start with :doc:`publications` and cite the
workflow paper that best matches your chemistry problem. The sections below
show the most common mappings.

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

If you are a student or researcher, start with :doc:`quickstart` and
:doc:`workflows` to see which command matches your chemistry problem. If you
want examples of published uses, see :doc:`publications`.

Which paper matches my problem?
--------------------------------

.. list-table::
   :header-rows: 1
   :widths: 30 22 48

   * - Chemistry problem
     - Workflow
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

* :doc:`quickstart`
* :doc:`workflows`
* :doc:`usage`
* :doc:`publications`
* :doc:`installation`
* :doc:`energy_table`
* :doc:`reaction_optimization`
* :doc:`molecule`
* :doc:`orientation_sampling`
* :doc:`xtb`

Developer and API details are still available later in the manual if you need
them:

* :doc:`api`
* :doc:`reference`
* :doc:`architecture`

How to cite PyAR
----------------

There is no single citation that covers every PyAR use case. In general:

* cite the workflow paper that best matches what you used
* if you used more than one workflow, cite the relevant workflow papers
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
   :caption: Contents

   quickstart
   workflows
   usage
   publications
   installation
   energy_table
   molecule
   orientation_sampling
   xtb
   reaction_optimization
   api
   reference
   architecture
