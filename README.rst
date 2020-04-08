PyAR
====

PyAR stands for "Python program for aggregation and reaction"

Contents
========

1. `Installation <#install>`__
2. `Features <#features>`__
3. `Interfaced with electronic structure theory
   programs: <#interface>`__
4. `Command Line Interface <#cli>`__
5. `Example usages <#examples>`__
6. `References <#references>`__
7. `Credits <#credits>`__
8. `License <#license>`__

Installation
============

Download the file pyar-master.zip. Unzip it. Go the the folder and
``sudo python3 setup.py install`` This will create python package in the
path **/usr/local/lib/python3.6/dist-packages/pyar/** and will create
the command line interface ``pyar-cli`` in **/usr/local/bin**

Features:
=========

-  Automated prediction of reactions between two reactants (A+B)
-  Automated prediction of the geometries of aggregates, atomic clusters
   etc.
-  Generate random orientations between two molecules.

Interfaced with electronic structure theory programs:
=====================================================

-  Turbomole
-  Mopac
-  Xtb
-  Orca
-  Psi4

Command line usage
==================

::

    usage: pyar-cli [-h] [-N HM_ORIENTATIONS] [--site SITE SITE]
                (-r | -a | -ba | -o | -mb MAKEBOND MAKEBOND) [-as AGGREGATE_SIZE]
                [-bs1 BAGGREGATOR_SIZE1] [-bs2 BAGGREGATOR_SIZE2]
                [-mns MAXIMUM_NUMBER_OF_SEEDS] [-gmin GMIN] [-gmax GMAX]
                [-gamma GAMMA] [-c CHARGE] [-m MULTIPLICITY] [--scftype {rhf,uhf}]
                --software {turbomole,obabel,mopac,xtb,xtb_turbo,orca,psi4}
                [-v {0,1,2,3,4}]
                files [files ...]

    is a program to predict aggregation, reaction, clustering. There are also
    modules for strochastic generation of orientations of two more molecules and
    atoms

    positional arguments:
      files                 input coordinate files

    optional arguments:
      -h, --help            show this help message and exit
      -N HM_ORIENTATIONS    how many orientations to be used
      --site SITE SITE      atom for site specific aggregation/solvation
      -r, --react           Run a reactor calculation
      -a, --aggregate       Run a aggregator calculation
      -ba, --baggregate     Run a baggregator calculation
      -o, --optimize        Optimize the molecules
      -mb MAKEBOND MAKEBOND, --makebond MAKEBOND MAKEBOND
                            make bonds between the given atoms of twofragments
      -v {0,1,2,3,4}, --verbosity {0,1,2,3,4}
                            increase output verbosity (0=Debug; 1=Info; 2:
                            Warning; 3: Error; 4: Critical)

    aggregator:
      Aggregator specific option

      -as AGGREGATE_SIZE, --aggregate-size AGGREGATE_SIZE
                            number of monomers in aggregate

    baggregator:
      Baggregator specific option

      -bs1 BAGGREGATOR_SIZE1, --baggregator-size1 BAGGREGATOR_SIZE1
                            size of atom1
      -bs2 BAGGREGATOR_SIZE2, --baggregator-size2 BAGGREGATOR_SIZE2
                            size of atom2
      -mns MAXIMUM_NUMBER_OF_SEEDS, --maximum-number-of-seeds MAXIMUM_NUMBER_OF_SEEDS
                            maximum number of seeds

    reactor:
      Reactor specific option

      -gmin GMIN            minimum value of gamma
      -gmax GMAX            maximum value of gamma

    optimizer:
      Optimizer specific option

      -gamma GAMMA          value of gamma

    chemistry:
      Options related to model chemistry

      -c CHARGE, --charge CHARGE
                            Charge of the system
      -m MULTIPLICITY, --multiplicity MULTIPLICITY
                            Multiplicity of the system
      --scftype {rhf,uhf}   specify rhf or uhf (defulat=rhf)
      --software {turbomole,obabel,mopac,xtb,xtb_turbo,orca,psi4}
                            Software

Examples
========

Reaction
--------

To study the reaction between two reactants A and B using ORCA software
interface, with force from 100 to 1000 using N=8 trial orientation, the
commandline line argument is,

``pyar-cli -r A.xyz B.xyz -N 8 -gmin 100 -gmax 1000 --ssoftware orca``

A.xyz and B.xyz are the cartesian coordinate files of the reactants.

Aggregation and clustering
--------------------------

``pyar-cli -a A.xyz A.xyz -N 8 -as 10 --software xtb``

This will generate a molecular aggregate/cluster of size upto 10 with
**XTB** package using 8 trial orientations.

References
==========

1. "A Global Optimizer for Nanoclusters ", Maya Khatun, Rajat Shubhro
   Majumdar, Anakuthil Anoop Frontiers in Chemistry 2019, 644
2. "A tabu-search based strategy for modeling molecular aggregates and
   binary reactions" S Nandi, SR McAnanama-Brereton, MP Waller, A Anoop,
   Computational and Theoretical Chemistry 2017, 1111, 69-81

Credits
=======

The idea and the initial work started with Dr. Mark P. Waller at
University of Muenster in 2011. Mark is currently at pending.ai. Later
Dr. Surajit Nandi who as PhD student in AnoopLab made immense
contributions. Another major help came from a Masters project student,
Deebankur Bhattacharyya.

License
=======

GNU GPL V3
