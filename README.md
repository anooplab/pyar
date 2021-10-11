# PyAR
PyAR stands for "Python program for aggregation and reaction"

# Contents

1.  [Installation](#installation)
1.  [Features](#features)
1.  [Interfaced with electronic structure theory programs:](#interfaced-with-electronic-structure-theory-programmes)
1.  [Command Line Interface](#command-line-usage)
1.  [Example usages](#examples)
1.  [References](#references)
1.  [Credits](#credits)
1.  [License](#license)

# Installation

## Download or clone

Download the file pyar-master.zip. Unzip it. 

or clone the repository

```git clone https://github.com/anooplab/pyar.git```

## Install requirements

```sudo python3 -m pip install -r requirements.txt```

## Install PyAR

Go the the folder and ```python3 setup.py install```
This will create python package in the path **$HOME/.local/lib/python3.xx/dist-packages/pyar/**
and will create the command line interface ```pyar-cli``` in **$HOME/.local/bin**

or

Run the following command in the pyar folder
```python3 -m Hpip install .```

# Features:
* Automated prediction of unknown reactions between two reactants (A+B)
* Automated prediction of the geometries of aggregates, atomic clusters etc.
* Automated search for reaction for bond forming between two atoms in two different molecules.

## Requirements
* python >= 3.6
* numpy>=1.18.4
* pandas>=1.0.5
* scipy>=1.5.2
* scikit-learn>=0.23.2
* autograd>=1.3

# Interfaced with electronic structure theory programmes
- Turbomole
- Mopac
- Xtb
- Orca
- Psi4

# Command line usage

The easy usage is by command line. This can also be used as library.
The detailed help can be obtained by typing ```command -h```

## pyar-cli
The main program can be used as below:

```
pyar-cli options files
```

There are other scripts for a few automation tasks.

## pyar-tabu
pyar-tabu can be used for
* for making different orientations of two molecules.
* Making a composite molecule containing a _seed_ molecule and __N__ number of monomer molecules.
* Orient two molecules such that _i_'th atom of one molecule and _j_ 'th atom of second molecule have shortest distance 
between them

## pyar-clustering
* for a clustering analysis of n input molecules to find unique molecules. 

## pyar-optimiser
* for the bulk optimisation of several molecules

# Examples

## Reaction

To study the reaction between two reactants A and B using ORCA software interface, with force from 100 to 1000 using N=8 trial orientation, the commandline line argument is,  

```pyar-cli -r A.xyz B.xyz -N 8 -gmin 100 -gmax 1000 --ssoftware orca```

A.xyz and B.xyz are the cartesian coordinate files of the reactants.


## Aggregation and clustering

```pyar-cli -a A.xyz -N 8 -as 10 --software xtb```

This will generate a molecular aggregate/cluster of size upto 10 with **XTB** package using 8 trial orientations.
The A.xyz is a standard cartesian coordinate file. For example the .xyz file for a carbon atom is the following:
```
1
Comment
C  0.0  0.0   0.0
```

# References

1. "A Global Optimizer for Nanoclusters ", Maya Khatun, Rajat Shubhro Majumdar, Anakuthil Anoop <a href="https://www.frontiersin.org/articles/10.3389/fchem.2019.00644/full">Frontiers in Chemistry 2019, 644</a>
1. "A tabu-search based strategy for modeling molecular aggregates and binary reactions" S Nandi, SR McAnanama-Brereton, MP Waller, A Anoop, <a href="https://www.sciencedirect.com/science/article/pii/S2210271X17301627">Computational and Theoretical Chemistry 2017, 1111, 69-81</a>  

# Credits

The idea and the initial work started with Dr. Mark P. Waller at University of Muenster in 2011. Mark is currently at <a href="http://pending.ai">pending.ai</a>.
Later Dr. Surajit Nandi who as PhD student in AnoopLab made immense contributions.  Another major help came from a Masters project student, Deebankur Bhattacharyya.

# License

GNU GPL V3
