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

Download the file pyar-master.zip. Unzip it. Go the the folder and ```sudo python3 setup.py install```
This will create python package in the path **/usr/local/lib/python3.6/dist-packages/pyar/**
and will create the command line interface ```pyar-cli``` in **/usr/local/bin**

# Features:
* Automated prediction of reactions between two reactants (A+B)
* Automated prediction of the geometries of aggregates, atomic clusters etc.
* Generate random orientations between two molecules.



# Interfaced with electronic structure theory programmes
- Turbomole
- Mopac
- Xtb
- Orca
- Psi4

# Command line usage
```
usage: pyar-cli options files

type pyar-cli -h for opttions
```


# Examples

## Reaction

To study the reaction between two reactants A and B using ORCA software interface, with force from 100 to 1000 using N=8 trial orientation, the commandline line argument is,  

```pyar-cli -r A.xyz B.xyz -N 8 -gmin 100 -gmax 1000 --ssoftware orca```

A.xyz and B.xyz are the cartesian coordinate files of the reactants.


## Aggregation and clustering

```pyar-cli -a A.xyz A.xyz -N 8 -as 10 --software xtb```

This will generate a molecular aggregate/cluster of size upto 10 with **XTB** package using 8 trial orientations.

# References

1. "A Global Optimizer for Nanoclusters ", Maya Khatun, Rajat Shubhro Majumdar, Anakuthil Anoop <a href="https://www.frontiersin.org/articles/10.3389/fchem.2019.00644/full">Frontiers in Chemistry 2019, 644</a>
1. "A tabu-search based strategy for modeling molecular aggregates and binary reactions" S Nandi, SR McAnanama-Brereton, MP Waller, A Anoop, <a href="https://www.sciencedirect.com/science/article/pii/S2210271X17301627">Computational and Theoretical Chemistry 2017, 1111, 69-81</a>  

# Credits

The idea and the initial work started with Dr. Mark P. Waller at University of Muenster in 2011. Mark is currently at <a href="http://pending.ai">pending.ai</a>.
Later Dr. Surajit Nandi who as PhD student in AnoopLab made immense contributions.  Another major help came from a Masters project student, Deebankur Bhattacharyya.

# License

GNU GPL V3