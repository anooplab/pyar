# PyAR
PyAR stands for "Python program for aggregation and reaction"

# Installation

Download the file pyar-master.zip. Unzip it. Go the the folder and ```sudo python3 setup.py install```
This will create python package in the path **/usr/local/lib/python3.6/dist-packages/pyar/**
and will create the command line interface ```pyar-cli``` in **/usr/local/bin**

or

Run the following command in the pyar folder
```sudo -H pip install .```
# Features:
* Automated prediction of unknown reactions between two reactants (A+B)
* Automated prediction of the geometries of aggregates, atomic clusters etc.
* Automated search for reaction for bond forming between two atoms in two different molecules.



## Setting Up Environment

To set up your environment for the tasks, you can create and edit your `.bashrc` or `.bash_profile` file, depending on your system configuration. After that source ~/.bashrc to run those changes.

```bash
# Create an alias for "mndo2020"
alias mndobin="mndo2020"

# Add Gaussian 16 to your PATH
export PATH=$PATH:/home/apps/g16

# Create an alias for GAUSS_EXEDIR
alias GAUSS_EXEDIR="g16"

# Create an alias for MLatom.py
alias mlatom="MLatom.py"

```

```bash
# Install dftd4 executable in this way
conda config --add channels conda-forge
conda install dftd4
conda install -c conda-forge openbabel
alias dftd4bin="dftd4"
```

```bash
#DBCV is not directly accessable via scikit-learn
pip install hdbscan
pip install git+https://github.com/christopherjenness/DBCV.git
```

## Requirements 
* python >= 3.6
* numpy>=1.18.4
* pandas>=1.0.5
* scipy>=1.5.2
* scikit-learn>=0.23.2
* autograd>=1.3

# Interfaced with electronic structure theory programmes
- mlatom_aiqm1
- aimnet2
- Mopac
- Turbomole
- Psi4
- Xtb
- Orca

# Molecule generations 

```pyar-cli -a c.xyz h.xyz -N 8 -as 6 6 --software aiqm1_mlatom  -m 1 2 ```

# Molecular clusters

## XTB
```pyar-cli -s water.xyz water.xyz --software xtb -ss 10  -N 16 -c 0 0 -m 1 1```
## AIMNet2
```pyar-cli -s water.xyz water.xyz --software aimnet_2  -ss 10  -N 16 -c 0 0 -m 1 1```

This will generate a molecules  upto maximum 6 carbon and 6 hydrogens with **mlatom_aiqm1** potential using 8 trial orientations.
Here c.xyz and h.xyz are standard cartesian coordinate files. 
```bash
1
carbon
C  0.0  0.0   0.0
```
```bash
1
hydrogen
H  0.0  0.0   0.0
```
# Reaction

To study the reaction between two reactants A and B using ORCA software interface, with force from 100 to 1000 using N=8 trial orientation, the commandline line argument is,

```pyar-cli -r A.xyz B.xyz -N 8 -gmin 100 -gmax 1000 --ssoftware orca```

A.xyz and B.xyz are the cartesian coordinate files of the reactants

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



# References

1. "A Global Optimizer for Nanoclusters ", Maya Khatun, Rajat Shubhro Majumdar, Anakuthil Anoop <a href="https://www.frontiersin.org/articles/10.3389/fchem.2019.00644/full">Frontiers in Chemistry 2019, 644</a>
1. "A tabu-search based strategy for modeling molecular aggregates and binary reactions" S Nandi, SR McAnanama-Brereton, MP Waller, A Anoop, <a href="https://www.sciencedirect.com/science/article/pii/S2210271X17301627">Computational and Theoretical Chemistry 2017, 1111, 69-81</a>  
