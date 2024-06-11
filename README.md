# PyMLGen
This is an integration of pyar aggregator module and MLAtom package, where the Tabu heuristic of aggregator module run with the interface of AIQM1 potential.

## Setting Up Environment

To set up your environment for the tasks, you can create and edit your `.bashrc` or `.bash_profile` file, depending on your system configuration. After that source ~/.bashrc to run those changes.

```bash
# Create an alias for "mndo2020"
alias mndobin="mndo2020"

# Add MLatom_v2.3.3 to your PATH
export PATH=$PATH:/home/20cy91r19/apps/MLatom_v2.3.3

# Add Gaussian 16 to your PATH
export PATH=$PATH:/home/apps/g16

# Create an alias for GAUSS_EXEDIR
alias GAUSS_EXEDIR="g16"

# Create an alias for MLatom.py test
alias mlatom-test="MLatom.py"

```

```bash
# Install dftd4 executable in this way
conda config --add channels conda-forge
conda install dftd4
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
- Xtb
- Orca

# Molecule generations 

```pymlgen-cli -a c.xyz h.xyz -N 8 -as 6 6 --software mlatom_aiqm1  -m 1 2 ```

# Molecular clusters

## XTB
```pymlgen-cli -s water.xyz water.xyz --software xtb -ss 10  -N 16 -c 0 0 -m 1 1```
## AIMNet2
```pymlgen-cli -s water.xyz water.xyz --software aimnet-2 --model -ss 10  -N 16 -c 0 0 -m 1 1```

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


# References

1. "A Global Optimizer for Nanoclusters ", Maya Khatun, Rajat Shubhro Majumdar, Anakuthil Anoop <a href="https://www.frontiersin.org/articles/10.3389/fchem.2019.00644/full">Frontiers in Chemistry 2019, 644</a>
1. "A tabu-search based strategy for modeling molecular aggregates and binary reactions" S Nandi, SR McAnanama-Brereton, MP Waller, A Anoop, <a href="https://www.sciencedirect.com/science/article/pii/S2210271X17301627">Computational and Theoretical Chemistry 2017, 1111, 69-81</a>  
