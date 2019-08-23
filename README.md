# PyAR
PyAR stands for "Python program for aggregation and reaction"

# Contents

1.  [Features](#features)
2.  [Interfaced with electronic structure theory programs:](#interface)
3.  [Command Line Interface](#cli)
4.  [Example usages](#examples)

# Features:
* Automated prediction of reactions between two reactants (A+B)
* Automated prediction of the geometries of aggregates, atomic clusters etc.
* Generate random orientations between two molecules.



# Interfaced with electronic structure theory programs:<a name="interface"></a>
- Turbomole
- Mopac
- Xtb
- Orca
- Psi4

# Command line usage<a name="cli"></a>
```
usage: PyAR [-h] [-N HM_ORIENTATIONS] [--site SITE SITE]
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
```

# Examples

## Reaction

To study the reaction between two reactants A and B using ORCA software interface, with force from 100 to 1999 using N=8 trial orientation, the commandline line argument is,  
```python3 <path>/pyar/cli -r A.xyz B.xyz -N 8 -gmin 100 -gmax 1000 --ssoftware orca```
