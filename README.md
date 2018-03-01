# PyAR
PyAR stands for "Python program for aggregation and reaction"

## Features:
* Automated prediction of reactions between two reactants (A+B)
* Automated prediction of the geometries of aggregates, atomic clusters etc.

## Interfaced with electronic structure theory programs:
- Turbomole
- Mopac
- Xtb
- Orca

## Command line usage
usage: PyAR [-h] [-N HM_ORIENTATIONS] [-cite CITE] (-r | -a | -o)
            [-as AGGREGATE_SIZE] [-gmin GMIN] [-gmax GMAX] [-gamma GAMMA]
            [-c CHARGE] [-m MULTIPLICITY] [--scftype {rhf,uhf}] --software
            {turbomole,OBabel,mopac,xtb,xtb_turbo,orca}
            files [files ...]

is a program to predict aggregation and reaction

positional arguments:
  files                 input coordinate files

optional arguments:
  -h, --help            show this help message and exit
  -N HM_ORIENTATIONS    how many orientations to be used
  -cite CITE, --cite CITE
                        atom for site specific aggregation/solvation
  -r, --react           Run a reactor calculation
  -a, --aggregate       Run a aggregator calculation
  -o, --optimize        Optimize the molecules

aggregator:
  Aggregator specific option

  -as AGGREGATE_SIZE, --aggregate-size AGGREGATE_SIZE
                        number of monomers in aggregate

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
  --software {turbomole,OBabel,mopac,xtb,xtb_turbo,orca}
                        Software

