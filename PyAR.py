import argparse
import copy
import sys
import time
import aggregator
import reactor
from optimiser import optimise
from Molecule import Molecule


def argument_parse():
    """ Parse command line arguments"""
    parser = argparse.ArgumentParser(prog='PyAR', description='is a \
             program to predict aggregation and reaction')
    parser.add_argument("input_files", metavar='files',
                        type=str, nargs='+', help='input coordinate files')
    parser.add_argument("-N", dest='hm_orientations',
                        help='how many orientations to be used')
    run_type_group = parser.add_mutually_exclusive_group(required=True)
    run_type_group.add_argument("-r", "--react",
                                help="Run a reactor calculation", action='store_true')
    run_type_group.add_argument("-a", "--aggregate",
                                help="Run a aggregator calculation", action='store_true')
    run_type_group.add_argument("-o", "--optimize",
                                help="Optimize the molecules", action='store_true')
    aggregator_group = parser.add_argument_group('aggregator', 'Aggregator specific option')
    aggregator_group.add_argument('-as', '--aggregate-size', type=int,
                                  help='number of monomers in aggregate')
    reactor_group = parser.add_argument_group('reactor', 'Reactor specific option')
    reactor_group.add_argument('-gmin', type=float,
                               help='minimum value of gamme')
    reactor_group.add_argument('-gmax', type=float,
                               help='maximum value of gamme')

    chemistry_group = parser.add_argument_group('chemistry', 'Options related\
                                           to model chemistry')
    chemistry_group.add_argument("-c", "--charge", type=int, default=0,
                                 help="Charge of the system")
    chemistry_group.add_argument("-m", "--multiplicity", type=int,
                                 default=1, help="Multiplicity of the system")
    chemistry_group.add_argument("--scftype", type=str, choices=['rhf', 'uhf'],
                                 default='rhf', help="specify rhf or uhf (defulat=rhf)")
    chemistry_group.add_argument("--software", type=str, choices=['turbomole', 
                                 'obabel', 'mopac', 'xtb', 'xtb_turbo', 'orca'], 
                                 required=True, help="Software")

    return parser.parse_args()


def setup_molecules(input_files):

    molecules = []
    for each_file in input_files:
        try:
            mol = Molecule.from_xyz(each_file)
            print(each_file)
            molecules.append(mol)
        except IOError:
            print("File {} does not exist".format(each_file))
            sys.exit()
    print("I've parsed these molecules as input: {}". format([i.name for i in molecules]))
    return molecules


def main():
    args = argument_parse()
    method_args = {
        'charge': args.charge,
        'multiplicity': args.multiplicity,
        'scftype': args.scftype,
        'software': args.software
    }
    N = args.hm_orientations
    size_of_aggregate = args.aggregate_size
    minimum_gamma = args.gmin
    maximum_gamma = args.gmax
    input_files = args.input_files

    input_molecules = setup_molecules(input_files)

    if args.aggregate:
        if size_of_aggregate is None:
            print('For an Aggregation run '
                  'specify the aggregate size '
                  '(number of monomers to be added) '
                  'using the argument\n -as <integer>')
            sys.exit()
        if N is None:
            print("For aggregation, specify how many orientations"
                  "are to be used, by the argument\n"
                  "-N <number of orientations>")
        if len(input_molecules) == 1:
            input_molecules.append(copy.copy(input_molecules[0]))
        monomer = input_molecules[-1]
        seeds = input_molecules[:-1]

        t1 = time.clock()
        aggregator.aggregate(seeds, monomer,
                             aggregate_size=size_of_aggregate,
                             hm_orientations=N,
                             method=method_args)
        print('Time:', time.clock() - t1)

    if args.react:
        if len(input_molecules) == 1:
            print('Reactor requires at least two molecules')
            sys.exit()
        if minimum_gamma is None or maximum_gamma is None:
            print('For a Reactor run specify the '
                  'values of gamma_min and gamma_max using \n'
                  '-gmin <integer> -gmax <integer>')
            sys.exit()
        if N is None:
            print("For aggregation, specify how many orientations"
                  "are to be used, by the argument\n"
                  "-N <number of orientations>")
        t1 = time.clock()
        reactor.react(input_molecules[0], input_molecules[1],
                      gamma_min=minimum_gamma, gamma_max=maximum_gamma,
                      hm_orientations=N, method=method_args)
        print('Time:', time.clock() - t1)
        return

    if args.optimize:
        import csv
        with open('energy.csv', 'w') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["Name", "Energy"])
            for each_mol in input_molecules:
                status = optimise(each_mol, method_args)
                if status is True:
                    writer.writerow([each_mol.name, each_mol.energy])
                else:
                    writer.writerow([each_mol.name, None])


if __name__ == "__main__":
    main()
