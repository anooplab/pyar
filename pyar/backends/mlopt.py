"""MLatom optimization backend command-line entrypoint."""

import argparse

import mlatom as ml


def main():
    parser = argparse.ArgumentParser(description='Optimize a molecule using a specified model.')
    parser.add_argument('input_molecule', metavar='input_molecule', type=str,
                        help='the path to the input molecule file')
    parser.add_argument('-c', '--charge', type=int, default=0,
                        help='the charge of the molecule')
    parser.add_argument('-m', '--multiplicity', type=int, default=1,
                        help='the multiplicity of the molecule')
    parser.add_argument('final_molecule', metavar='final_molecule', type=str,
                        help='the path to the final molecule file')
    args = parser.parse_args()

    initmol = ml.data.molecule.from_xyz_file(args.input_molecule)
    initmol.charge = args.charge
    initmol.multiplicity = args.multiplicity
    mymodel = ml.models.methods(method='AIQM1')
    geomopt = ml.optimize_geometry(model=mymodel, initial_molecule=initmol,
                                   program='Gaussian',
                                   maximum_number_of_steps=10000)
    final_mol = geomopt.optimized_molecule
    final_mol.write_file_with_xyz_coordinates(filename=args.final_molecule)
    print("Final energy: ", final_mol.energy)


if __name__ == '__main__':
    main()
