import argparse
import pyar.mlatom as ml

# Create the parser
parser = argparse.ArgumentParser(description='Optimize a molecule using a specified model.')

# Add the arguments
parser.add_argument('input_molecule', metavar='input_molecule', type=str, help='the path to the input molecule file')
parser.add_argument('-c', '--charge', type=int, default=0, help='the charge of the molecule')
parser.add_argument('-m', '--multiplicity', type=int, default=1, help='the multiplicity of the molecule')
parser.add_argument('final_molecule', metavar='final_molecule', type=str, help='the path to the final molecule file')

# Parse the arguments
args = parser.parse_args()

# Get the initial guess for the molecules to optimize
initmol = ml.data.molecule.from_xyz_file(args.input_molecule)
initmol.charge = args.charge
initmol.multiplicity = args.multiplicity

# Choose one of the predefined (automatically recognized) methods
mymodel = ml.models.methods(method='AIQM1')

# Optimize the geometry with the chosen optimizer
geomopt = ml.optimize_geometry(model=mymodel, initial_molecule=initmol, program='Gaussian', maximum_number_of_steps=10000)

# Get the final geometry
final_mol = geomopt.optimized_molecule

# Check and save the final geometry
final_mol.write_file_with_xyz_coordinates(filename=args.final_molecule)
print("Final energy: ", final_mol.energy)