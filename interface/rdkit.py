"""
Not completed
"""

from rdkit.Chem import AllChem

from rdkit import Chem

# babel.xyz_to_sdf_file(xyz_input_files, sdf_output_file):

sdf_output_file = 'test.sdf'
suppl = Chem.SDMolSupplier(sdf_output_file)
ms = [x for x in suppl if x is not None]
for m in ms:
    AllChem.UFFOptimizeMolecule(m, ignoreInterfragInteractions=False)
    ff = AllChem.UFFGetMoleculeForceField(m)
    energy = ff.CalcEnergy()
# print energy
