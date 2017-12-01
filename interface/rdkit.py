"""
mopac.py - interface to mopac program
"""
'''
Copyright (C) 2016 by Surajit Nandi, Anoop Ayyappan, and Mark P. Waller
Indian Institute of Technology Kharagpur, India and Westfaelische Wilhelms
Universitaet Muenster, Germany

This file is part of the PyAR project.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
'''

import os
import subprocess as subp

import numpy as np

import babel

from rdkit import Chem
from rdkit.Chem import AllChem




#babel.xyz_to_sdf_file(xyz_input_files, sdf_output_file):

sdf_output_file = 'test.sdf'
suppl = Chem.SDMolSupplier(sdf_output_file)
ms = [x for x in suppl if x is not None]
for m in ms:
    AllChem.UFFOptimizeMolecule(m,ignoreInterfragInteractions=False)
    ff=AllChem.UFFGetMoleculeForceField(m)
    energy=ff.CalcEnergy()
#    print energy
