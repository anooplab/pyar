"""Scientific constants owned by PyAR, with explicit units and source profiles.

Distinct historical precisions and radius models intentionally retain their
original values. Consolidation does not authorize changing a Hamiltonian,
connectivity policy or conversion convention. ASE-derived tables originate in
pyar.data.atomic_data; AFIR radii follow Pyykko and Atsumi (2009).
"""
import numpy as np
from numpy import nan

BOHR_ANGSTROM_LEGACY = 0.52917726
BOHR_ANGSTROM_CODATA = 0.529177210903
BOHR_ANGSTROM_BACKEND = 0.52917721092
BOHR_ANGSTROM_ROUNDED = 0.529177
HARTREE_EV = 27.211386245988
EV_HARTREE_ROUNDED = 0.0367493
HARTREE_KJ_MOL = 2625.5
HARTREE_KCAL_MOL_LEGACY = 627.509541
HARTREE_KCAL_MOL_SCAN = 627.509474
HARTREE_KCAL_MOL_REPORT = 627.51
HARTREE_KCAL_MOL_ROUNDED = 627.509
BOND_DISTANCE_FACTOR = 1.3
GRAPH_BOND_FACTOR = 1.15
TRACE_BOND_OFF_FACTOR = 1.35

# Preserved from pyar/data/atomic_data.py; values and table coverage unchanged.
ASE_COVALENT_RADII_ANGSTROM = np.array([
    0.2,  # X
    0.31,  # H
    0.28,  # He
    1.28,  # Li
    0.96,  # Be
    0.84,  # B
    0.76,  # C
    0.71,  # N
    0.66,  # O
    0.57,  # F
    0.58,  # Ne
    1.66,  # Na
    1.41,  # Mg
    1.21,  # Al
    1.11,  # Si
    1.07,  # P
    1.05,  # S
    1.02,  # Cl
    1.06,  # Ar
    2.03,  # K
    1.76,  # Ca
    1.70,  # Sc
    1.60,  # Ti
    1.53,  # V
    1.39,  # Cr
    1.39,  # Mn
    1.32,  # Fe
    1.26,  # Co
    1.24,  # Ni
    1.32,  # Cu
    1.22,  # Zn
    1.22,  # Ga
    1.20,  # Ge
    1.19,  # As
    1.20,  # Se
    1.20,  # Br
    1.16,  # Kr
    2.20,  # Rb
    1.95,  # Sr
    1.90,  # Y
    1.75,  # Zr
    1.64,  # Nb
    1.54,  # Mo
    1.47,  # Tc
    1.46,  # Ru
    1.42,  # Rh
    1.39,  # Pd
    1.45,  # Ag
    1.44,  # Cd
    1.42,  # In
    1.39,  # Sn
    1.39,  # Sb
    1.38,  # Te
    1.39,  # I
    1.40,  # Xe
    2.44,  # Cs
    2.15,  # Ba
    2.07,  # La
    2.04,  # Ce
    2.03,  # Pr
    2.01,  # Nd
    1.99,  # Pm
    1.98,  # Sm
    1.98,  # Eu
    1.96,  # Gd
    1.94,  # Tb
    1.92,  # Dy
    1.92,  # Ho
    1.89,  # Er
    1.90,  # Tm
    1.87,  # Yb
    1.87,  # Lu
    1.75,  # Hf
    1.70,  # Ta
    1.62,  # W
    1.51,  # Re
    1.44,  # Os
    1.41,  # Ir
    1.36,  # Pt
    1.36,  # Au
    1.32,  # Hg
    1.45,  # Tl
    1.46,  # Pb
    1.48,  # Bi
    1.40,  # Po
    1.50,  # At
    1.50,  # Rn
    2.60,  # Fr
    2.21,  # Ra
    2.15,  # Ac
    2.06,  # Th
    2.00,  # Pa
    1.96,  # U
    1.90,  # Np
    1.87,  # Pu
    1.80,  # Am
    1.69,  # Cm
    0.2,  # Bk
    0.2,  # Cf
    0.2,  # Es
    0.2,  # Fm
    0.2,  # Md
    0.2,  # No
    0.2,  # Lr
])

# Preserved from pyar/data/atomic_data.py; values and table coverage unchanged.
ASE_VDW_RADII_ANGSTROM = np.array([
    np.nan,  # X
    1.20,  # H
    1.40,  # He [1]
    1.82,  # Li [1]
    1.53,  # Be [5]
    1.92,  # B [5]
    1.70,  # C [1]
    1.55,  # N [1]
    1.52,  # O [1]
    1.47,  # F [1]
    1.54,  # Ne [1]
    2.27,  # Na [1]
    1.73,  # Mg [1]
    1.84,  # Al [5]
    2.10,  # Si [1]
    1.80,  # P [1]
    1.80,  # S [1]
    1.75,  # Cl [1]
    1.88,  # Ar [1]
    2.75,  # K [1]
    2.31,  # Ca [5]
    np.nan,  # Sc
    np.nan,  # Ti
    np.nan,  # V
    np.nan,  # Cr
    np.nan,  # Mn
    np.nan,  # Fe
    np.nan,  # Co
    1.63,  # Ni [1]
    1.40,  # Cu [1]
    1.39,  # Zn [1]
    1.87,  # Ga [1]
    2.11,  # Ge [5]
    1.85,  # As [1]
    1.90,  # Se [1]
    1.85,  # Br [1]
    2.02,  # Kr [1]
    3.03,  # Rb [5]
    2.49,  # Sr [5]
    np.nan,  # Y
    np.nan,  # Zr
    np.nan,  # Nb
    np.nan,  # Mo
    np.nan,  # Tc
    np.nan,  # Ru
    np.nan,  # Rh
    1.63,  # Pd [1]
    1.72,  # Ag [1]
    1.58,  # Cd [1]
    1.93,  # In [1]
    2.17,  # Sn [1]
    2.06,  # Sb [5]
    2.06,  # Te [1]
    1.98,  # I [1]
    2.16,  # Xe [1]
    3.43,  # Cs [5]
    2.49,  # Ba [5]
    np.nan,  # La
    np.nan,  # Ce
    np.nan,  # Pr
    np.nan,  # Nd
    np.nan,  # Pm
    np.nan,  # Sm
    np.nan,  # Eu
    np.nan,  # Gd
    np.nan,  # Tb
    np.nan,  # Dy
    np.nan,  # Ho
    np.nan,  # Er
    np.nan,  # Tm
    np.nan,  # Yb
    np.nan,  # Lu
    np.nan,  # Hf
    np.nan,  # Ta
    np.nan,  # W
    np.nan,  # Re
    np.nan,  # Os
    np.nan,  # Ir
    1.75,  # Pt [1]
    1.66,  # Au [1]
    1.55,  # Hg [1]
    1.96,  # Tl [1]
    2.02,  # Pb [1]
    2.07,  # Bi [5]
    1.97,  # Po [5]
    2.02,  # At [5]
    2.20,  # Rn [5]
    3.48,  # Fr [5]
    2.83,  # Ra [5]
    np.nan,  # Ac
    np.nan,  # Th
    np.nan,  # Pa
    1.86,  # U [1]
    np.nan,  # Np
    np.nan,  # Pu
    np.nan,  # Am
    np.nan,  # Cm
    np.nan,  # Bk
    np.nan,  # Cf
    np.nan,  # Es
    np.nan,  # Fm
    np.nan,  # Md
    np.nan,  # No
    np.nan])  # Lr

# Preserved from pyar/data/new_atomic_data.py; values and table coverage unchanged.
COVALENT_RADII_ANGSTROM  = {
  'X': 0.0,
 'Ac': 1.86,
 'Ag': 1.28,
 'Al': 1.26,
 'Ar': 0.96,
 'As': 1.21,
 'At': 1.47,
 'Au': 1.24,
 'B': 0.85,
 'Ba': 1.96,
 'Be': 1.02,
 'Bi': 1.51,
 'Br': 1.14,
 'C': 0.75,
 'Ca': 1.71,
 'Cd': 1.36,
 'Ce': 1.63,
 'Cl': 0.99,
 'Co': 1.11,
 'Cr': 1.22,
 'Cs': 2.32,
 'Cu': 1.12,
 'Dy': 1.67,
 'Er': 1.65,
 'Eu': 1.68,
 'F': 0.64,
 'Fe': 1.16,
 'Fr': 2.23,
 'Ga': 1.24,
 'Gd': 1.69,
 'Ge': 1.21,
 'H': 0.32,
 'He': 0.46,
 'Hf': 1.52,
 'Hg': 1.33,
 'Ho': 1.66,
 'I': 1.33,
 'In': 1.42,
 'Ir': 1.22,
 'K': 1.96,
 'Kr': 1.17,
 'La': 1.8,
 'Li': 1.33,
 'Lu': 1.62,
 'Mg': 1.39,
 'Mn': 1.19,
 'Mo': 1.38,
 'N': 0.71,
 'Na': 1.55,
 'Nb': 1.47,
 'Nd': 1.74,
 'Ne': 0.67,
 'Ni': 1.1,
 'O': 0.63,
 'Os': 1.29,
 'P': 1.11,
 'Pa': nan,
 'Pb': 1.44,
 'Pd': 1.2,
 'Pm': 1.73,
 'Po': 1.45,
 'Pr': 1.76,
 'Pt': 1.23,
 'Ra': 2.01,
 'Rb': 2.1,
 'Re': 1.31,
 'Rh': 1.25,
 'Rn': 1.42,
 'Ru': 1.25,
 'S': 1.03,
 'Sb': 1.4,
 'Sc': 1.48,
 'Se': 1.16,
 'Si': 1.16,
 'Sm': 1.72,
 'Sn': 1.4,
 'Sr': 1.85,
 'Ta': 1.46,
 'Tb': 1.68,
 'Tc': 1.28,
 'Te': 1.36,
 'Th': 1.75,
 'Ti': 1.36,
 'Tl': 1.44,
 'Tm': 1.64,
 'V': 1.34,
 'W': 1.37,
 'Xe': 1.31,
 'Y': 1.63,
 'Yb': 1.7,
 'Zn': 1.18,
 'Zr': 1.54}

# Preserved from pyar/biases/afir.py; values and table coverage unchanged.
AFIR_COVALENT_RADII_ANGSTROM = {'x ': 0.00,
                  'h': 0.000001, 'he': 0.46, 'li': 1.20, 'be': 0.94, 'b': 0.77,
                  'c': 0.75, 'n': 0.71, 'o': 0.63, 'f': 0.64, 'ne': 0.67,
                  'na': 1.40, 'mg': 1.25, 'al': 1.13, 'si': 1.04, 'p': 1.10,
                  's': 1.02, 'cl': 0.99, 'ar': 0.96, 'k': 1.76, 'ca': 1.54,
                  'sc': 1.33, 'ti': 1.22, 'v': 1.21, 'cr': 1.10, 'mn': 1.07,
                  'fe': 1.04, 'co': 1.00, 'ni': 0.99, 'cu': 1.01, 'zn': 1.09,
                  'ga': 1.12, 'ge': 1.09, 'as': 1.15, 'se': 1.10, 'br': 1.14,
                  'kr': 1.17, 'rb': 1.89, 'sr': 1.67, 'y': 1.47, 'zr': 1.39,
                  'nb': 1.32, 'mo': 1.24, 'tc': 1.15, 'ru': 1.13, 'rh': 1.13,
                  'pd': 1.08, 'ag': 1.15, 'cd': 1.23, 'in': 1.28, 'sn': 1.26,
                  'sb': 1.26, 'te': 1.23, 'i': 1.32, 'xe': 1.31, 'cs': 2.09,
                  'ba': 1.76, 'la': 1.62, 'ce': 1.47, 'pr': 1.58, 'nd': 1.57,
                  'pm': 1.56, 'sm': 1.55, 'eu': 1.51, 'gd': 1.52, 'tb': 1.51,
                  'dy': 1.50, 'ho': 1.49, 'er': 1.49, 'tm': 1.48, 'yb': 1.53,
                  'lu': 1.46, 'hf': 1.37, 'ta': 1.31, 'w': 1.23, 're': 1.18,
                  'os': 1.16, 'ir': 1.11, 'pt': 1.12, 'au': 1.13, 'hg': 1.32,
                  'tl': 1.30, 'pb': 1.30, 'bi': 1.36, 'po': 1.31, 'at': 1.38,
                  'rn': 1.42, 'fr': 2.01, 'ra': 1.81, 'ac': 1.67, 'th': 1.58,
                  'pa': 1.52, 'u': 1.53, 'np': 1.54, 'pu': 1.55}

# Preserved element-keyed van der Waals model.
VDW_RADII_ANGSTROM  = {
  'X': 0.0,
 'Ac': 2.47,
 'Ag': 2.11,
 'Al': 1.84,
 'Ar': 1.88,
 'As': 1.85,
 'At': 2.02,
 'Au': 2.14,
 'B': 1.92,
 'Ba': 2.68,
 'Be': 1.53,
 'Bi': 2.07,
 'Br': 1.85,
 'C': 1.7,
 'Ca': 2.31,
 'Cd': 2.18,
 'Ce': 2.42,
 'Cl': 1.75,
 'Co': 2.0,
 'Cr': 2.06,
 'Cs': 3.43,
 'Cu': 1.96,
 'Dy': 2.31,
 'Er': 2.29,
 'Eu': 2.35,
 'F': 1.47,
 'Fe': 2.04,
 'Fr': 3.48,
 'Ga': 1.87,
 'Gd': 2.34,
 'Ge': 2.11,
 'H': 1.1,
 'He': 1.4,
 'Hf': 2.23,
 'Hg': 2.23,
 'Ho': 2.3,
 'I': 1.98,
 'In': 1.93,
 'Ir': 2.13,
 'K': 2.75,
 'Kr': 2.02,
 'La': 2.43,
 'Li': 1.82,
 'Lu': 2.24,
 'Mg': 1.73,
 'Mn': 2.05,
 'Mo': 2.17,
 'N': 1.55,
 'Na': 2.27,
 'Nb': 2.18,
 'Nd': 2.39,
 'Ne': 1.54,
 'Ni': 1.97,
 'O': 1.52,
 'Os': 2.16,
 'P': 1.8,
 'Pa': nan,
 'Pb': 2.02,
 'Pd': 2.1,
 'Pm': 2.38,
 'Po': 1.97,
 'Pr': 2.4,
 'Pt': 2.13,
 'Ra': 2.83,
 'Rb': 3.03,
 'Re': 2.16,
 'Rh': 2.1,
 'Rn': 2.2,
 'Ru': 2.13,
 'S': 1.8,
 'Sb': 2.06,
 'Sc': 2.15,
 'Se': 1.9,
 'Si': 2.1,
 'Sm': 2.36,
 'Sn': 2.17,
 'Sr': 2.49,
 'Ta': 2.22,
 'Tb': 2.33,
 'Tc': 2.16,
 'Te': 2.06,
 'Th': 2.45,
 'Ti': 2.11,
 'Tl': 1.96,
 'Tm': 2.27,
 'V': 2.07,
 'W': 2.18,
 'Xe': 2.16,
 'Y': 2.32,
 'Yb': 2.26,
 'Zn': 2.01,
 'Zr': 2.23}
