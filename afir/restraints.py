from __future__ import print_function
import sys
import numpy as np
from itertools import product
from math import sqrt
from units import *
import os

# covalent radii (taken from Pyykko and Atsumi, Chem. Eur. J. 15, 2009, 188-197)
# values for metals decreased by 10 %
from units import kjoules2atomicunits

covalent_radii = {'x ': 0.00,
                  'h': 0.32, 'he': 0.46, 'li': 1.20, 'be': 0.94, 'b': 0.77, \
                  'c': 0.75, 'n': 0.71, 'o': 0.63, 'f': 0.64, 'ne': 0.67, \
                  'na': 1.40, 'mg': 1.25, 'al': 1.13, 'si': 1.04, 'p': 1.10, \
                  's': 1.02, 'cl': 0.99, 'ar': 0.96, 'k': 1.76, 'ca': 1.54, \
                  'sc': 1.33, 'ti': 1.22, 'v': 1.21, 'cr': 1.10, 'mn': 1.07, \
                  'fe': 1.04, 'co': 1.00, 'ni': 0.99, 'cu': 1.01, 'zn': 1.09, \
                  'ga': 1.12, 'ge': 1.09, 'as': 1.15, 'se': 1.10, 'br': 1.14, \
                  'kr': 1.17, 'rb': 1.89, 'sr': 1.67, 'y': 1.47, 'zr': 1.39, \
                  'nb': 1.32, 'mo': 1.24, 'tc': 1.15, 'ru': 1.13, 'rh': 1.13, \
                  'pd': 1.08, 'ag': 1.15, 'cd': 1.23, 'in': 1.28, 'sn': 1.26, \
                  'sb': 1.26, 'te': 1.23, 'i': 1.32, 'xe': 1.31, 'cs': 2.09, \
                  'ba': 1.76, 'la': 1.62, 'ce': 1.47, 'pr': 1.58, 'nd': 1.57, \
                  'pm': 1.56, 'sm': 1.55, 'eu': 1.51, 'gd': 1.52, 'tb': 1.51, \
                  'dy': 1.50, 'ho': 1.49, 'er': 1.49, 'tm': 1.48, 'yb': 1.53, \
                  'lu': 1.46, 'hf': 1.37, 'ta': 1.31, 'w': 1.23, 're': 1.18, \
                  'os': 1.16, 'ir': 1.11, 'pt': 1.12, 'au': 1.13, 'hg': 1.32, \
                  'tl': 1.30, 'pb': 1.30, 'bi': 1.36, 'po': 1.31, 'at': 1.38, \
                  'rn': 1.42, 'fr': 2.01, 'ra': 1.81, 'ac': 1.67, 'th': 1.58, \
                  'pa': 1.52, 'u': 1.53, 'np': 1.54, 'pu': 1.55}


def get_covalent_radius(z):
    return angstrom2bohr(covalent_radii[z.lower()])


def calculate_gradient(ac, ac2, alpha, v, w, p):
    zn, cn = ac
    alo, aco = ac2
    gr = np.zeros(3)
    for zi, ci in zip(alo, aco):
        Rni = get_covalent_radius(zn) + get_covalent_radius(zi)
        rni = np.linalg.norm(cn-ci)
        dq = cn-ci
        gr += (-alpha/w*(((1-p)*dq*Rni**p*rni**(-1-p))+(p*v/w*dq*Rni**p*rni**(-2-p))))
    return np.array(gr)


def isotropic(atoms_in_fragment, atoms_list, coords, force):

    assert len(atoms_in_fragment) == 2

    idwp = 6.0
    covalent_radii['h'] = 0.000001

    fragment1, fragment2 = [coords[flist,:] for flist in atoms_in_fragment]
    al1, al2 = [np.array(atoms_list)[flist] for flist in atoms_in_fragment]
    radius_one = [get_covalent_radius(z) for z in al1]
    radius_two = [get_covalent_radius(z) for z in al2]

    r_ij = np.array([np.linalg.norm(a - b) for a,b in product(fragment1, fragment2)])
    R_ij = np.array([(a+b) for a, b in product(radius_one, radius_two)])

    if bohr2angstrom(np.min(r_ij)) > 5.0:
        print('Warning! Fragments may be flying away!\nCheck the geometry in', os.getcwd())

    nom = np.sum((R_ij/r_ij)**idwp*r_ij)
    den = np.sum((R_ij/r_ij)**idwp)

    alpha = calculate_alpha(force)

    total_restraint_energy = alpha * nom/den
    gradients = np.zeros((len(coords),3))

    for index, ac in enumerate(zip(atoms_list,coords)):
        if index in atoms_in_fragment[0]:
            gradients[index] = calculate_gradient(ac, (al2, fragment2), alpha, nom, den, idwp)
        elif index in atoms_in_fragment[1]:
            gradients[index] = calculate_gradient(ac, (al2, fragment1), alpha, nom, den, idwp)
    total_restraint_gradients = np.sqrt(np.sum(gradients**2))
    return total_restraint_energy, total_restraint_gradients, gradients


def calculate_alpha(force):
    gamma = kjoules2atomicunits(force)
    epsilon = kjoules2atomicunits(1.0061)
    r_zero = angstrom2bohr(3.8164)
    #    eqn. 3 JCTC 2011,7,2335
    alpha = gamma / ((2 ** (-1.0 / 6.0) - (1 + sqrt(1 + gamma / epsilon)) ** (-1.0 / 6.0)) * r_zero)
    return alpha


def main():
    pass


if __name__ == '__main__':
    main()
