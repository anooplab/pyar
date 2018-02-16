def angstrom2bohr(input):
    # 1 bohr = 0.52917726 angstrom
    return input/0.52917726


def bohr2angstrom(input):
    # 1 bohr = 0.52917726 angstrom
    return input * 0.52917726


def kjoules2atomicunits(input):
    # 1 a. u. =  2625.5 kj/mol
    return float(input) / 2625.5


def atomicunit2kcal(input):
    # 1 a.u. = 627.509541 kcal/mol
    return float(input) * 627.509541
