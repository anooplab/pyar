def angstrom2bohr(input_value):
    # 1 bohr = 0.52917726 angstrom
    return input_value / 0.52917726


def bohr2angstrom(input_value):
    # 1 bohr = 0.52917726 angstrom
    return input_value * 0.52917726


def kilojoules2atomic_units(input_value):
    # 1 a. u. =  2625.5 kj/my_mol
    return float(input_value) / 2625.5


def atomic_unit2kilo_calories(input_value):
    # 1 a.u. = 627.509541 kcal/my_mol
    return float(input_value) * 627.509541
