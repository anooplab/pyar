"""Golden values at 9f1e567: relocation must not change any radius model."""
import hashlib
import json
import numpy as np
import pytest
from pyar import constants

RADIUS_HASHES = {'ASE_COVALENT_RADII_ANGSTROM': 'adfc463937e9eacafc9ded0da8ad1b1911529ee42028bc4c087089c03334178f', 'ASE_VDW_RADII_ANGSTROM': '065b718b1b4c00eb1f1bdcfa9d4fd7f069a36ff60ce661011878a14e308b86a4', 'COVALENT_RADII_ANGSTROM': '208313ceb42f9709944053526a33af2795d78d9c30f8a314b7f696ad48c9df38', 'VDW_RADII_ANGSTROM': 'ca2400d2943c13cdd80c0c29c4112d21836b84ea27488893daf73187639e029b', 'AFIR_COVALENT_RADII_ANGSTROM': '0c5e3481747156bfd61bf7803dd2f018dffc8a11821e8793b48c49e09c480e07'}

@pytest.mark.parametrize('name,expected', RADIUS_HASHES.items())
def test_radius_profile_matches_scientific_baseline(name, expected):
    value = getattr(constants, name)
    normalized = value.tolist() if isinstance(value, np.ndarray) else value
    assert hashlib.sha256(json.dumps(normalized, sort_keys=True).encode()).hexdigest() == expected


def test_conversion_profiles_retain_original_precision():
    from pyar.data.units import angstrom2bohr, bohr2angstrom
    assert bohr2angstrom(1.) == .52917726
    assert angstrom2bohr(1.) == 1. / .52917726
    assert constants.BOHR_ANGSTROM_CODATA == .529177210903
    assert constants.HARTREE_KCAL_MOL_SCAN == 627.509474
    assert constants.HARTREE_KCAL_MOL_REPORT == 627.51
