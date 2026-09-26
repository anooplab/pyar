"""Fixtures and contract tests for optional bond-order analysis."""

import pytest

from pyar.bonding_analysis import (
    parse_bonding_analysis,
    parse_orca_mayer_bond_orders,
    parse_xtb_wiberg_bond_orders,
)


ORCA_OUTPUT = """
*****************************
* MAYER POPULATION ANALYSIS *
*****************************
Mayer bond orders larger than 0.1
B( 0-C , 1-C ) : 1.3850 B( 0-C , 2-H ) : 0.9612
B( 1-C , 0-C ) : 1.3850
"""

XTB_OUTPUT = """
Wiberg/Mayer (AO) data.
largest (>0.10) Wiberg bond orders for each atom
     1  O   1.782        H    2 0.891    H    3 0.891
     2  H   0.892        O    1 0.891
     3  H   0.892        O    1 0.891
molecular dipole moment from electron density (au)
"""


def test_orca_parser_normalizes_symmetric_pairs_and_metadata():
    result = parse_orca_mayer_bond_orders(ORCA_OUTPUT, expected_atoms=3)

    assert result.available is True
    assert result.scheme == "mayer"
    assert result.bond_order(1, 0) == pytest.approx(1.3850)
    assert result.bond_order(0, 2) == pytest.approx(0.9612)
    assert result.metadata["source"] == "Mayer population analysis"


def test_xtb_parser_converts_one_based_indices():
    result = parse_xtb_wiberg_bond_orders(XTB_OUTPUT, expected_atoms=3)

    assert result.scheme == "wiberg_ao"
    assert result.bond_order(0, 1) == pytest.approx(0.891)
    assert result.bond_order(2, 0) == pytest.approx(0.891)


@pytest.mark.parametrize("parser", [parse_orca_mayer_bond_orders, parse_xtb_wiberg_bond_orders])
def test_parsers_reject_missing_analysis(parser):
    with pytest.raises(ValueError, match="section is missing"):
        parser("normal backend output", expected_atoms=2)


def test_parser_rejects_invalid_atom_indices():
    with pytest.raises(ValueError, match="exceeds"):
        parse_orca_mayer_bond_orders(
            ORCA_OUTPUT.replace("0-C , 1-C", "0-C , 4-C"),
            expected_atoms=3,
        )


def test_unavailable_backend_is_nonfatal_and_explicit():
    result = parse_bonding_analysis("gaussian", "Gaussian output", expected_atoms=2)

    assert result.available is False
    assert result.bond_orders == {}
    assert result.reason == "no_stable_parser_in_current_integration"


def test_malformed_supported_output_is_not_reported_as_unavailable():
    with pytest.raises(ValueError, match="entries are missing or malformed"):
        parse_bonding_analysis("orca", "MAYER POPULATION ANALYSIS", expected_atoms=2)
