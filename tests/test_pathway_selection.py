"""Aggregation addition orders must be reproducible across fresh runs."""
from pyar.core.molecule import Molecule
from pyar.growth.service import select_pathways


def _monomers():
    a = Molecule(['H'], [[0, 0, 0]], name='a')
    b = Molecule(['He'], [[0, 0, 0]], name='b')
    return [a, a, b, b]


def _labels(pathways):
    return [''.join(monomer.name for monomer in path) for path in pathways]


def test_select_pathways_is_seeded_ordered_and_distinct():
    first = _labels(select_pathways(_monomers(), 4))
    second = _labels(select_pathways(_monomers(), 4))
    assert first == second
    assert len(first) == len(set(first)) == 4


def test_select_pathways_stops_at_the_number_of_distinct_orders():
    # aabb has six distinct orders; asking for more must not loop or repeat.
    pathways = _labels(select_pathways(_monomers(), 50))
    assert sorted(pathways) == sorted({'aabb', 'abab', 'abba', 'baab', 'baba', 'bbaa'})
