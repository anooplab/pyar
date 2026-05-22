import unittest
from types import SimpleNamespace
from unittest import mock

from pyar.data_analysis import clustering


class ClusteringTests(unittest.TestCase):
    def test_choose_geometries_falls_back_when_mbtr_fails(self):
        molecules = [
            SimpleNamespace(name="a", atoms_list=["C"], coordinates=[[0.0, 0.0, 0.0]], energy=0.0),
            SimpleNamespace(name="b", atoms_list=["H"], coordinates=[[1.0, 0.0, 0.0]], energy=1.0),
        ]

        with mock.patch("pyar.representations.mbtr_descriptor", side_effect=TypeError("boom")):
            with mock.patch("pyar.data_analysis.clustering.remove_similar", return_value=molecules) as remover:
                result = clustering.choose_geometries(molecules, maximum_number_of_seeds=1)

        remover.assert_called_once_with(molecules)
        self.assertEqual(result, molecules)


if __name__ == "__main__":
    unittest.main()
