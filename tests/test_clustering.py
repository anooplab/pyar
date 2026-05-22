import unittest
from types import SimpleNamespace
from unittest import mock

from pyar.data_analysis import clustering


class ClusteringTests(unittest.TestCase):
    def setUp(self):
        clustering._MBTR_RUNTIME_DISABLED = False
        clustering._MBTR_DISABLE_REASON = None

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

    def test_mbtr_is_disabled_after_first_runtime_failure(self):
        molecules = [
            SimpleNamespace(name="a", atoms_list=["C"], coordinates=[[0.0, 0.0, 0.0]], energy=0.0),
            SimpleNamespace(name="b", atoms_list=["H"], coordinates=[[1.0, 0.0, 0.0]], energy=1.0),
        ]

        with mock.patch("pyar.representations.mbtr_descriptor", side_effect=TypeError("mbtr boom")) as mbtr_mock:
            with mock.patch("pyar.data_analysis.clustering.remove_similar", return_value=molecules):
                clustering.choose_geometries(molecules, maximum_number_of_seeds=1)

        self.assertEqual(mbtr_mock.call_count, 1)

        with mock.patch("pyar.representations.mbtr_descriptor") as mbtr_mock:
            with mock.patch("pyar.data_analysis.clustering.remove_similar", return_value=molecules):
                result = clustering.choose_geometries(molecules, maximum_number_of_seeds=1)

        mbtr_mock.assert_not_called()
        self.assertEqual(result, molecules)


if __name__ == "__main__":
    unittest.main()
