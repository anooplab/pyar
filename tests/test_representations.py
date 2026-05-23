import unittest
from unittest import mock

from dscribe.core.system import System

from pyar import representations


class RepresentationTests(unittest.TestCase):
    def test_mbtr_descriptor_uses_dscribe_system(self):
        captured = {}

        class FakeMBTR:
            def __init__(self, *args, **kwargs):
                pass

            def create(self, molecule):
                captured["is_system"] = isinstance(molecule, System)
                captured["symbols"] = molecule.get_chemical_symbols()
                return [0.0]

        with mock.patch.object(representations, "MBTR", FakeMBTR):
            result = representations.mbtr_descriptor(["H"], [[0.0, 0.0, 0.0]])

        self.assertTrue(captured["is_system"])
        self.assertEqual(captured["symbols"], ["H"])
        self.assertEqual(result, [0.0])

    def test_mbtr_descriptor_computes_with_real_dscribe_path(self):
        result = representations.mbtr_descriptor(
            ["H", "H"],
            [[0.0, 0.0, 0.0], [0.0, 0.0, 0.74]],
        )

        self.assertGreater(len(result), 0)
        self.assertTrue(any(value != 0.0 for value in result))

    def test_fingerprint_handles_overlapping_atoms(self):
        result = representations.fingerprint(
            ["H", "H"],
            [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
        )

        self.assertEqual(len(result), 2)
        self.assertTrue(all(value == value for value in result))
        self.assertTrue(all(abs(value) != float("inf") for value in result))


if __name__ == "__main__":
    unittest.main()
