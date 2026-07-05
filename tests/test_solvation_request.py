import unittest

from pyar.solvation import SolvationRequest, SolvationRequestError


class DummyMolecule:
    def __init__(self, name="seed", n_atoms=1):
        self.name = name
        self.atoms_list = ["H"] * n_atoms
        self.coordinates = [[0.0, 0.0, float(index)] for index in range(n_atoms)]
        self.charge = 0
        self.multiplicity = 1
        self.scftype = "rhf"
        self.fragments = [[index] for index in range(n_atoms)]


class SolvationRequestTests(unittest.TestCase):
    def _request_kwargs(self, **overrides):
        kwargs = {
            "seeds": [DummyMolecule("seed", n_atoms=2)],
            "monomer": DummyMolecule("solvent", n_atoms=1),
            "aggregate_size": 2,
            "hm_orientations": 8,
            "qc_params": {"software": "xtb"},
            "maximum_number_of_seeds": 4,
            "site": None,
            "connectivity_policy": "off",
        }
        kwargs.update(overrides)
        return kwargs

    def test_from_options_preserves_persisted_state_shape(self):
        request = SolvationRequest.from_options(
            **self._request_kwargs(
                aggregate_size="3",
                maximum_number_of_seeds="5",
                site=(0, 3),
                connectivity_policy="auto",
            )
        )

        state = request.to_state_dict()
        self.assertEqual(state["aggregate_size"], 3)
        self.assertEqual(state["orientations"], 8)
        self.assertEqual(state["backend_parameters"], {"software": "xtb"})
        self.assertEqual(state["maximum_number_of_seeds"], 5)
        self.assertEqual(state["site"], [0, 3])
        self.assertEqual(state["connectivity_policy"], "off")
        self.assertEqual(state["seeds"][0]["name"], "seed")
        self.assertEqual(state["seeds"][0]["atoms"], ["H", "H"])
        self.assertEqual(state["seeds"][0]["fragment_definition"], [[0], [1]])
        self.assertEqual(state["monomer"]["name"], "solvent")
        self.assertEqual(state["monomer"]["atoms"], ["H"])
        self.assertEqual(state["monomer"]["scftype"], "rhf")

    def test_from_options_accepts_auto_and_off_connectivity(self):
        for policy in ("auto", "off", None):
            with self.subTest(policy=policy):
                request = SolvationRequest.from_options(
                    **self._request_kwargs(connectivity_policy=policy)
                )
                self.assertEqual(request.connectivity_policy, "off")

    def test_from_options_rejects_invalid_values(self):
        invalid_values = [
            ({"seeds": []}, "Solvation requires at least one seed molecule"),
            ({"monomer": None}, "Solvation requires one solvent molecule"),
            ({"aggregate_size": 0}, "--solvation-size must be at least 1"),
            ({"maximum_number_of_seeds": 0}, "--maximum-number-of-seeds must be at least 1"),
            ({"connectivity_policy": "prefer"}, "require connectivity policy 'off'"),
            ({"connectivity_policy": "strict"}, "require connectivity policy 'off'"),
            ({"connectivity_policy": "sometimes"}, "Unknown connectivity policy"),
        ]

        for overrides, message in invalid_values:
            with self.subTest(overrides=overrides):
                with self.assertRaisesRegex(SolvationRequestError, message):
                    SolvationRequest.from_options(**self._request_kwargs(**overrides))


if __name__ == "__main__":
    unittest.main()
