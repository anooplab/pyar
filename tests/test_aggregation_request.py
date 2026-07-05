import unittest

from pyar.aggregation import AggregateRequest, AggregateRequestError


class DummyMolecule:
    def __init__(self, name="seed"):
        self.name = name
        self.atoms_list = ["C", "H"]
        self.coordinates = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]
        self.charge = 0
        self.multiplicity = 1
        self.scftype = "rhf"
        self.fragments = [0, 0]


class AggregateRequestTests(unittest.TestCase):
    def _request_kwargs(self, **overrides):
        kwargs = {
            "molecules": [DummyMolecule()],
            "aggregate_sizes": [1],
            "hm_orientations": 8,
            "qc_params": {"software": "xtb"},
            "maximum_number_of_seeds": 4,
            "first_pathway": 0,
            "number_of_pathways": 2,
            "site": None,
            "connectivity_policy": "auto",
        }
        kwargs.update(overrides)
        return kwargs

    def test_from_options_preserves_persisted_state_shape(self):
        request = AggregateRequest.from_options(
            **self._request_kwargs(
                aggregate_sizes=["2"],
                maximum_number_of_seeds="5",
                first_pathway="1",
                number_of_pathways="3",
                site=(0, 2),
                connectivity_policy=None,
            )
        )

        state = request.to_state_dict()
        self.assertEqual(state["aggregate_sizes"], [2])
        self.assertEqual(state["orientations"], 8)
        self.assertEqual(state["backend_parameters"], {"software": "xtb"})
        self.assertEqual(state["maximum_number_of_seeds"], 5)
        self.assertEqual(state["first_pathway"], 1)
        self.assertEqual(state["number_of_pathways"], 3)
        self.assertEqual(state["site"], [0, 2])
        self.assertEqual(state["connectivity_policy"], "auto")
        self.assertEqual(state["fragments"][0]["atoms"], ["C", "H"])
        self.assertEqual(state["fragments"][0]["coordinates"], [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
        self.assertEqual(state["fragments"][0]["charge"], 0)
        self.assertEqual(state["fragments"][0]["multiplicity"], 1)
        self.assertEqual(state["fragments"][0]["scftype"], "rhf")
        self.assertEqual(state["fragments"][0]["fragment_definition"], [0, 0])

    def test_from_options_accepts_explicit_connectivity_policies(self):
        for policy in ("auto", "off", "prefer", "strict"):
            with self.subTest(policy=policy):
                request = AggregateRequest.from_options(
                    **self._request_kwargs(connectivity_policy=policy)
                )
                self.assertEqual(request.connectivity_policy, policy)

    def test_from_options_rejects_invalid_values(self):
        invalid_values = [
            ({"molecules": []}, "Aggregation requires at least one molecule"),
            ({"aggregate_sizes": []}, "Aggregate sizes must be specified for every molecule"),
            ({"aggregate_sizes": [0]}, "Aggregate sizes must be positive integers"),
            ({"maximum_number_of_seeds": 0}, "--maximum-number-of-seeds must be at least 1"),
            ({"first_pathway": -1}, "--first-pathway must be non-negative"),
            ({"number_of_pathways": 0}, "--number-of-pathways must be at least 1"),
            ({"connectivity_policy": "sometimes"}, "Unknown connectivity policy"),
        ]

        for overrides, message in invalid_values:
            with self.subTest(overrides=overrides):
                with self.assertRaisesRegex(AggregateRequestError, message):
                    AggregateRequest.from_options(**self._request_kwargs(**overrides))


if __name__ == "__main__":
    unittest.main()
