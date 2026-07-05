import unittest

from pyar.conformer import ConformerRequest, ConformerRequestError


class ConformerRequestTests(unittest.TestCase):
    def _request_kwargs(self, **overrides):
        kwargs = {
            "input_format": "auto",
            "num_conformers": 5,
            "top_n": 2,
            "backend_top_n": None,
            "num_seeds": 1,
            "diversity_fraction": 0.2,
            "compactness_fraction": 0.2,
            "rms_threshold": 0.25,
            "use_random_coords": True,
            "torsion_kicks": True,
            "torsion_rounds": 2,
            "torsion_mode": "random",
            "torsion_kicks_per_conformer": 6,
            "torsion_max_bonds": 3,
            "torsion_dedup_rms": 0.5,
            "force_field": "auto",
            "seed": 1,
            "num_threads": 0,
            "max_iterations": 200,
            "charge": None,
            "multiplicity": 1,
            "scftype": "rhf",
            "qc_params": None,
        }
        kwargs.update(overrides)
        return kwargs

    def test_from_options_normalizes_persisted_state_shape(self):
        request = ConformerRequest.from_options(
            "CCO",
            **self._request_kwargs(
                num_conformers="5",
                top_n="2",
                backend_top_n="7",
                num_seeds="3",
                diversity_fraction="0.25",
                compactness_fraction="0.5",
                rms_threshold="0.2",
                use_random_coords=1,
                torsion_kicks=0,
                torsion_mode="torsion-kick",
                torsion_kicks_per_conformer="4",
                torsion_max_bonds="2",
                torsion_dedup_rms="0.3",
                force_field="mmff",
                seed="10",
                num_threads="2",
                max_iterations="150",
                qc_params={"software": "xtb"},
            ),
        )

        self.assertEqual(request.num_conformers, 5)
        self.assertEqual(request.backend_top_n, 7)
        self.assertEqual(request.seed_values, (10, 11, 12))
        self.assertEqual(request.torsion_mode, "random")
        self.assertFalse(request.torsion_kicks)
        self.assertTrue(request.use_random_coords)
        self.assertTrue(request.backend_requested)
        state = request.to_state_dict()
        self.assertEqual(state["seed_values"], [10, 11, 12])
        self.assertEqual(state["generation_dedup_rms"], 0.5)
        self.assertEqual(state["backend_parameters"], {"software": "xtb"})

    def test_from_options_rejects_invalid_values(self):
        invalid_values = [
            ("num_conformers", 0, "--num-conformers must be at least 1"),
            ("top_n", 0, "--top-n must be at least 1"),
            ("num_seeds", 0, "--num-seeds must be at least 1"),
            ("rms_threshold", -0.1, "--rms-threshold must be non-negative"),
            ("backend_top_n", 0, "--backend-top-n must be at least 1"),
            ("diversity_fraction", 1.5, "--diversity-fraction must be between 0 and 1"),
            ("compactness_fraction", -0.1, "--compactness-fraction must be between 0 and 1"),
            ("torsion_mode", "grid", "--torsion-mode must be 'random'"),
        ]

        for option, value, message in invalid_values:
            with self.subTest(option=option):
                with self.assertRaisesRegex(ConformerRequestError, message):
                    ConformerRequest.from_options(
                        "CCO",
                        **self._request_kwargs(**{option: value}),
                    )


if __name__ == "__main__":
    unittest.main()
