import unittest

from pyar.backend_capabilities import (
    BackendCapabilities,
    backend_family,
    backend_supports_geometry_optimization,
    backend_supports_staged_optimization,
    get_backend_capabilities,
    register_backend_capabilities,
    supported_geometry_backends,
    unsupported_qc_options,
)


class BackendCapabilitiesTests(unittest.TestCase):
    def test_geometry_backends_are_explicitly_registered(self):
        self.assertEqual(supported_geometry_backends(), ("xtb", "aimnet_2", "orca", "gaussian"))
        self.assertTrue(backend_supports_geometry_optimization("xtb"))
        self.assertTrue(backend_supports_geometry_optimization("aimnet_2"))
        self.assertTrue(backend_supports_geometry_optimization("orca"))
        self.assertTrue(backend_supports_geometry_optimization("gaussian"))
        self.assertFalse(backend_supports_geometry_optimization("psi4"))

    def test_backend_family_and_qc_option_filtering_remain_centralized(self):
        self.assertEqual(backend_family("xtb"), "semiempirical")
        self.assertTrue(backend_supports_staged_optimization("xtb"))
        self.assertEqual(
            unsupported_qc_options("xtb", {"method", "opt_threshold", "nprocs"}),
            ["method"],
        )

    def test_energy_gradient_flag_requires_a_registered_provider_for_route_support(self):
        original = get_backend_capabilities("fake")
        try:
            register_backend_capabilities(
                "fake",
                BackendCapabilities(
                    family="test",
                    energy_gradient=True,
                    supported_options=frozenset({"nprocs"}),
                ),
            )
            self.assertTrue(get_backend_capabilities("fake").energy_gradient)
            self.assertFalse(backend_supports_geometry_optimization("fake"))
        finally:
            register_backend_capabilities("fake", original)


if __name__ == "__main__":
    unittest.main()
