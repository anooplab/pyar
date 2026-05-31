import unittest

from pyar.backend_capabilities import (
    BackendCapabilities,
    backend_family,
    backend_supports_biased_optimization,
    backend_supports_energy_gradient,
    backend_supports_geometry_optimization,
    backend_supports_method_basis_options,
    backend_supports_native_optimization,
    backend_supports_staged_optimization,
    get_backend_capabilities,
    normalize_backend_name,
    register_backend_capabilities,
    validate_backend_capability,
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
        self.assertTrue(backend_supports_native_optimization("xtb"))
        self.assertTrue(backend_supports_biased_optimization("xtb"))
        self.assertFalse(backend_supports_method_basis_options("xtb"))
        self.assertTrue(backend_supports_energy_gradient("xtb"))
        self.assertTrue(backend_supports_staged_optimization("xtb"))
        self.assertEqual(
            unsupported_qc_options("xtb", {"method", "opt_threshold", "nprocs"}),
            ["method"],
        )

    def test_backend_aliases_resolve_to_canonical_entries(self):
        self.assertEqual(normalize_backend_name("mlatom_aiqm1"), "aiqm1_mlatom")
        self.assertEqual(normalize_backend_name("xtb_aiqm1"), "xtb-aiqm1")
        self.assertEqual(normalize_backend_name("xtbturbo"), "xtb_turbo")
        self.assertTrue(get_backend_capabilities("mlatom_aiqm1").supports_method_basis_options is False)

    def test_backend_validation_reports_missing_feature_cleanly(self):
        with self.assertRaisesRegex(ValueError, "does not support method/basis options"):
            validate_backend_capability("xtb", ("method_basis_options",), context="backend routing")

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
