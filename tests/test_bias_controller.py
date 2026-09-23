"""Tests for fixed, scheduled, and resistance-based bias controllers."""

import unittest

import numpy as np

from pyar.biases.controller import (
    BiasController,
    FixedBiasController,
    resolve_controller_policy,
)


class BiasControllerTests(unittest.TestCase):
    def test_cli_controller_policy_is_inferred_or_validated_from_tuning_options(self):
        self.assertEqual(resolve_controller_policy(), "fixed")
        self.assertEqual(resolve_controller_policy(alpha_min=0.1), "adaptive")
        self.assertEqual(resolve_controller_policy(scheduled_alpha=0.5), "scheduled")
        with self.assertRaisesRegex(ValueError, "require --bias-controller adaptive"):
            resolve_controller_policy("fixed", smoothing=0.5)
        with self.assertRaisesRegex(ValueError, "cannot be combined"):
            resolve_controller_policy(smoothing=0.5, scheduled_alpha=0.5)

    def test_explicit_fixed_controller_matches_legacy_fixed_policy(self):
        gradient = np.asarray([[0.2, -0.1, 0.4]])
        coordinate = np.asarray([[0.0, 1.0, 0.0]])
        explicit = FixedBiasController()
        compatible = BiasController("fixed")
        for controller in (explicit, compatible):
            decision = controller.start_segment(gradient, coordinate, 0.75, 1.25)
            self.assertEqual(decision.alpha, 1.25)
            self.assertEqual(decision.alpha_target, 1.25)
            self.assertIsNone(decision.alpha_critical)
            self.assertEqual(controller.state_dict()["alpha_max"], 1.25)
        self.assertEqual(explicit.state_dict(), compatible.state_dict())

    def test_fixed_and_scheduled_policies_respect_bounds(self):
        gradient = np.zeros((2, 3))
        coordinate = np.ones((2, 3))
        self.assertEqual(BiasController("fixed").select(gradient, coordinate, 2.0).alpha, 2.0)
        decision = BiasController("scheduled", scheduled_alpha=3.0, alpha_min=0.5).select(
            gradient, coordinate, 2.0
        )
        self.assertEqual(decision.alpha, 2.0)

    def test_adaptive_policy_uses_resistance_estimate_and_smoothing(self):
        controller = BiasController("adaptive", safety_margin=0.1, smoothing=0.5)
        coordinate = np.asarray([[1.0, 0.0, 0.0]])
        first = controller.start_segment(np.asarray([[-0.5, 0.0, 0.0]]), coordinate, 2.0, 1.0)
        second = controller.start_segment(np.asarray([[-0.9, 0.0, 0.0]]), coordinate, 1.5, 1.0)
        self.assertAlmostEqual(first.alpha_critical, 0.5)
        self.assertAlmostEqual(first.alpha, 0.6)
        self.assertAlmostEqual(second.alpha_target, 1.0)
        self.assertAlmostEqual(second.alpha, 0.8)

    def test_proposals_do_not_advance_history(self):
        controller = BiasController("adaptive", smoothing=0.5)
        qgrad = np.array([[1., 0., 0.]])
        controller.start_segment(-0.5*qgrad, qgrad, 2., 1.)
        state = controller.state_dict()
        first = controller.select(-0.9*qgrad, qgrad, 1.)
        self.assertEqual(first, controller.select(-0.9*qgrad, qgrad, 1.))
        self.assertEqual(state, controller.state_dict())

    def test_checkpoint_restores_smoothing_and_energy_continuity(self):
        controller = BiasController("adaptive", smoothing=0.5)
        qgrad = np.array([[1., 0., 0.]])
        first = controller.start_segment(-0.5*qgrad, qgrad, 2., 1.)
        second = controller.start_segment(-0.9*qgrad, qgrad, 1.5, 1.)
        self.assertAlmostEqual(first.alpha*1.5, second.alpha*1.5 + controller.energy_offset)
        restored = BiasController("adaptive", smoothing=0.5)
        restored.load_state_dict(controller.state_dict(), 1.)
        self.assertEqual(controller.state_dict(), restored.state_dict())
        for active in (controller, restored):
            active.start_segment(-0.8*qgrad, qgrad, 1., 1.)
        self.assertEqual(controller.state_dict(), restored.state_dict())
        with self.assertRaisesRegex(ValueError, "configuration"):
            BiasController("adaptive").load_state_dict(controller.state_dict(), 1.)

    def test_adaptive_synthetic_sign_zero_gradient_and_bound_cases(self):
        coordinate = np.asarray([[1.0, 0.0, 0.0]])
        cases = (
            (np.asarray([[-3.0, 0.0, 0.0]]), 3.0),  # PES opposes progress
            (np.asarray([[3.0, 0.0, 0.0]]), 0.0),   # PES supports progress
            (np.asarray([[0.0, 2.0, 0.0]]), 0.0),   # orthogonal
            (np.zeros((1, 3)), 0.0),
        )
        for physical_gradient, expected in cases:
            decision = BiasController("adaptive").select(physical_gradient, coordinate, 2.0)
            self.assertAlmostEqual(decision.alpha_critical, expected)
        zero_coordinate = np.zeros((1, 3))
        decision = BiasController("adaptive", epsilon=1e-8).select(
            np.ones((1, 3)), zero_coordinate, 2.0
        )
        self.assertEqual(decision.alpha_critical, 0.0)
        self.assertTrue(np.isfinite(decision.alpha))
        clipped = BiasController("adaptive").select(
            np.asarray([[-100.0, 0.0, 0.0]]), coordinate, 0.25
        )
        self.assertEqual(clipped.alpha, 0.25)
        lower_clipped = BiasController("adaptive", alpha_min=0.2).select(
            np.asarray([[100.0, 0.0, 0.0]]), coordinate, 1.0
        )
        self.assertEqual(lower_clipped.alpha, 0.2)
        self.assertEqual(BiasController("adaptive").select(
            np.asarray([[-100.0, 0.0, 0.0]]), coordinate, 0.0
        ).alpha, 0.0)

    def test_adaptive_rejects_nonfinite_and_mismatched_gradients(self):
        controller = BiasController("adaptive")
        with self.assertRaisesRegex(ValueError, "finite arrays"):
            controller.select(np.asarray([[np.nan, 0., 0.]]), np.ones((1, 3)), 1.0)
        with self.assertRaisesRegex(ValueError, "identical shapes"):
            controller.select(np.ones((1, 3)), np.ones((2, 3)), 1.0)

    def test_smoothing_must_be_positive_and_at_most_one(self):
        for smoothing in (0.0, 1.01):
            with self.subTest(smoothing=smoothing), self.assertRaisesRegex(
                ValueError, "smoothing must be greater than 0"
            ):
                BiasController("adaptive", smoothing=smoothing)

    def test_checkpoint_rejects_changed_alpha_ceiling(self):
        controller = BiasController("adaptive")
        coordinate = np.asarray([[1.0, 0.0, 0.0]])
        controller.start_segment(np.zeros((1, 3)), coordinate, 1.0, 1.0)
        restored = BiasController("adaptive")
        with self.assertRaisesRegex(ValueError, "alpha_max"):
            restored.load_state_dict(controller.state_dict(), 0.5)


if __name__ == "__main__":
    unittest.main()
