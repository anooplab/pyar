"""Tests for fixed, scheduled, and resistance-based bias controllers."""

import unittest

import numpy as np

from pyar.biases.controller import BiasController


class BiasControllerTests(unittest.TestCase):
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


if __name__ == "__main__":
    unittest.main()
