import io
import unittest
from types import SimpleNamespace
from unittest import mock

import numpy as np

from pyar.sampling.metrics import quaternion_coverage_metrics, sphere_coverage_metrics
from pyar.sampling.rotation import halton_quaternions
from pyar.sampling.sphere import fibonacci_sphere, generate_directions, maximin_select
from pyar.sampling.trial_generator import generate_trial_vectors
from pyar.scripts import benchmark_orientations

orientation_sampling = SimpleNamespace(
    fibonacci_sphere=fibonacci_sphere,
    maximin_select=maximin_select,
    sphere_coverage_metrics=sphere_coverage_metrics,
    generate_directions=generate_directions,
    halton_quaternions=halton_quaternions,
    quaternion_coverage_metrics=quaternion_coverage_metrics,
    generate_trial_vectors=generate_trial_vectors,
)


class OrientationSamplingTests(unittest.TestCase):
    def test_fibonacci_points_are_deterministic_unit_vectors(self):
        first = orientation_sampling.fibonacci_sphere(8)
        second = orientation_sampling.fibonacci_sphere(8)

        np.testing.assert_allclose(first, second)
        np.testing.assert_allclose(np.linalg.norm(first, axis=1), np.ones(8))

    def test_maximin_reduces_candidate_pool_to_diverse_points(self):
        candidates = np.array([
            [0.0, 0.0, 1.0],
            [0.01, 0.0, 1.0],
            [0.0, 0.0, -1.0],
            [1.0, 0.0, 0.0],
        ])

        selected = orientation_sampling.maximin_select(candidates, 2)

        np.testing.assert_allclose(selected[0], [0.0, 0.0, 1.0])
        np.testing.assert_allclose(selected[1], [0.0, 0.0, -1.0])

    def test_coverage_metrics_reward_antipodal_separation(self):
        well_separated = orientation_sampling.sphere_coverage_metrics(
            [[0.0, 0.0, 1.0], [0.0, 0.0, -1.0]],
            evaluation_points=1000,
        )
        close = orientation_sampling.sphere_coverage_metrics(
            [[0.0, 0.0, 1.0], [0.01, 0.0, 1.0]],
            evaluation_points=1000,
        )

        self.assertGreater(
            well_separated.minimum_separation_degrees,
            close.minimum_separation_degrees,
        )
        self.assertLess(
            well_separated.covering_radius_degrees,
            close.covering_radius_degrees,
        )

    def test_generate_directions_supports_candidate_methods(self):
        for method in (
            "random",
            "lhs",
            "lhs_maximin",
            "fibonacci",
            "fibonacci_maximin",
        ):
            points = orientation_sampling.generate_directions(method, 6, seed=2)
            self.assertEqual(points.shape, (6, 3))

    def test_halton_quaternions_are_deterministic_unit_quaternions(self):
        first = orientation_sampling.halton_quaternions(8)
        second = orientation_sampling.halton_quaternions(8)

        np.testing.assert_allclose(first, second)
        np.testing.assert_allclose(np.linalg.norm(first, axis=1), np.ones(8))

    def test_rotation_metrics_reward_distinct_orientations(self):
        well_separated = orientation_sampling.quaternion_coverage_metrics(
            orientation_sampling.halton_quaternions(4),
            evaluation_points=500,
        )
        close = orientation_sampling.quaternion_coverage_metrics(
            np.repeat([[0.0, 0.0, 0.0, 1.0]], 4, axis=0),
            evaluation_points=500,
        )

        self.assertGreater(
            well_separated.minimum_separation_degrees,
            close.minimum_separation_degrees,
        )
        self.assertLess(
            well_separated.covering_radius_degrees,
            close.covering_radius_degrees,
        )

    def test_generate_trial_vectors_uses_direction_and_rotation_defaults(self):
        trial_vectors = orientation_sampling.generate_trial_vectors(6)

        self.assertEqual(trial_vectors.shape, (6, 6))
        np.testing.assert_allclose(np.linalg.norm(trial_vectors[:, :3], axis=1), np.ones(6))
        self.assertFalse(np.allclose(trial_vectors[:, 3:], 0.0))

    def test_default_trial_vectors_preserve_full_sphere_fibonacci_design(self):
        expected = orientation_sampling.fibonacci_sphere(12)

        generated = orientation_sampling.generate_trial_vectors(12, use_angles=False)

        np.testing.assert_allclose(generated[:, :3], expected)

    def test_sequence_offset_generates_reproducible_distinct_variants(self):
        first = orientation_sampling.generate_trial_vectors(8, sequence_offset=2)
        second = orientation_sampling.generate_trial_vectors(8, sequence_offset=2)
        different = orientation_sampling.generate_trial_vectors(8, sequence_offset=3)

        np.testing.assert_allclose(first, second)
        self.assertFalse(np.allclose(first, different))

    def test_benchmark_cli_prints_quality_columns(self):
        stdout = io.StringIO()
        with mock.patch(
            "sys.argv",
            [
                "pyar-benchmark-orientations",
                "-N",
                "4",
                "-m",
                "fibonacci",
                "--evaluation-points",
                "100",
                "-r",
                "1",
            ],
        ):
            with mock.patch("sys.stdout", stdout):
                benchmark_orientations.main()

        output = stdout.getvalue()
        self.assertIn("min_sep_deg", output)
        self.assertIn("cover_radius_deg", output)
        self.assertIn("4\tfibonacci", output)


if __name__ == "__main__":
    unittest.main()
