import io
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import numpy as np

from pyar.sampling import (
    fibonacci_sphere,
    generate_directions,
    generate_trial_vectors,
    halton_quaternions,
    maximin_select,
    quaternion_coverage_metrics,
    random_quaternions,
    random_sphere,
    sphere_coverage_metrics,
)
from pyar.sampling.trial_generator import write_trial_vectors
from pyar.scripts import benchmark_orientations

orientation_sampling = SimpleNamespace(
    fibonacci_sphere=fibonacci_sphere,
    maximin_select=maximin_select,
    sphere_coverage_metrics=sphere_coverage_metrics,
    generate_directions=generate_directions,
    halton_quaternions=halton_quaternions,
    quaternion_coverage_metrics=quaternion_coverage_metrics,
    generate_trial_vectors=generate_trial_vectors,
    random_sphere=random_sphere,
    random_quaternions=random_quaternions,
)


class OrientationSamplingTests(unittest.TestCase):
    def test_fibonacci_points_are_deterministic_unit_vectors(self):
        first = orientation_sampling.fibonacci_sphere(8)
        second = orientation_sampling.fibonacci_sphere(8)

        np.testing.assert_allclose(first, second)
        np.testing.assert_allclose(np.linalg.norm(first, axis=1), np.ones(8))

    def test_generate_directions_fibonacci_is_deterministic_and_offset_sensitive(self):
        baseline_first = orientation_sampling.generate_directions("fibonacci", 8)
        baseline_second = orientation_sampling.generate_directions("fibonacci", 8)
        offset_first = orientation_sampling.generate_directions("fibonacci", 8, sequence_offset=2)
        offset_second = orientation_sampling.generate_directions("fibonacci", 8, sequence_offset=2)
        different_offset = orientation_sampling.generate_directions("fibonacci", 8, sequence_offset=3)

        np.testing.assert_allclose(baseline_first, baseline_second)
        np.testing.assert_allclose(offset_first, offset_second)
        self.assertFalse(np.allclose(baseline_first, different_offset))
        np.testing.assert_allclose(np.linalg.norm(baseline_first, axis=1), np.ones(8))
        np.testing.assert_allclose(np.linalg.norm(offset_first, axis=1), np.ones(8))

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
        for method in ("random", "lhs", "lhs_maximin", "fibonacci", "fibonacci_maximin"):
            points = orientation_sampling.generate_directions(method, 6, seed=2)
            self.assertEqual(points.shape, (6, 3))
            np.testing.assert_allclose(np.linalg.norm(points, axis=1), np.ones(6), atol=1e-12)

    def test_random_direction_and_rotation_sampling_are_seeded(self):
        first_directions = orientation_sampling.random_sphere(6, seed=17)
        second_directions = orientation_sampling.random_sphere(6, seed=17)
        first_quaternions = orientation_sampling.random_quaternions(6, seed=17)
        second_quaternions = orientation_sampling.random_quaternions(6, seed=17)

        np.testing.assert_allclose(first_directions, second_directions)
        np.testing.assert_allclose(first_quaternions, second_quaternions)
        np.testing.assert_allclose(np.linalg.norm(first_directions, axis=1), np.ones(6), atol=1e-12)
        np.testing.assert_allclose(np.linalg.norm(first_quaternions, axis=1), np.ones(6), atol=1e-12)

    def test_halton_quaternions_are_deterministic_unit_quaternions(self):
        first = orientation_sampling.halton_quaternions(8, seed=11)
        second = orientation_sampling.halton_quaternions(8, seed=11)
        different = orientation_sampling.halton_quaternions(8, seed=12)

        np.testing.assert_allclose(first, second)
        self.assertFalse(np.allclose(first, different))
        np.testing.assert_allclose(np.linalg.norm(first, axis=1), np.ones(8), atol=1e-12)

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

    def test_metrics_return_expected_fields_and_reject_invalid_inputs(self):
        sphere_metrics = orientation_sampling.sphere_coverage_metrics(
            [[0.0, 0.0, 1.0], [0.0, 0.0, -1.0]],
            evaluation_points=500,
        )
        quaternion_metrics = orientation_sampling.quaternion_coverage_metrics(
            orientation_sampling.halton_quaternions(4, seed=3),
            evaluation_points=500,
        )

        self.assertTrue(hasattr(sphere_metrics, "minimum_separation_degrees"))
        self.assertTrue(hasattr(sphere_metrics, "centroid_norm"))
        self.assertTrue(hasattr(quaternion_metrics, "covering_radius_degrees"))
        self.assertEqual(
            sphere_metrics,
            orientation_sampling.sphere_coverage_metrics(
                [[0.0, 0.0, 1.0], [0.0, 0.0, -1.0]],
                evaluation_points=500,
            ),
        )
        self.assertEqual(
            quaternion_metrics,
            orientation_sampling.quaternion_coverage_metrics(
                orientation_sampling.halton_quaternions(4, seed=3),
                evaluation_points=500,
            ),
        )

        with self.assertRaisesRegex(ValueError, "points must have shape"):
            orientation_sampling.sphere_coverage_metrics([], evaluation_points=10)
        with self.assertRaisesRegex(ValueError, "quaternions must have shape"):
            orientation_sampling.quaternion_coverage_metrics([], evaluation_points=10)

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
        self.assertFalse(np.allclose(first[:, :3], different[:, :3]))
        self.assertFalse(np.allclose(first[:, 3:], different[:, 3:]))

    def test_write_trial_vectors_uses_stable_overwrite_format(self):
        vectors = np.array(
            [
                [1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 1.5, -2.25, 3.75],
            ],
            dtype=float,
        )
        expected_text = (
            "1.0000000000000000e+00 0.0000000000000000e+00 0.0000000000000000e+00 "
            "0.0000000000000000e+00 0.0000000000000000e+00 0.0000000000000000e+00\n"
            "0.0000000000000000e+00 1.0000000000000000e+00 0.0000000000000000e+00 "
            "1.5000000000000000e+00 -2.2500000000000000e+00 3.7500000000000000e+00\n"
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "trial_vectors.dat"
            write_trial_vectors(vectors, output_path)
            first = output_path.read_text()
            write_trial_vectors(vectors, output_path)
            second = output_path.read_text()

        self.assertEqual(first, expected_text)
        self.assertEqual(second, expected_text)

    def test_generate_and_write_trial_vectors_are_reproducible(self):
        vectors = orientation_sampling.generate_trial_vectors(4, sequence_offset=5)

        with tempfile.TemporaryDirectory() as tmpdir:
            first_path = Path(tmpdir) / "first.dat"
            second_path = Path(tmpdir) / "second.dat"
            write_trial_vectors(vectors, first_path)
            write_trial_vectors(vectors, second_path)

            first_text = first_path.read_text()
            second_text = second_path.read_text()

        self.assertEqual(first_text, second_text)
        self.assertEqual(first_text.count("\n"), 4)

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

    def test_invalid_sampling_method_names_fail_clearly(self):
        with self.assertRaisesRegex(ValueError, "Unknown direction method: 'foo'"):
            orientation_sampling.generate_directions("foo", 4)
        with self.assertRaisesRegex(ValueError, "Valid methods are: halton, random"):
            orientation_sampling.generate_trial_vectors(4, rotation_method="foo")


if __name__ == "__main__":
    unittest.main()
