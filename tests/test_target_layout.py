"""Import tests for the PyAR 2.0 public API."""

import unittest


class TargetLayoutImportTests(unittest.TestCase):
    def test_public_api_entrypoints_are_importable(self):
        import pyar
        from pyar import Molecule
        from pyar import backends, sampling
        from pyar.backends import get_backend_capabilities
        from pyar.biases import afir, softmin
        from pyar.io import AggregateResult, ConformerResult, ReactionResult, WorkflowResult
        from pyar.io import results
        from pyar.selection import choose_geometries, print_energy_table, read_energy_from_xyz_file, reports
        from pyar.aggregation.layout_v2 import AggregateRunStateV2 as AggregateRunState
        from pyar.state.reaction_v2 import ReactionRunStateV2 as ReactionRunState
        from pyar.aggregation.layout_v2 import AggregateRunStateV2 as AggregateRunStateImpl
        from pyar.state.reaction_v2 import ReactionRunStateV2 as ReactionRunStateImpl
        from pyar.workflows import aggregate, conformer_search, react, microsolvate

        self.assertEqual(pyar.__all__, ["Molecule", "__version__"])
        self.assertTrue(hasattr(Molecule, "from_xyz"))
        self.assertTrue(hasattr(results, "WorkflowResult"))
        self.assertTrue(issubclass(AggregateResult, WorkflowResult))
        self.assertTrue(issubclass(ConformerResult, WorkflowResult))
        self.assertTrue(issubclass(ReactionResult, WorkflowResult))
        self.assertTrue(hasattr(backends, "SF"))
        self.assertTrue(callable(get_backend_capabilities))
        self.assertTrue(hasattr(sampling, "generate_trial_vectors"))
        self.assertTrue(hasattr(sampling, "generate_directions"))
        self.assertTrue(hasattr(sampling, "fibonacci_sphere"))
        self.assertTrue(hasattr(sampling, "latin_hypercube_sphere"))
        self.assertTrue(hasattr(sampling, "random_sphere"))
        self.assertTrue(hasattr(sampling, "maximin_select"))
        self.assertTrue(hasattr(sampling, "halton_quaternions"))
        self.assertTrue(hasattr(sampling, "random_quaternions"))
        self.assertTrue(hasattr(sampling, "sphere_coverage_metrics"))
        self.assertTrue(hasattr(sampling, "quaternion_coverage_metrics"))
        self.assertTrue(callable(choose_geometries))
        self.assertTrue(callable(print_energy_table))
        self.assertTrue(callable(read_energy_from_xyz_file))
        self.assertTrue(callable(reports.print_energy_table))
        self.assertTrue(callable(aggregate))
        self.assertTrue(callable(conformer_search))
        self.assertTrue(callable(react))
        self.assertTrue(callable(microsolvate))
        self.assertTrue(callable(afir.isotropic))
        self.assertTrue(callable(softmin.softmin))
        self.assertIs(AggregateRunState, AggregateRunStateImpl)
        self.assertIs(ReactionRunState, ReactionRunStateImpl)

    def test_retired_import_aliases_are_absent(self):
        from importlib.util import find_spec
        for name in ('pyar.Molecule', 'pyar.molecule_geometry', 'pyar.molecule_io',
                     'pyar.orientation_sampling', 'pyar.trial_generation', 'pyar.aggregator',
                     'pyar.reactor', 'pyar.modern_cli'):
            with self.subTest(name=name):
                self.assertIsNone(find_spec(name))


if __name__ == "__main__":
    unittest.main()
