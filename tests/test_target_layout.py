"""Import tests for the PyAR 2.0 public API."""

import importlib
import unittest


class TargetLayoutImportTests(unittest.TestCase):
    def test_public_api_entrypoints_are_importable(self):
        import pyar
        from pyar import Molecule
        from pyar import backends, sampling
        from pyar.backends import get_backend_capabilities
        from pyar.biases import afir
        from pyar.io import AggregateResult, ReactionResult, SolvationResult, WorkflowResult
        from pyar.io import results
        from pyar.selection import choose_geometries
        from pyar.state import AggregateRunState, ReactionRunState, SolvationRunState
        from pyar.state.aggregate import AggregateRunState as AggregateRunStateImpl
        from pyar.state.reaction import ReactionRunState as ReactionRunStateImpl
        from pyar.state.solvation import SolvationRunState as SolvationRunStateImpl
        from pyar.workflows import aggregate, react, solvate

        self.assertEqual(pyar.__all__, ["Molecule", "__version__"])
        self.assertTrue(hasattr(Molecule, "from_xyz"))
        self.assertTrue(hasattr(results, "WorkflowResult"))
        self.assertTrue(issubclass(AggregateResult, WorkflowResult))
        self.assertTrue(issubclass(SolvationResult, WorkflowResult))
        self.assertTrue(issubclass(ReactionResult, WorkflowResult))
        self.assertTrue(hasattr(backends, "SF"))
        self.assertTrue(callable(get_backend_capabilities))
        self.assertTrue(hasattr(sampling, "generate_trial_vectors"))
        self.assertFalse(hasattr(sampling, "fibonacci_sphere"))
        self.assertTrue(callable(choose_geometries))
        self.assertTrue(callable(aggregate))
        self.assertTrue(callable(react))
        self.assertTrue(callable(solvate))
        self.assertTrue(callable(afir.isotropic))
        self.assertIs(AggregateRunState, AggregateRunStateImpl)
        self.assertIs(ReactionRunState, ReactionRunStateImpl)
        self.assertIs(SolvationRunState, SolvationRunStateImpl)

    def test_legacy_import_paths_are_removed(self):
        legacy_paths = (
            "pyar.Molecule",
            "pyar.orientation_sampling",
            "pyar.trial_generation",
            "pyar.aggregator",
            "pyar.reactor",
            "pyar.interface",
            "pyar.afir.restraints",
            "pyar.data_analysis.clustering",
        )

        for legacy_path in legacy_paths:
            with self.subTest(legacy_path=legacy_path):
                with self.assertRaises(ModuleNotFoundError):
                    importlib.import_module(legacy_path)


if __name__ == "__main__":
    unittest.main()
