"""Import tests for the target package layout scaffolding."""

import unittest


class TargetLayoutImportTests(unittest.TestCase):
    def test_core_and_io_entrypoints_are_importable(self):
        from pyar.core import Molecule
        from pyar.core import geometry
        from pyar.io import AggregateResult, ReactionResult, SolvationResult, WorkflowResult
        from pyar.io import results
        from pyar import backends
        from pyar import reactor as package_legacy_reaction_workflow
        from pyar.interface import xtb as legacy_xtb
        from pyar import sampling
        from pyar.state.aggregate import AggregateRunState as AggregateRunStateImpl
        from pyar.state.reaction import ReactionRunState as ReactionRunStateImpl
        from pyar.state.solvation import SolvationRunState as SolvationRunStateImpl
        from pyar.state import AggregateRunState, ReactionRunState, SolvationRunState
        from pyar.workflows.reaction import react
        from pyar.workflows import reaction as reaction_workflow
        import pyar.reactor as legacy_reaction_workflow

        self.assertTrue(hasattr(Molecule, "from_xyz"))
        self.assertTrue(hasattr(geometry, "rotate_3d"))
        self.assertTrue(hasattr(results, "WorkflowResult"))
        self.assertTrue(issubclass(AggregateResult, WorkflowResult))
        self.assertTrue(issubclass(SolvationResult, WorkflowResult))
        self.assertTrue(issubclass(ReactionResult, WorkflowResult))
        self.assertTrue(hasattr(backends, "SF"))
        from pyar.backends import geometric, xtb, orca

        self.assertTrue(hasattr(geometric, "Geometric"))
        self.assertTrue(hasattr(xtb, "Xtb"))
        self.assertTrue(hasattr(orca, "Orca"))
        self.assertIs(legacy_xtb, xtb)
        self.assertTrue(hasattr(sampling, "generate_trial_vectors"))
        self.assertIsNotNone(AggregateRunState)
        self.assertIsNotNone(ReactionRunState)
        self.assertIsNotNone(SolvationRunState)
        self.assertTrue(callable(react))
        self.assertIs(legacy_reaction_workflow, reaction_workflow)
        self.assertIs(package_legacy_reaction_workflow, reaction_workflow)
        self.assertIs(AggregateRunState, AggregateRunStateImpl)
        self.assertIs(ReactionRunState, ReactionRunStateImpl)
        self.assertIs(SolvationRunState, SolvationRunStateImpl)


if __name__ == "__main__":
    unittest.main()
