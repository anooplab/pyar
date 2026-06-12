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
        from pyar.io import AggregateResult, ConformerResult, ReactionResult, SolvationResult, WorkflowResult
        from pyar.io import results
        from pyar.selection import choose_geometries, print_energy_table, read_energy_from_xyz_file, reports
        from pyar.state import AggregateRunState, ReactionRunState, SolvationRunState
        from pyar.state.aggregate import AggregateRunState as AggregateRunStateImpl
        from pyar.state.reaction import ReactionRunState as ReactionRunStateImpl
        from pyar.state.solvation import SolvationRunState as SolvationRunStateImpl
        from pyar.workflows import aggregate, conformer_search, react, solvate

        self.assertEqual(pyar.__all__, ["Molecule", "__version__"])
        self.assertTrue(hasattr(Molecule, "from_xyz"))
        self.assertTrue(hasattr(results, "WorkflowResult"))
        self.assertTrue(issubclass(AggregateResult, WorkflowResult))
        self.assertTrue(issubclass(ConformerResult, WorkflowResult))
        self.assertTrue(issubclass(SolvationResult, WorkflowResult))
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
        self.assertTrue(callable(solvate))
        self.assertTrue(callable(afir.isotropic))
        self.assertIs(AggregateRunState, AggregateRunStateImpl)
        self.assertIs(ReactionRunState, ReactionRunStateImpl)
        self.assertIs(SolvationRunState, SolvationRunStateImpl)

    def test_legacy_import_paths_forward_to_canonical_modules(self):
        legacy_to_canonical = {
            "pyar.Molecule": {
                "pyar.core.molecule": ["Molecule", "parse_xyz", "read_xyz"],
            },
            "pyar.molecule_geometry": {
                "pyar.core.geometry": [
                    "align",
                    "moments_of_inertia_tensor",
                    "move_to_centre_of_mass",
                    "move_to_origin",
                    "rotate_3d",
                    "rotation_matrix",
                    "translate",
                ],
            },
            "pyar.molecule_io": {
                "pyar.io.xyz": [
                    "XYZParseError",
                    "as_coordinates_array",
                    "parse_xyz",
                    "read_xyz",
                    "validate_fragments",
                    "write_turbomole_coord",
                    "write_xyz",
                ],
            },
            "pyar.orientation_sampling": {
                "pyar.sampling.sphere": [
                    "fibonacci_sphere",
                    "generate_directions",
                    "latin_hypercube_sphere",
                    "maximin_select",
                    "random_sphere",
                ],
                "pyar.sampling.rotation": [
                    "canonicalize_quaternion",
                    "halton_quaternions",
                    "quaternions_to_euler_zxz",
                    "random_quaternions",
                ],
                "pyar.sampling.metrics": [
                    "RotationCoverageMetrics",
                    "SphereCoverageMetrics",
                    "quaternion_coverage_metrics",
                    "sphere_coverage_metrics",
                ],
            },
            "pyar.trial_generation": {
                "pyar.sampling.trial_generator": [
                    "broken",
                    "check_close_contact",
                    "create_composite_molecule",
                    "create_trial_geometries",
                    "ellipsoid_wall_potential",
                    "generate_points",
                    "generate_trial_vectors",
                    "get_connectivity",
                    "merge_two_molecules",
                    "minimum_separation",
                    "plot_points",
                    "polar_to_cartesian",
                    "spherical_to_cartesian",
                    "write_trial_vectors",
                ],
            },
            "pyar.aggregator": {
                "pyar.workflows.aggregate": [
                    "aggregate",
                    "aggregate_from_formulas",
                    "expand_formula_to_aggregate_inputs",
                    "generate_molecule_from_formula",
                ],
                "pyar.workflows.solvation": ["solvate"],
            },
            "pyar.reactor": {
                "pyar.workflows.reaction": [
                    "build_gamma_schedule",
                    "build_reaction_request",
                    "format_gamma_id",
                    "initialize_reaction_run",
                    "optimize_all",
                    "print_header",
                    "react",
                    "relax_without_afir_bias",
                    "saved_product_identities",
                    "with_gamma",
                    "without_afir_bias",
                ],
            },
            "pyar.interface": {
                "pyar.backends": ["SF", "which", "write_xyz", "require_executable"],
            },
            "pyar.interface.mlopt": {
                "pyar.backends.mlopt": ["main"],
            },
            "pyar.afir.restraints": {
                "pyar.biases.afir": [
                    "covalent_radii",
                    "get_covalent_radius",
                    "get_data_structure",
                    "isotropic",
                    "main",
                    "resolve_gamma",
                ],
            },
            "pyar.data_analysis.clustering": {
                "pyar.selection.clustering": [
                    "calc_fingerprint_distance",
                    "choose_geometries",
                    "cluster_logger",
                    "plot_energy_histogram",
                    "read_energy_from_xyz_file",
                    "print_energy_table",
                ],
            },
            "pyar.selection.reports": {
                "pyar.selection.reports": [
                    "plot_energy_histogram",
                    "print_energy_table",
                    "read_energy_from_xyz_file",
                ],
            },
        }

        for legacy_path, canonical_modules in legacy_to_canonical.items():
            with self.subTest(legacy_path=legacy_path):
                legacy_module = importlib.import_module(legacy_path)
                for canonical_path, symbols in canonical_modules.items():
                    canonical_module = importlib.import_module(canonical_path)
                    for symbol in symbols:
                        with self.subTest(canonical_path=canonical_path, symbol=symbol):
                            self.assertIs(
                                getattr(legacy_module, symbol),
                                getattr(canonical_module, symbol),
                            )


if __name__ == "__main__":
    unittest.main()
