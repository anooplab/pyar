"""Tests for XYZ-only structural analysis and reporting."""

import io
import json
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

from pyar.structure_comparison.coordinate_graph import (
    analyze_coordinate_structure,
    analyze_growth_transition,
    infer_coordinate_graph,
)


def _molecule(symbols, coordinates, name="fixture"):
    return SimpleNamespace(atoms_list=symbols, coordinates=coordinates, name=name, fragments=[])


class CoordinateGraphTests(unittest.TestCase):
    def test_coordinate_graph_assigns_adjacency_only(self):
        molecule = _molecule(["C", "H"], [[0, 0, 0], [0.7, 0, 0]])

        graph = infer_coordinate_graph(molecule)
        report = analyze_coordinate_structure(molecule)

        self.assertEqual(graph.number_of_edges(), 1)
        self.assertIn("normalized_distance", graph.edges[0, 1])
        self.assertEqual(report["geometric_component_pattern"], "single-connected-component")
        self.assertIn("not a bond-order", report["interpretation"])
        self.assertNotIn("charge", report)

    def test_none_model_leaves_each_atom_as_a_component(self):
        molecule = _molecule(["H", "H"], [[0, 0, 0], [0.7, 0, 0]])

        report = analyze_coordinate_structure(molecule, model="none")

        self.assertEqual(report["edges"], [])
        self.assertEqual(report["component_count"], 2)

    def test_explicit_distance_cutoff_is_coordinate_only(self):
        molecule = _molecule(["Au", "Au"], [[0, 0, 0], [2.5, 0, 0]])

        report = analyze_coordinate_structure(
            molecule, model="distance-cutoff", cutoff=2.6
        )

        self.assertEqual(len(report["edges"]), 1)
        self.assertEqual(report["cutoff_angstrom"], 2.6)

    def test_growth_transition_reports_coordinate_edge_change(self):
        carbon = _molecule(["C"], [[0, 0, 0]], "carbon")
        hydrogen = _molecule(["H"], [[5, 0, 0]], "hydrogen")
        product = _molecule(["C", "H"], [[0, 0, 0], [0.7, 0, 0]], "CH")

        report = analyze_growth_transition([carbon, hydrogen], product)

        self.assertEqual(report["input_to_output"]["added_edges"], [[0, 1]])
        self.assertEqual(report["input_to_output"]["input_component_count"], 2)
        self.assertEqual(report["input_to_output"]["output_component_count"], 1)
        self.assertIn("not proof", report["limitations"][1])

    def test_cli_analyze_writes_json_without_chemical_identity_fields(self):
        from pyar.scripts import clustering as clustering_script

        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir, "fixture.xyz")
            path.write_text("2\nfixture: energy=-1.0\nC 0 0 0\nH 0.7 0 0\n")
            with mock.patch(
                "sys.argv", ["pyar-clustering", str(path), "--mode", "analyze"]
            ), redirect_stdout(io.StringIO()) as stdout:
                clustering_script.main()

        report = json.loads(stdout.getvalue())
        self.assertEqual(report["method"], "coordinate-only-adjacency")
        self.assertEqual(report["structures"][0]["component_count"], 1)
        self.assertNotIn("smiles", report["structures"][0])

    def test_cli_analyze_accepts_explicit_atomic_cluster_cutoff(self):
        from pyar.scripts import clustering as clustering_script

        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir, "cluster.xyz")
            path.write_text("2\ncluster: energy=-1.0\nAu 0 0 0\nAu 2.5 0 0\n")
            argv = [
                "pyar-clustering", str(path), "--mode", "analyze",
                "--coordinate-model", "distance-cutoff", "--bond-cutoff", "2.6",
            ]
            with mock.patch("sys.argv", argv), redirect_stdout(io.StringIO()) as stdout:
                clustering_script.main()

        report = json.loads(stdout.getvalue())
        self.assertEqual(report["structures"][0]["model"], "distance-cutoff")
        self.assertEqual(len(report["structures"][0]["edges"]), 1)


if __name__ == "__main__":
    unittest.main()
