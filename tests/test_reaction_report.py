import json
import hashlib

import numpy as np

from pyar.backends import write_xyz
from pyar.workflows.reaction_report import (
    assemble_validated_path,
    build_reaction_characterization_report,
    write_reaction_characterization_report,
)
from pyar.xyz_trajectory import parse_xyz_trajectory


def _endpoint_stage(reactant, product, *, assignment=(0, 1)):
    return {
        "endpoint_assignment": list(assignment), "endpoints_are_minima": True,
        "forward_frequency": {"evaluated_coordinates_angstrom": reactant,
                              "energy_hartree": -1.0, "minimum_confirmed": True},
        "backward_frequency": {"evaluated_coordinates_angstrom": product,
                               "energy_hartree": -0.7, "minimum_confirmed": True},
    }


def test_validated_path_is_assembled_in_assigned_chemical_direction(tmp_path):
    symbols = ["H", "H"]
    reactant = [[0., 0., 0.], [1., 0., 0.]]
    ts = [[0., 0., 0.], [1.2, 0., 0.]]
    product = [[0., 0., 0.], [1.5, 0., 0.]]
    execution = tmp_path / "execution"
    execution.mkdir()
    write_xyz(symbols, ts, execution / "frequency_geometry.xyz")
    def trajectory(name, rows):
        with (execution / name).open("w") as stream:
            for index, (coords, energy) in enumerate(rows):
                stream.write(f"2\nimage={index} energy_hartree={energy}\n")
                for symbol, xyz in zip(symbols, coords):
                    stream.write(f"{symbol} {xyz[0]} {xyz[1]} {xyz[2]}\n")
    trajectory("irc_forward_path.xyz", ((ts, -.8), (reactant, -1.0)))
    trajectory("irc_forward_relaxed_path.xyz", ((reactant, -1.0),))
    trajectory("irc_backward_path.xyz", ((ts, -.8), (product, -.7)))
    trajectory("irc_backward_relaxed_path.xyz", ((product, -.7),))
    endpoint_stage = _endpoint_stage(reactant, product)
    report = assemble_validated_path(execution, symbols, endpoint_stage, [0, 1],
                                     tmp_path / "validated_reaction_path.xyz")
    frames = parse_xyz_trajectory(tmp_path / "validated_reaction_path.xyz")
    assert report["frame_count"] == 3
    assert np.allclose(frames[0][1], reactant)
    assert np.allclose(frames[1][1], ts)
    assert np.allclose(frames[2][1], product)
    swapped = assemble_validated_path(execution, symbols, endpoint_stage, [1, 0],
                                      tmp_path / "swapped.xyz")
    swapped_frames = parse_xyz_trajectory(tmp_path / "swapped.xyz")
    assert swapped["reactant_branch"] == "backward"
    assert np.allclose(swapped_frames[0][1], product)
    assert np.allclose(swapped_frames[-1][1], reactant)


def test_report_validates_electronic_energy_arithmetic_and_separated_reference(tmp_path):
    symbols = ["H", "H"]
    reactant = [[0., 0., 0.], [1., 0., 0.]]
    ts = [[0., 0., 0.], [1.2, 0., 0.]]
    product = [[0., 0., 0.], [1.5, 0., 0.]]
    execution = tmp_path / "execution"
    execution.mkdir()
    write_xyz(symbols, ts, execution / "frequency_geometry.xyz")
    def trajectory(name, rows):
        with (execution / name).open("w") as stream:
            for coords, energy in rows:
                stream.write(f"2\nenergy_hartree={energy}\n")
                for symbol, xyz in zip(symbols, coords):
                    stream.write(f"{symbol} {xyz[0]} {xyz[1]} {xyz[2]}\n")
    trajectory("irc_forward_path.xyz", ((ts, -.8), (reactant, -1.0)))
    trajectory("irc_forward_relaxed_path.xyz", ((reactant, -1.0),))
    trajectory("irc_backward_path.xyz", ((ts, -.8), (product, -.7)))
    trajectory("irc_backward_relaxed_path.xyz", ((product, -.7),))
    with (execution / "irc_path.xyz").open("w") as stream:
        for coords, energy in ((reactant, -1.0), (ts, -.8), (product, -.7)):
            stream.write(f"2\nenergy_hartree={energy}\n")
            for symbol, xyz in zip(symbols, coords):
                stream.write(f"{symbol} {xyz[0]} {xyz[1]} {xyz[2]}\n")
    method = {"software": "xtb", "xtb_model": "gfn2"}
    stage_records = {}
    for name in ("relax", "neb", "ts", "frequency", "irc", "endpoint-relax", "endpoint-frequency"):
        artifact = execution / f"{name}.artifact"
        artifact.write_text(name)
        stage_records[name] = {"qc_params": method, "artifacts": {
            str(artifact): hashlib.sha256(artifact.read_bytes()).hexdigest()}}
    path_result = {
        "status": "complete", "execution_directory": str(execution),
        **stage_records,
        "frequency": {**stage_records["frequency"], "first_order_saddle_confirmed": True,
                      "imaginary_frequencies_cm-1": [-100.], "energy_hartree": -.8},
        "irc": {**stage_records["irc"], "irc_converged": True},
        "endpoint-frequency": {**stage_records["endpoint-frequency"], **_endpoint_stage(reactant, product)},
    }
    route = {"status": "complete", "through": "all", "reactant_product_connection_confirmed": True}
    provenance = {"product_id": "product_001", "pathway_id": "pathway_001"}
    report = build_reaction_characterization_report(
        route_root=tmp_path, product={"symbols": symbols}, provenance=provenance,
        path_result=path_result, route_result=route, alternatives=[],
        reference_mode="separated-reactants",
        reactant_reference={"reference_energy_hartree": -1.2,
                            "request": {"physical_method": {"software": "xtb", "xtb_model": "gfn2"}}},
    )
    energies = report["electronic_energetics"]
    assert energies["validated"]
    assert energies["E_R"] == -1.2
    assert energies["delta_E_rxn"] == .5
    assert np.isclose(energies["delta_E_forward_activation"], .4)
    assert np.isclose(energies["delta_E_reverse_activation"], -.1)
    assert report["validated_reaction_path"]["frame_count"] == 3
    output = write_reaction_characterization_report(tmp_path / "reaction_characterization.json", report)
    assert json.loads(output.read_text())["electronic_energetics"]["validated"]


def test_failed_connection_suppresses_energetics_and_validated_trace(tmp_path):
    method = {"software": "xtb"}
    path_result = {"status": "scientific_gate_failed", "execution_directory": str(tmp_path),
                   "frequency": {"qc_params": method, "first_order_saddle_confirmed": True,
                                 "energy_hartree": -.8},
                   "irc": {"qc_params": method, "irc_converged": True},
                   "endpoint-frequency": {"qc_params": method, "endpoints_are_minima": True}}
    report = build_reaction_characterization_report(
        route_root=tmp_path, product={}, provenance={}, path_result=path_result,
        route_result={"status": "failed", "reactant_product_connection_confirmed": False},
        alternatives=[],
    )
    assert report["electronic_energetics"]["validated"] is False
    assert report["electronic_energetics"]["E_R"] is None
    assert report["validated_reaction_path"] is None
    assert not (tmp_path / "validated_reaction_path.xyz").exists()
