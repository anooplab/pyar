from pathlib import Path
import json
from types import SimpleNamespace

import numpy as np
import pytest
from ase.calculators.calculator import Calculator, all_changes

from pyar.neb import (
    _aligned_rmsd,
    _bond_set,
    _relax_endpoint,
    build_neb_images,
    read_xyz,
    _build_parser,
    _frequency,
    _load_stage,
    _match_endpoints,
    _save_stage,
    _write_xyz_trajectory,
    run_neb,
)


DATA = Path(__file__).parent / "data" / "neb"


def test_cli_exposes_individual_workflow_stages():
    parser = _build_parser()
    defaults = vars(parser.parse_args(["--software", "xtb"]))
    assert defaults["stage"] == "all"
    assert defaults["interpolation"] == "linear"
    assert defaults["idpp_fmax"] == 0.1
    assert defaults["idpp_steps"] == 100
    for stage in ("relax", "neb", "ts", "frequency", "irc", "endpoints"):
        args = vars(parser.parse_args(["--stage", stage, "--software", "xtb"]))
        assert args["stage"] == stage
    idpp = vars(parser.parse_args([
        "--software", "xtb", "--interpolation", "idpp",
        "--idpp-fmax", "0.02", "--idpp-steps", "25",
    ]))
    assert (idpp["interpolation"], idpp["idpp_fmax"], idpp["idpp_steps"]) == (
        "idpp", 0.02, 25,
    )


class HarmonicCalculator(Calculator):
    implemented_properties = ["energy", "forces"]

    def __init__(self, target):
        super().__init__()
        self.target = np.asarray(target, dtype=float)

    def calculate(self, atoms=None, properties=("energy", "forces"), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        displacement = atoms.get_positions() - self.target
        self.results["energy"] = 0.5 * float(np.sum(displacement**2))
        self.results["forces"] = -displacement


def test_neb_images_include_ts_guess_as_middle_waypoint():
    symbols, images = build_neb_images(
        DATA / "hcn.xyz", DATA / "hnc.xyz", DATA / "guess.xyz", 5
    )

    assert symbols == ["H", "C", "N"]
    assert len(images) == 5
    np.testing.assert_allclose(images[0], read_xyz(DATA / "hcn.xyz")[1])
    np.testing.assert_allclose(images[2], read_xyz(DATA / "guess.xyz")[1])
    np.testing.assert_allclose(images[-1], read_xyz(DATA / "hnc.xyz")[1])


def test_linear_initialization_is_the_existing_piecewise_cartesian_path():
    symbols, default_images = build_neb_images(
        DATA / "hcn.xyz", DATA / "hnc.xyz", DATA / "guess.xyz", 7,
    )
    explicit_symbols, explicit_images = build_neb_images(
        DATA / "hcn.xyz", DATA / "hnc.xyz", DATA / "guess.xyz", 7,
        interpolation="linear",
    )
    start = read_xyz(DATA / "hcn.xyz")[1]
    end = read_xyz(DATA / "hnc.xyz")[1]
    ts = read_xyz(DATA / "guess.xyz")[1]
    expected = [
        ((1 - i / 3) * start + (i / 3) * ts) if i <= 3
        else ((1 - (i - 3) / 3) * ts + ((i - 3) / 3) * end)
        for i in range(7)
    ]
    assert symbols == explicit_symbols
    for actual, explicit, old_path in zip(default_images, explicit_images, expected):
        np.testing.assert_array_equal(actual, explicit)
        np.testing.assert_allclose(actual, old_path, rtol=0, atol=0)


def _write_xyz(path, symbols, coordinates):
    path.write_text(
        f"{len(symbols)}\nfixture\n" + "".join(
            f"{symbol} {x:.12f} {y:.12f} {z:.12f}\n"
            for symbol, (x, y, z) in zip(symbols, coordinates)
        )
    )


def test_idpp_initialization_preserves_waypoints_and_regularizes_fixture(tmp_path):
    symbols = ["H", "H", "H"]
    start = np.asarray([[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 2.0, 0.0]])
    ts = np.asarray([[0.8, 0.0, 0.0], [-0.8, 0.0, 0.0], [0.2, 1.5, 0.5]])
    end = np.asarray([[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 1.5, 1.0]])
    paths = [tmp_path / f"{name}.xyz" for name in ("start", "end", "ts")]
    for path, coordinates in zip(paths, (start, end, ts)):
        _write_xyz(path, symbols, coordinates)

    _, linear = build_neb_images(*paths, 5)
    returned_symbols, idpp = build_neb_images(
        *paths, 5, interpolation="idpp", idpp_fmax=0.01, idpp_steps=500,
    )

    assert returned_symbols == symbols
    assert len(idpp) == 5
    np.testing.assert_array_equal(idpp[0], start)
    np.testing.assert_array_equal(idpp[2], ts)
    np.testing.assert_array_equal(idpp[-1], end)
    assert all(np.all(np.isfinite(frame)) for frame in idpp)
    # The fixed atom mapping makes the linear midpoint an unphysical 0.2 A H-H contact.
    linear_contact = np.linalg.norm(linear[1][0] - linear[1][1])
    assert linear_contact == pytest.approx(0.2)
    idpp_contact = np.linalg.norm(idpp[1][0] - idpp[1][1])
    assert idpp_contact > linear_contact


def test_invalid_neb_interpolation_is_rejected(tmp_path):
    with pytest.raises(ValueError, match="expected 'linear' or 'idpp'"):
        build_neb_images(
            DATA / "hcn.xyz", DATA / "hnc.xyz", DATA / "guess.xyz", 5,
            interpolation="spline",
        )


def test_idpp_failure_is_reported_without_linear_fallback(tmp_path, monkeypatch):
    import ase.mep.neb

    def fail(*args, **kwargs):
        raise RuntimeError("fixture optimizer failure")

    monkeypatch.setattr(ase.mep.neb, "idpp_interpolate", fail)
    with pytest.raises(RuntimeError, match="ASE IDPP initialization failed.*fixture optimizer failure"):
        build_neb_images(
            DATA / "hcn.xyz", DATA / "hnc.xyz", DATA / "guess.xyz", 5,
            interpolation="idpp",
        )


def test_neb_requires_odd_image_count_and_matching_atom_order(tmp_path):
    with pytest.raises(ValueError, match="odd integer"):
        build_neb_images(DATA / "hcn.xyz", DATA / "hnc.xyz", DATA / "guess.xyz", 4)

    mismatch = tmp_path / "mismatch.xyz"
    mismatch.write_text("3\nwrong order\nC 0 0 0\nH 1 0 0\nN 2 0 0\n")
    with pytest.raises(ValueError, match="ordered elements"):
        build_neb_images(DATA / "hcn.xyz", mismatch, DATA / "guess.xyz", 5)


def test_hnc_product_has_expected_connectivity():
    symbols, coordinates, _ = read_xyz(DATA / "hnc.xyz")

    assert _bond_set(symbols, coordinates) == {(1, 2), (0, 2)}


def test_endpoint_relaxation_rejects_a_product_that_loses_its_bond(tmp_path):
    with pytest.raises(RuntimeError, match="did not survive unbiased relaxation"):
        _relax_endpoint(
            ["H", "H"],
            np.asarray([[0.0, 0.0, 0.0], [0.7, 0.0, 0.0]]),
            HarmonicCalculator([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]),
            tmp_path,
            "product",
            0.01,
            200,
        )


def test_aligned_rmsd_ignores_rigid_rotation_and_translation():
    geometry = np.asarray([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.2, 0.8, 0.0]])
    rotation = np.asarray([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    transformed = geometry @ rotation + np.asarray([4.0, -2.0, 1.0])

    assert _aligned_rmsd(geometry, transformed) < 1e-12


@pytest.mark.parametrize("contents", ["0\nempty\n", "-1\ninvalid\n", "1\nnan\nH nan 0 0\n", "1\nelement\nZz 0 0 0\n"])
def test_xyz_rejects_empty_nonfinite_and_invalid_elements(tmp_path, contents):
    path = tmp_path / "bad.xyz"
    path.write_text(contents)
    with pytest.raises(ValueError, match="Invalid XYZ"):
        read_xyz(path)


def test_two_irc_branches_cannot_both_match_the_same_reference():
    symbols, hcn, _ = read_xyz(DATA / "hcn.xyz")
    _, hnc, _ = read_xyz(DATA / "hnc.xyz")
    result = _match_endpoints(symbols, [hcn, hcn], [hcn, hnc], 0.5)
    assert not result["irc_endpoint_connectivities_match"]
    assert not result["irc_endpoint_geometries_match"]
    reverse = _match_endpoints(symbols, [hnc, hcn], [hcn, hnc], 0.5)
    assert reverse["endpoint_assignment"] == [1, 0]


def test_identical_observed_minima_fail_with_overlapping_reference_tolerances():
    symbols = ["C"] * 8
    references = [np.array([[i * spacing, 0., 0.] for i in range(8)])
                  for spacing in (1.4, 1.7)]
    shared_minimum = (references[0] + references[1]) / 2
    assert _aligned_rmsd(*references) > 0.5
    assert _bond_set(symbols, references[0]) == _bond_set(symbols, references[1])
    result = _match_endpoints(symbols, [shared_minimum, shared_minimum], references, 0.5)
    assert not result["observed_endpoints_distinct"]
    assert not result["irc_endpoint_geometries_match"]
    assert result["endpoint_assignment"] is None


def test_stage_handoff_rejects_modified_artifacts_and_different_methods(tmp_path):
    calculator = SimpleNamespace(qc_params={"software": "xtb", "charge": 0, "nprocs": 1})
    geometry = tmp_path / "geometry.xyz"
    geometry.write_text((DATA / "hcn.xyz").read_text())
    _save_stage(tmp_path, "ts", {"ts_optimization_converged": True}, calculator, [], [geometry])
    parallel = SimpleNamespace(qc_params=dict(calculator.qc_params, nprocs=8))
    assert _load_stage(tmp_path, "ts", parallel)["ts_optimization_converged"]
    charged = SimpleNamespace(qc_params=dict(calculator.qc_params, charge=1))
    with pytest.raises(ValueError, match="different backend"):
        _load_stage(tmp_path, "ts", charged)
    geometry.write_text((DATA / "hnc.xyz").read_text())
    with pytest.raises(ValueError, match="Stale ts artifact"):
        _load_stage(tmp_path, "ts", calculator)


@pytest.fixture
def fake_path_backend(monkeypatch):
    """Replace only expensive numerical stages; exercise real routing and gates."""
    import geometric.neb
    import pyar.neb as neb

    calls = []
    failures = set()
    symbols, hcn, _ = read_xyz(DATA / "hcn.xyz")
    _, hnc, _ = read_xyz(DATA / "hnc.xyz")
    _, guess, _ = read_xyz(DATA / "guess.xyz")
    monkeypatch.setattr("pyar.backends.geometric.PyarGeometricCalculator",
                        lambda qc_params: SimpleNamespace(qc_params=qc_params))

    def optimize(symbols, coordinates, calculator, output, label, max_cycles, **kwargs):
        calls.append(label)
        final = {"ts": guess, "irc_forward": hcn, "irc_backward": hnc}.get(label, coordinates)
        if "same_endpoint" in failures and label == "irc_backward":
            final = hcn
        frames, energies = [np.asarray(coordinates), np.asarray(final)], [-1., -2.]
        _write_xyz_trajectory(output / f"{label}_path.xyz", symbols, frames, energies)
        return frames, energies, label not in failures

    def frequency(symbols, coordinates, calculator, output, label, threshold):
        calls.append(f"frequency_{label}")
        np.savetxt(output / f"{label}_hessian.txt", np.eye(3 * len(symbols)))
        (output / f"{label}_frequencies.vdata").write_text("frequency data")
        return {"first_order_saddle_confirmed": label == "ts" and "frequency_ts" not in failures,
                "minimum_confirmed": label != "ts" and f"frequency_{label}" not in failures,
                "evaluated_coordinates_angstrom": np.asarray(coordinates).tolist(), "energy_hartree": -2.}

    def chain(molecule, engine, scratch, params, plain):
        calls.append("neb")
        structures = [SimpleNamespace(M=SimpleNamespace(xyzs=[frame]), energy=float(i == len(molecule.xyzs)//2))
                      for i, frame in enumerate(molecule.xyzs)]
        if "no_interior_maximum" in failures:
            structures[0].energy = 2.
        return SimpleNamespace(Structures=structures, avgg=0., maxg=1. if "neb" in failures else 0.)

    monkeypatch.setattr(neb, "_optimize_geometry", optimize)
    monkeypatch.setattr(neb, "_frequency", frequency)
    monkeypatch.setattr(geometric.neb, "ElasticBand", chain)
    monkeypatch.setattr(geometric.neb, "OptimizeChain", lambda chain, engine, params: (chain, 3))
    return calls, failures


def test_complete_workflow_includes_final_endpoint_optimization_and_frequencies(tmp_path, fake_path_backend):
    result = run_neb(DATA / "hcn.xyz", DATA / "hnc.xyz", DATA / "guess.xyz",
                     software="xtb", output=tmp_path)
    calls, _ = fake_path_backend
    assert calls == ["reactant", "product", "neb", "ts", "frequency_ts", "irc_forward", "irc_backward",
                     "irc_forward_relaxed", "frequency_irc_forward_relaxed",
                     "irc_backward_relaxed", "frequency_irc_backward_relaxed"]
    assert result["reactant_product_connection_confirmed"]
    assert result["endpoints"]["endpoints_are_minima"]
    assert json.loads((tmp_path / "workflow_summary.json").read_text())["status"] == "complete"


@pytest.mark.parametrize("failure", ["irc_forward", "same_endpoint", "irc_backward_relaxed", "frequency_irc_forward_relaxed"])
def test_endpoint_confirmation_requires_both_irc_branches_and_two_verified_minima(tmp_path, fake_path_backend, failure):
    _, failures = fake_path_backend
    failures.add(failure)
    result = run_neb(DATA / "hcn.xyz", DATA / "hnc.xyz", DATA / "guess.xyz",
                     software="xtb", output=tmp_path)
    assert not result["reactant_product_connection_confirmed"]


def test_independent_stages_do_not_repeat_previous_calculations(tmp_path, fake_path_backend):
    calls, _ = fake_path_backend
    run_neb(DATA / "hcn.xyz", DATA / "hnc.xyz", software="xtb", stage="relax", output=tmp_path)
    assert calls == ["reactant", "product"]
    run_neb(ts_guess=DATA / "guess.xyz", software="xtb", stage="neb", output=tmp_path)
    assert calls[-1] == "neb"
    run_neb(software="xtb", stage="ts", output=tmp_path)
    assert calls[-1] == "ts"
    assert "irc_forward" not in calls
    run_neb(software="xtb", stage="frequency", output=tmp_path)
    assert calls[-1] == "frequency_ts"
    run_neb(software="xtb", stage="irc", output=tmp_path)
    assert calls[-2:] == ["irc_forward", "irc_backward"]
    run_neb(software="xtb", stage="endpoints", output=tmp_path)
    assert len(calls) == 11


@pytest.mark.parametrize(
    ("stored", "requested"),
    [
        ({"interpolation": "linear"}, {"interpolation": "idpp"}),
        ({"interpolation": "idpp", "idpp_fmax": 0.1}, {"interpolation": "idpp", "idpp_fmax": 0.2}),
        ({"interpolation": "idpp", "idpp_steps": 100}, {"interpolation": "idpp", "idpp_steps": 101}),
    ],
)
def test_neb_initialization_options_are_restart_parameters(tmp_path, fake_path_backend, stored, requested):
    run_neb(DATA / "hcn.xyz", DATA / "hnc.xyz", software="xtb",
            stage="relax", output=tmp_path)
    stored_options = {"interpolation": "linear", "idpp_fmax": 0.1, "idpp_steps": 100}
    stored_options.update(stored)
    run_neb(ts_guess=DATA / "guess.xyz", software="xtb", stage="neb",
            images=3, output=tmp_path, **stored_options)
    parameters = json.loads((tmp_path / "neb_summary.json").read_text())["parameters"]
    assert parameters["interpolation"] == stored_options["interpolation"]
    requested_options = dict(stored_options, **requested)
    with pytest.raises(ValueError, match="Stage neb used different stage-specific parameters"):
        run_neb(software="xtb", stage="ts", output=tmp_path, **requested_options)


def test_pre_interpolation_schema_two_neb_summary_reuses_historical_linear_default(
    tmp_path, fake_path_backend,
):
    run_neb(DATA / "hcn.xyz", DATA / "hnc.xyz", software="xtb",
            stage="relax", output=tmp_path)
    run_neb(ts_guess=DATA / "guess.xyz", software="xtb", stage="neb",
            images=3, output=tmp_path)
    summary_path = tmp_path / "neb_summary.json"
    summary = json.loads(summary_path.read_text())
    for key in ("interpolation", "idpp_fmax", "idpp_steps"):
        summary["parameters"].pop(key)
    summary_path.write_text(json.dumps(summary))

    result = run_neb(software="xtb", stage="ts", images=3, output=tmp_path)
    assert result["ts_optimization_converged"]


def test_neb_ts_stage_consumes_optimized_highest_energy_image(tmp_path, fake_path_backend, monkeypatch):
    import pyar.neb as neb

    original_optimize = neb._optimize_geometry
    ts_input = {}

    def record_ts_input(symbols, coordinates, calculator, output, label, max_cycles, **kwargs):
        if label == "ts":
            ts_input["coordinates"] = np.asarray(coordinates).copy()
        return original_optimize(symbols, coordinates, calculator, output, label, max_cycles, **kwargs)

    monkeypatch.setattr(neb, "_optimize_geometry", record_ts_input)
    run_neb(DATA / "hcn.xyz", DATA / "hnc.xyz", DATA / "guess.xyz",
            software="xtb", output=tmp_path)
    assert "coordinates" in ts_input
    np.testing.assert_array_equal(
        ts_input["coordinates"], read_xyz(tmp_path / "ts_guess.xyz")[1],
    )


def test_standalone_ts_requires_no_endpoint_files(tmp_path, fake_path_backend):
    calls, _ = fake_path_backend
    result = run_neb(ts_geometry=DATA / "guess.xyz", software="xtb", stage="ts", output=tmp_path)
    assert result["ts_optimization_converged"]
    assert calls == ["ts"]


def test_failed_frequency_gate_prevents_irc_and_preserves_failed_status(tmp_path, fake_path_backend):
    calls, failures = fake_path_backend
    failures.add("frequency_ts")
    previous = Path.cwd()
    with pytest.raises(ValueError, match="stationary TS"):
        run_neb(DATA / "hcn.xyz", DATA / "hnc.xyz", DATA / "guess.xyz", software="xtb", output=tmp_path)
    assert Path.cwd() == previous
    assert "irc_forward" not in calls
    assert json.loads((tmp_path / "irc_summary.json").read_text())["status"] == "failed"


@pytest.mark.parametrize("failure", ["neb", "no_interior_maximum"])
def test_invalid_neb_does_not_launch_ts(tmp_path, fake_path_backend, failure):
    calls, failures = fake_path_backend
    failures.add(failure)
    with pytest.raises(ValueError, match="converged NEB with an interior maximum"):
        run_neb(DATA / "hcn.xyz", DATA / "hnc.xyz", DATA / "guess.xyz", software="xtb", output=tmp_path)
    assert "ts" not in calls


def test_upstream_rerun_invalidates_downstream_hessian(tmp_path, fake_path_backend):
    run_neb(DATA / "hcn.xyz", DATA / "hnc.xyz", DATA / "guess.xyz", software="xtb", output=tmp_path)
    (tmp_path / "ts_optimized.xyz").write_text((DATA / "hcn.xyz").read_text())
    with pytest.raises(ValueError, match="Stale frequency artifact"):
        run_neb(software="xtb", stage="irc", output=tmp_path)


@pytest.mark.parametrize("fail_validation", [False, True])
def test_frequency_summary_change_invalidates_irc_with_identical_artifacts(
        tmp_path, fake_path_backend, fail_validation):
    run_neb(DATA / "hcn.xyz", DATA / "hnc.xyz", DATA / "guess.xyz",
            software="xtb", output=tmp_path)
    original = json.loads((tmp_path / "frequency_summary.json").read_text())
    if fail_validation:
        fake_path_backend[1].add("frequency_ts")
    run_neb(software="xtb", stage="frequency", output=tmp_path,
            imaginary_frequency_threshold=40.)
    current = json.loads((tmp_path / "frequency_summary.json").read_text())
    assert original["artifacts"] == current["artifacts"]
    assert current["parameters"]["imaginary_frequency_threshold"] == 40.
    with pytest.raises(ValueError, match="Stale irc dependency: frequency"):
        run_neb(software="xtb", stage="endpoints", output=tmp_path)


def test_dependency_scientific_gate_is_checked_even_with_matching_hash(tmp_path):
    calculator = SimpleNamespace(qc_params={"software": "xtb"})
    _save_stage(tmp_path, "frequency", {"first_order_saddle_confirmed": False},
                calculator, [], [])
    _save_stage(tmp_path, "irc", {"irc_converged": True}, calculator, [], [], ["frequency"])
    with pytest.raises(ValueError, match="scientifically invalid frequency"):
        _load_stage(tmp_path, "irc", calculator)


def test_frequency_index_alone_does_not_establish_stationarity(tmp_path):
    symbols, coordinates, _ = read_xyz(DATA / "guess.xyz")
    calculator = HarmonicCalculator(coordinates + [0.1, 0., 0.])
    result = _frequency(symbols, coordinates, calculator, tmp_path, "test", 20.)
    assert result["imaginary_frequency_count"] == 0
    assert not result["stationary"]
    assert not result["minimum_confirmed"]


def test_numerically_linear_geometry_retains_both_bending_modes(tmp_path):
    symbols, geometry, _ = read_xyz(DATA / "hcn.xyz")
    distorted = geometry.copy()
    distorted[0, 1] += 1e-6
    result = _frequency(symbols, distorted, HarmonicCalculator(geometry), tmp_path, "test", 20.)
    assert 0 < result["linear_geometry_correction_angstrom"] < 1e-5
    assert len(result["frequencies_cm-1"]) == 4
