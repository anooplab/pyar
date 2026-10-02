from pathlib import Path
from importlib.util import find_spec
import json
import subprocess
import sys
from types import SimpleNamespace

import numpy as np
import pytest
from ase.calculators.calculator import Calculator, all_changes

from pyar.neb import (
    _aligned_rmsd,
    _align_rigid_frame,
    _ase_idpp_api,
    _bond_set,
    canonical_neb_parameters,
    canonical_ts_parameters,
    _relax_endpoint,
    build_neb_images,
    read_xyz,
    _build_parser,
    _frequency,
    _load_stage,
    _match_endpoints,
    _physical_settings,
    _select_neb_ts_candidate,
    _save_stage,
    _write_xyz_trajectory,
    run_neb,
)


DATA = Path(__file__).parent / "data" / "neb"


def test_cli_exposes_individual_workflow_stages():
    parser = _build_parser()
    defaults = vars(parser.parse_args(["--software", "xtb"]))
    assert defaults["stage"] == "all"
    assert defaults["xtb_model"] == "gxtb"
    assert defaults["interpolation"] == "linear"
    assert defaults["idpp_fmax"] == 0.1
    assert defaults["idpp_steps"] == 100
    assert defaults["geodesic_tol"] == 0.002
    assert defaults["geodesic_max_iter"] == 15
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
    geodesic = vars(parser.parse_args([
        "--software", "xtb", "--interpolation", "geodesic",
        "--geodesic-tol", "0.005", "--geodesic-max-iter", "8",
    ]))
    assert (geodesic["interpolation"], geodesic["geodesic_tol"],
            geodesic["geodesic_max_iter"]) == ("geodesic", 0.005, 8)
    assert defaults["ts_optimizer"] == "geometric"
    assert defaults["ts_fmax"] is None
    assert defaults["sella_fmax"] == 0.05
    shared_fmax = vars(parser.parse_args(["--software", "xtb", "--ts-fmax", "0.02"]))
    assert shared_fmax["ts_fmax"] == 0.02
    sella = vars(parser.parse_args([
        "--software", "xtb", "--ts-optimizer", "sella", "--sella-fmax", "0.03",
    ]))
    assert (sella["ts_optimizer"], sella["sella_fmax"]) == ("sella", 0.03)
    gfn2 = vars(parser.parse_args(["--software", "xtb", "--xtb-model", "gfn2"]))
    assert gfn2["xtb_model"] == "gfn2"
    with pytest.raises(SystemExit):
        parser.parse_args(["--software", "xtb", "--xtb-model", "gfn1"])


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


class AnalyticFirstOrderSaddleCalculator(Calculator):
    """Two-atom radial potential with a first-order saddle at a target bond length."""

    implemented_properties = ["energy", "forces"]

    def calculate(self, atoms=None, properties=("energy", "forces"), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        displacement = atoms.get_positions()[1] - atoms.get_positions()[0]
        distance = np.linalg.norm(displacement)
        delta = distance - 2.0
        derivative = -delta
        gradient = derivative * displacement / distance
        self.results["energy"] = -0.5 * delta**2
        self.results["forces"] = np.array([gradient, -gradient])


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
        idpp_fmax=float("nan"), idpp_steps=0,
        geodesic_tol=float("nan"), geodesic_max_iter=0,
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
        geodesic_tol=float("nan"), geodesic_max_iter=0,
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
    with pytest.raises(ValueError, match="expected 'linear', 'idpp', or 'geodesic'"):
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


def test_idpp_api_falls_back_to_legacy_ase_module_path(monkeypatch):
    import pyar.neb as neb

    interpolate = object()
    calls = []

    def import_module(name):
        calls.append(name)
        if name == "ase.mep.neb":
            raise ImportError("modern path unavailable")
        return SimpleNamespace(idpp_interpolate=interpolate)

    monkeypatch.setattr(neb.importlib, "import_module", import_module)
    atoms, actual_interpolate = _ase_idpp_api()

    from ase import Atoms

    assert atoms is Atoms
    assert actual_interpolate is interpolate
    assert calls == ["ase.mep.neb", "ase.neb"]


def test_idpp_api_unavailable_raises_clear_error(monkeypatch):
    import pyar.neb as neb

    def import_module(name):
        raise ImportError(f"{name} unavailable")

    monkeypatch.setattr(neb.importlib, "import_module", import_module)
    with pytest.raises(RuntimeError, match="does not provide IDPP interpolation"):
        _ase_idpp_api()


def test_geodesic_request_without_optional_dependency_is_clear_and_no_fallback(monkeypatch):
    import pyar.neb as neb

    monkeypatch.setitem(sys.modules, "geodesic_interpolate", None)
    with pytest.raises(RuntimeError, match=r"optional.*geodesic-interpolate.*pyar-chem\[geodesic\]"):
        neb._geodesic_api()


@pytest.mark.parametrize("stage", ["all", "neb"])
def test_geodesic_dependency_is_required_for_neb_initialization(tmp_path, monkeypatch, stage):
    import pyar.neb as neb

    monkeypatch.setitem(sys.modules, "geodesic_interpolate", None)
    with pytest.raises(RuntimeError, match=r"optional.*geodesic-interpolate.*pyar-chem\[geodesic\]"):
        run_neb(
            DATA / "hcn.xyz", DATA / "hnc.xyz", DATA / "guess.xyz",
            software="xtb", stage=stage, interpolation="geodesic", output=tmp_path,
        )


@pytest.mark.parametrize("stage", ["ts", "frequency"])
def test_standalone_geometry_stages_do_not_load_geodesic_dependency(
    tmp_path, fake_path_backend, monkeypatch, stage,
):
    import pyar.neb as neb

    def unexpected_geodesic_load():
        raise AssertionError("standalone stage must not load geodesic-interpolate")

    monkeypatch.setattr(neb, "_geodesic_api", unexpected_geodesic_load)
    monkeypatch.setattr(neb, "_installed_geodesic_version", lambda: None)
    result = run_neb(
        ts_geometry=DATA / "guess.xyz", software="xtb", stage=stage,
        interpolation="geodesic", output=tmp_path,
    )
    if stage == "ts":
        assert result["ts_optimization_converged"]
        assert fake_path_backend[0] == ["ts"]
    else:
        assert result["first_order_saddle_confirmed"]
        assert fake_path_backend[0] == ["frequency_ts"]
        assert result["backend_energy_gradient_evaluations"] == 2
        assert result["wall_seconds"] >= 0


@pytest.mark.skipif(
    find_spec("geodesic_interpolate") is None,
    reason="geodesic-interpolate is an optional dependency",
)
def test_completed_geodesic_neb_can_feed_downstream_stages_without_dependency(
    tmp_path, fake_path_backend, monkeypatch,
):
    import pyar.neb as neb

    run_neb(
        DATA / "hcn.xyz", DATA / "hnc.xyz", software="xtb",
        stage="relax", output=tmp_path,
    )
    run_neb(
        ts_guess=DATA / "guess.xyz", software="xtb", stage="neb", images=3,
        interpolation="geodesic",
        output=tmp_path,
    )
    neb_summary_path = tmp_path / "neb_summary.json"
    neb_summary = json.loads(neb_summary_path.read_text())
    neb_summary["geodesic_interpolate_version"] = "1.0.0"
    neb_summary_path.write_text(json.dumps(neb_summary))
    assert neb_summary["geodesic_interpolate_version"] == "1.0.0"

    def unexpected_geodesic_access(*args, **kwargs):
        raise AssertionError("consuming a completed stage must not inspect geodesic-interpolate")

    monkeypatch.setattr(neb, "_geodesic_api", unexpected_geodesic_access)
    monkeypatch.setattr(neb, "_installed_geodesic_version", unexpected_geodesic_access)
    monkeypatch.setattr(neb, "package_version", lambda package: "9.9.9")
    for stage in ("ts", "frequency", "irc", "endpoints"):
        result = run_neb(
            software="xtb", stage=stage, images=3,
            interpolation="geodesic", output=tmp_path,
        )
        assert result


def test_neb_import_and_help_do_not_load_optional_geodesic_package():
    code = """import sys
import pyar.neb
assert 'geodesic_interpolate' not in sys.modules
assert 'sella' not in sys.modules
from pyar.neb import main
try:
    main(['--help'])
except SystemExit as exc:
    assert exc.code == 0
assert 'geodesic_interpolate' not in sys.modules
assert 'sella' not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True, capture_output=True, text=True)


def test_geodesic_alignment_preserves_internal_geometry_without_reflection():
    reference = np.asarray([
        [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0],
    ])
    reflected = reference.copy()
    reflected[:, 0] *= -1.0
    restored = _align_rigid_frame(reflected, reference)
    before = np.linalg.norm(reflected[:, None, :] - reflected[None, :, :], axis=2)
    after = np.linalg.norm(restored[:, None, :] - restored[None, :, :], axis=2)
    np.testing.assert_allclose(after, before, atol=1e-12, rtol=0)
    assert not np.allclose(restored, reference, atol=1e-8, rtol=0)


@pytest.mark.skipif(
    find_spec("geodesic_interpolate") is None,
    reason="geodesic-interpolate is an optional dependency",
)
def test_geodesic_initialization_preserves_fixed_waypoints_frames_and_determinism(tmp_path):
    symbols, start, _ = read_xyz(DATA / "hcn.xyz")
    _, end, _ = read_xyz(DATA / "hnc.xyz")
    _, ts, _ = read_xyz(DATA / "guess.xyz")

    def rotate_translate(coords, angle, shift):
        rotation = np.asarray([
            [np.cos(angle), -np.sin(angle), 0.0],
            [np.sin(angle), np.cos(angle), 0.0],
            [0.0, 0.0, 1.0],
        ])
        return np.asarray(coords) @ rotation + np.asarray(shift)

    start = rotate_translate(start, 0.7, [8.0, -5.0, 3.0])
    ts = rotate_translate(ts, -1.1, [-4.0, 7.0, 1.5])
    end = rotate_translate(end, 2.0, [3.0, 6.0, -2.0])
    paths = [tmp_path / f"geo_{name}.xyz" for name in ("start", "end", "ts")]
    for path, coordinates in zip(paths, (start, end, ts)):
        _write_xyz(path, symbols, coordinates)
    start = read_xyz(paths[0])[1]
    end = read_xyz(paths[1])[1]
    ts = read_xyz(paths[2])[1]

    options = dict(
        interpolation="geodesic", idpp_fmax=0.02, idpp_steps=40,
        geodesic_tol=0.01, geodesic_max_iter=4,
    )
    returned_symbols, first = build_neb_images(*paths, 7, **options)
    _, second = build_neb_images(*paths, 7, **options)

    assert returned_symbols == symbols
    assert len(first) == 7
    np.testing.assert_array_equal(first[0], start)
    np.testing.assert_array_equal(first[3], ts)
    np.testing.assert_array_equal(first[-1], end)
    assert all(np.all(np.isfinite(frame)) for frame in first)
    for actual, repeated in zip(first, second):
        np.testing.assert_allclose(actual, repeated, atol=1e-10, rtol=0)


def test_neb_parameter_canonicalization_is_method_aware():
    common = dict(
        images=11, max_cycles=100, max_gradient=0.05, average_gradient=0.025,
        spring=1.0, climb=0.5, align=False,
    )
    linear_a = canonical_neb_parameters(dict(
        common, interpolation="linear", idpp_fmax=-1, idpp_steps=0,
        geodesic_tol=-1, geodesic_max_iter=0,
    ))
    linear_b = canonical_neb_parameters(dict(
        common, interpolation="linear", idpp_fmax=9, idpp_steps=80,
        geodesic_tol=8, geodesic_max_iter=99,
    ))
    assert linear_a == linear_b
    assert set(linear_a) == {*common, "interpolation"}

    idpp_a = canonical_neb_parameters(dict(
        common, interpolation="idpp", idpp_fmax=0.1, idpp_steps=100,
        geodesic_tol=0.002, geodesic_max_iter=15,
    ))
    idpp_b = canonical_neb_parameters(dict(
        common, interpolation="idpp", idpp_fmax=0.1, idpp_steps=100,
        geodesic_tol=9, geodesic_max_iter=1,
    ))
    assert idpp_a == idpp_b
    assert set(idpp_a) == {*common, "interpolation", "idpp_fmax", "idpp_steps"}

    geo_a = canonical_neb_parameters(dict(
        common, interpolation="geodesic", idpp_fmax=0.1, idpp_steps=100,
        geodesic_tol=0.002, geodesic_max_iter=15,
    ))
    geo_b = dict(geo_a, geodesic_tol=0.003)
    assert geo_a != geo_b


def test_ts_parameter_canonicalization_is_optimizer_aware():
    geometric_a = canonical_ts_parameters({
        "ts_max_cycles": 200, "ts_optimizer": "geometric", "sella_fmax": -1,
    })
    geometric_b = canonical_ts_parameters({
        "ts_max_cycles": 200, "ts_optimizer": "geometric", "sella_fmax": float("nan"),
    })
    assert geometric_a == geometric_b == {
        "ts_max_cycles": 200, "ts_optimizer": "geometric", "ts_fmax": None,
    }
    assert canonical_ts_parameters({"ts_max_cycles": 50, "ts_optimizer": "sella"}) == {
        "ts_max_cycles": 50, "ts_optimizer": "sella", "ts_fmax": None,
        "sella_fmax": 0.05,
    }
    assert canonical_ts_parameters({"ts_optimizer": "geometric", "ts_fmax": 0.02}) == {
        "ts_optimizer": "geometric", "ts_fmax": 0.02,
    }
    assert canonical_ts_parameters({"ts_optimizer": "sella", "ts_fmax": 0.02,
                                    "sella_fmax": 0.03}) == {
        "ts_optimizer": "sella", "ts_fmax": 0.02,
    }
    with pytest.raises(ValueError, match="ts_optimizer must"):
        canonical_ts_parameters({"ts_optimizer": "other"})


def test_sella_missing_dependency_is_clear(monkeypatch):
    import pyar.neb as neb

    monkeypatch.setitem(sys.modules, "sella", None)
    with pytest.raises(RuntimeError, match=r"optional 'sella'.*pyar-chem\[sella\]"):
        neb._sella_api()


def test_sella_request_fails_clearly_when_optional_package_is_missing(
    tmp_path, fake_path_backend, monkeypatch,
):
    import pyar.neb as neb

    monkeypatch.setitem(sys.modules, "sella", None)
    with pytest.raises(RuntimeError, match=r"optional 'sella'.*pyar-chem\[sella\]"):
        run_neb(
            ts_geometry=DATA / "guess.xyz", software="xtb", stage="ts",
            ts_optimizer="sella", output=tmp_path,
        )
    assert json.loads((tmp_path / "ts_summary.json").read_text())["status"] == "failed"


@pytest.mark.parametrize("failure", ["exception", "empty", "nonfinite_coordinates", "nonfinite_energy"])
def test_sella_optimizer_failure_paths_are_explicit(tmp_path, monkeypatch, failure):
    import pyar.neb as neb

    symbols, coordinates, _ = read_xyz(DATA / "guess.xyz")

    class FakeSella:
        def __init__(self, atoms, **kwargs):
            self.atoms = atoms

        def attach(self, observer, interval=1):
            self.observer = observer

        def run(self, **kwargs):
            if failure == "exception":
                raise RuntimeError("optimizer fixture failure")
            if failure == "empty":
                return False
            if failure == "nonfinite_coordinates":
                self.atoms.set_positions(np.full_like(self.atoms.positions, np.nan))
            self.observer()
            return False

    monkeypatch.setattr(neb, "_sella_api", lambda: (FakeSella, "fixture-version"))
    calculator = HarmonicCalculator(coordinates)
    if failure == "nonfinite_energy":
        original_calculate = calculator.calculate

        def nonfinite_energy(*args, **kwargs):
            original_calculate(*args, **kwargs)
            calculator.results["energy"] = float("nan")

        monkeypatch.setattr(calculator, "calculate", nonfinite_energy)
    expected = {
        "exception": (RuntimeError, "Sella TS optimization failed"),
        "empty": (RuntimeError, "empty TS optimization trajectory"),
        "nonfinite_coordinates": (ValueError, "non-finite coordinates"),
        "nonfinite_energy": (ValueError, "non-finite energy"),
    }[failure]
    with pytest.raises(expected[0], match=expected[1]):
        neb._optimize_sella(symbols, coordinates, calculator, 2, 0.05)


def test_sella_adapter_passes_ase_force_threshold_and_step_limit(monkeypatch):
    import pyar.neb as neb

    symbols, coordinates, _ = read_xyz(DATA / "guess.xyz")
    calls = {}

    class RecordingSella:
        def __init__(self, atoms, **kwargs):
            calls["constructor"] = kwargs
            self.nsteps = 4

        def attach(self, observer, interval=1):
            self.observer = observer

        def run(self, *, fmax, steps):
            calls.update(fmax=fmax, steps=steps)
            self.observer()
            return False

    monkeypatch.setattr(neb, "_sella_api", lambda: (RecordingSella, "test-version"))
    calculator = HarmonicCalculator(coordinates)
    frames, energies, converged, _ = neb._optimize_sella(
        symbols, coordinates, calculator, 23, 0.037,
    )
    assert calls["constructor"] == {"logfile": None, "order": 1, "internal": False}
    assert calls["fmax"] == 0.037
    assert calls["steps"] == 23
    assert converged is False
    assert calculator._pyar_optimizer_steps == 4
    assert len(frames) == len(energies) == 1


def test_sella_ts_stage_records_optimizer_only_and_downstream_reuses_artifacts(
    tmp_path, fake_path_backend, monkeypatch,
):
    import pyar.neb as neb

    symbols, coordinates, _ = read_xyz(DATA / "guess.xyz")
    observed = {}

    def optimize_sella(actual_symbols, actual_coordinates, calculator, steps, fmax):
        observed.update(steps=steps, fmax=fmax, calculator=calculator)
        final = np.asarray(actual_coordinates, dtype=float) + 0.01
        return [np.asarray(actual_coordinates), final], [-1.0, -2.0], False, "2.6.0"

    monkeypatch.setattr(neb, "_optimize_sella", optimize_sella)
    result = run_neb(
        ts_geometry=DATA / "guess.xyz", software="xtb", stage="ts", output=tmp_path,
        ts_optimizer="sella", sella_fmax=0.025, ts_max_cycles=17,
    )
    assert result["ts_optimization_converged"] is False
    assert result["sella_version"] == "2.6.0"
    assert "first_order_saddle_confirmed" not in result
    assert observed["steps"] == 17
    assert observed["fmax"] == 0.025
    written_symbols, written_coordinates, _ = read_xyz(tmp_path / "ts_optimized.xyz")
    assert written_symbols == symbols
    np.testing.assert_allclose(written_coordinates, coordinates + 0.01)
    path_symbols, path = _read_xyz_trajectory(tmp_path / "ts_path.xyz")
    assert path_symbols == symbols
    assert len(path) == 2
    np.testing.assert_allclose(path[0], coordinates)
    np.testing.assert_allclose(path[-1], written_coordinates)
    assert result["ts_energy_hartree"] == -2.0
    assert result["backend_energy_gradient_evaluations"] == 0
    assert result["optimizer_steps"] is None
    assert np.isfinite(result["wall_seconds"]) and result["wall_seconds"] >= 0
    assert result["effective_ts_convergence"] == {
        "max_steps": 17, "fmax_ev_per_angstrom": 0.025,
        "geometric_convergence_set": None,
    }

    def unexpected_sella_access(*args, **kwargs):
        raise AssertionError("downstream stages must not inspect Sella installation")

    monkeypatch.setattr(neb, "_sella_api", unexpected_sella_access)
    monkeypatch.setattr(neb, "package_version", unexpected_sella_access)
    with pytest.raises(ValueError, match="TS optimization did not converge"):
        run_neb(software="xtb", stage="frequency", output=tmp_path,
                ts_optimizer="geometric")


@pytest.mark.parametrize("optimizer", ["geometric", "sella"])
def test_shared_ts_fmax_is_passed_to_selected_optimizer_and_recorded(
    tmp_path, fake_path_backend, monkeypatch, optimizer,
):
    import pyar.neb as neb

    observed = {}
    coordinates = read_xyz(DATA / "guess.xyz")[1]

    if optimizer == "geometric":
        def optimize_geometry(symbols, xyz, calculator, output, label, max_cycles, **kwargs):
            observed.update(fmax=kwargs["fmax"], max_cycles=max_cycles)
            frame = np.asarray(xyz, dtype=float).copy()
            _write_xyz_trajectory(output / f"{label}_path.xyz", symbols, [frame], [-1.0])
            return [frame], [-1.0], True

        monkeypatch.setattr(neb, "_optimize_geometry", optimize_geometry)
    else:
        def optimize_sella(symbols, xyz, calculator, max_steps, fmax):
            observed.update(fmax=fmax, max_steps=max_steps)
            frame = np.asarray(xyz, dtype=float).copy()
            return [frame], [-1.0], True, "test-version"

        monkeypatch.setattr(neb, "_optimize_sella", optimize_sella)

    result = run_neb(
        ts_geometry=DATA / "guess.xyz", software="xtb", stage="ts", output=tmp_path,
        ts_optimizer=optimizer, ts_fmax=0.02, sella_fmax=0.04,
    )
    assert observed["fmax"] == 0.02
    assert result["effective_ts_convergence"]["fmax_ev_per_angstrom"] == 0.02
    assert result["parameters"]["ts_fmax"] == 0.02
    assert "sella_fmax" not in result["parameters"]
    _, optimized, _ = read_xyz(tmp_path / "ts_optimized.xyz")
    np.testing.assert_array_equal(optimized, coordinates)


def test_ts_fmax_is_part_of_restart_identity(tmp_path, fake_path_backend):
    run_neb(
        ts_geometry=DATA / "guess.xyz", software="xtb", stage="ts", output=tmp_path,
        ts_fmax=0.02,
    )
    summary = json.loads((tmp_path / "ts_summary.json").read_text())
    calculator = SimpleNamespace(qc_params=summary["qc_params"])
    same_parameters = canonical_ts_parameters({
        "ts_optimizer": "geometric", "ts_max_cycles": 200, "ts_fmax": 0.02,
        "sella_fmax": 0.05,
    })
    reused = _load_stage(
        tmp_path, "ts", calculator, expected_parameters=same_parameters,
    )
    assert reused["effective_ts_convergence"]["fmax_ev_per_angstrom"] == 0.02
    with pytest.raises(ValueError, match="different stage-specific parameters"):
        _load_stage(tmp_path, "ts", calculator, expected_parameters=canonical_ts_parameters({
            "ts_optimizer": "geometric", "ts_max_cycles": 200, "ts_fmax": 0.03,
            "sella_fmax": 0.05,
        }))


def test_inactive_sella_controls_do_not_validate_or_affect_geometric_ts(
    tmp_path, fake_path_backend,
):
    result = run_neb(
        ts_geometry=DATA / "guess.xyz", software="xtb", stage="ts", output=tmp_path,
        ts_optimizer="geometric", sella_fmax=float("nan"),
    )
    assert result["ts_optimization_converged"]
    assert result["optimizer_steps"] is None
    assert result["effective_ts_convergence"] == {
        "max_steps": 200, "fmax_ev_per_angstrom": None,
        "geometric_convergence_set": "GAU_TIGHT",
    }
    assert result["parameters"] == {
        "ts_max_cycles": 200, "ts_optimizer": "geometric", "ts_fmax": None,
    }

    with pytest.raises(ValueError, match="sella_fmax must be positive and finite"):
        run_neb(
            ts_geometry=DATA / "guess.xyz", software="xtb", stage="ts", output=tmp_path,
            ts_optimizer="sella", sella_fmax=float("nan"),
        )


def test_completed_sella_workflow_supports_all_downstream_stages_without_sella(
    tmp_path, fake_path_backend, monkeypatch,
):
    import pyar.neb as neb

    def optimize_sella(symbols, coordinates, calculator, max_steps, fmax):
        coordinates = np.asarray(coordinates, dtype=float).copy()
        return [coordinates], [-2.0], True, "2.6.0"

    monkeypatch.setattr(neb, "_optimize_sella", optimize_sella)
    result = run_neb(
        DATA / "hcn.xyz", DATA / "hnc.xyz", DATA / "guess.xyz", software="xtb",
        ts_optimizer="sella", output=tmp_path,
    )
    assert result["ts"]["sella_version"] == "2.6.0"

    def unexpected_sella_access(*args, **kwargs):
        raise AssertionError("downstream stages must not inspect Sella installation")

    monkeypatch.setattr(neb, "_sella_api", unexpected_sella_access)
    monkeypatch.setattr(neb, "package_version", unexpected_sella_access)
    for stage in ("frequency", "irc", "endpoints"):
        assert run_neb(software="xtb", stage=stage, output=tmp_path)


def _read_xyz_trajectory(path):
    """Small test reader for the multi-frame XYZ written by the TS stage."""
    lines = Path(path).read_text().splitlines()
    frames, symbols = [], None
    index = 0
    while index < len(lines):
        count = int(lines[index])
        frame_symbols, coordinates = [], []
        for line in lines[index + 2:index + 2 + count]:
            fields = line.split()
            frame_symbols.append(fields[0])
            coordinates.append([float(value) for value in fields[1:4]])
        symbols = frame_symbols if symbols is None else symbols
        frames.append(np.asarray(coordinates))
        index += count + 2
    return symbols, frames


@pytest.mark.skipif(find_spec("sella") is None, reason="Sella is an optional dependency")
def test_real_sella_package_runs_with_ase_calculator_and_returns_finite_trajectory(tmp_path):
    from importlib.metadata import version
    import pyar.neb as neb

    symbols, coordinates, _ = read_xyz(DATA / "guess.xyz")
    frames, energies, converged, installed_version = neb._optimize_sella(
        symbols, coordinates, HarmonicCalculator(coordinates), max_steps=1, fmax=1e6,
    )
    assert installed_version == version("sella")
    assert converged is True
    assert frames and len(frames) == len(energies)
    assert all(np.all(np.isfinite(frame)) for frame in frames)
    assert all(np.isfinite(energy) for energy in energies)
    np.testing.assert_allclose(frames[-1], coordinates, atol=1e-8)


@pytest.mark.skipif(find_spec("sella") is None, reason="Sella is an optional dependency")
def test_real_sella_refines_displaced_geometry_to_analytic_first_order_saddle():
    import pyar.neb as neb

    initial = np.array([[0.0, 0.0, 0.0], [2.4, 0.0, 0.0]])
    frames, energies, converged, _ = neb._optimize_sella(
        ["H", "H"], initial, AnalyticFirstOrderSaddleCalculator(), max_steps=100, fmax=1e-5,
    )

    assert converged
    assert len(frames) == len(energies) and len(frames) > 1
    final_distance = np.linalg.norm(frames[-1][1] - frames[-1][0])
    assert abs(final_distance - 2.0) < abs(2.4 - 2.0)
    assert abs(final_distance - 2.0) < 1e-4
    assert np.isfinite(energies[-1])


def test_linear_restart_ignores_inactive_interpolator_controls(tmp_path, fake_path_backend):
    run_neb(DATA / "hcn.xyz", DATA / "hnc.xyz", software="xtb",
            stage="relax", output=tmp_path)
    run_neb(ts_guess=DATA / "guess.xyz", software="xtb", stage="neb",
            images=3, interpolation="linear", output=tmp_path)
    stored = json.loads((tmp_path / "neb_summary.json").read_text())
    assert "idpp_fmax" not in stored["parameters"]
    assert "geodesic_tol" not in stored["parameters"]

    result = run_neb(software="xtb", stage="ts", images=3, interpolation="linear",
                     idpp_fmax=float("nan"), idpp_steps=0,
                     geodesic_tol=float("nan"), geodesic_max_iter=0, output=tmp_path)
    assert result["ts_optimization_converged"]


def test_pr21_flat_linear_summary_ignores_its_inactive_idpp_defaults(
    tmp_path, fake_path_backend,
):
    run_neb(DATA / "hcn.xyz", DATA / "hnc.xyz", software="xtb",
            stage="relax", output=tmp_path)
    run_neb(ts_guess=DATA / "guess.xyz", software="xtb", stage="neb",
            images=3, output=tmp_path)
    summary_path = tmp_path / "neb_summary.json"
    summary = json.loads(summary_path.read_text())
    summary["parameters"].update(
        interpolation="linear", idpp_fmax=0.1, idpp_steps=100,
    )
    summary_path.write_text(json.dumps(summary))

    result = run_neb(software="xtb", stage="ts", images=3, output=tmp_path)
    assert result["ts_optimization_converged"]


def test_neb_restart_rejects_unrecognized_recorded_parameter(tmp_path, fake_path_backend):
    run_neb(DATA / "hcn.xyz", DATA / "hnc.xyz", software="xtb",
            stage="relax", output=tmp_path)
    run_neb(ts_guess=DATA / "guess.xyz", software="xtb", stage="neb",
            images=3, interpolation="linear", output=tmp_path)
    summary_path = tmp_path / "neb_summary.json"
    summary = json.loads(summary_path.read_text())
    summary["parameters"]["future_initializer_control"] = 0.25
    summary_path.write_text(json.dumps(summary))

    with pytest.raises(ValueError, match="Unrecognized NEB stage parameter.*future_initializer_control"):
        run_neb(software="xtb", stage="ts", images=3, interpolation="linear",
                output=tmp_path)


@pytest.mark.skipif(
    find_spec("geodesic_interpolate") is None,
    reason="geodesic-interpolate is an optional dependency",
)
def test_geodesic_restart_tracks_active_options_and_package_version(tmp_path, fake_path_backend):
    from importlib.metadata import version

    run_neb(DATA / "hcn.xyz", DATA / "hnc.xyz", software="xtb",
            stage="relax", output=tmp_path)
    run_neb(ts_guess=DATA / "guess.xyz", software="xtb", stage="neb", images=3,
            interpolation="geodesic", idpp_fmax=0.1, idpp_steps=20,
            geodesic_tol=0.002, geodesic_max_iter=4, output=tmp_path)
    summary_path = tmp_path / "neb_summary.json"
    summary = json.loads(summary_path.read_text())
    assert summary["parameters"] == canonical_neb_parameters(summary["parameters"])
    assert summary["geodesic_interpolate_version"] == version("geodesic-interpolate")

    same = run_neb(software="xtb", stage="ts", images=3, interpolation="geodesic",
                   idpp_fmax=0.1, idpp_steps=20, geodesic_tol=0.002,
                   geodesic_max_iter=4, output=tmp_path)
    assert same["ts_optimization_converged"]
    with pytest.raises(ValueError, match="different stage-specific parameters"):
        run_neb(software="xtb", stage="ts", images=3, interpolation="geodesic",
                idpp_fmax=0.1, idpp_steps=20, geodesic_tol=0.003,
                geodesic_max_iter=4, output=tmp_path)

    summary = json.loads(summary_path.read_text())
    summary["geodesic_interpolate_version"] = "0.0.0"
    summary_path.write_text(json.dumps(summary))
    reused = run_neb(software="xtb", stage="ts", images=3, interpolation="geodesic",
                     idpp_fmax=0.1, idpp_steps=20, geodesic_tol=0.002,
                     geodesic_max_iter=4, output=tmp_path)
    assert reused["ts_optimization_converged"]


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
    summary_path = tmp_path / "ts_summary.json"
    historical_summary = json.loads(summary_path.read_text())
    historical_summary["qc_params"].pop("xtb_model")
    historical_summary["qc_params"].update({"method": "BP86", "basis": "def2-SVP"})
    summary_path.write_text(json.dumps(historical_summary))
    assert _load_stage(tmp_path, "ts", calculator)["ts_optimization_converged"]
    parallel = SimpleNamespace(qc_params=dict(calculator.qc_params, nprocs=8))
    assert _load_stage(tmp_path, "ts", parallel)["ts_optimization_converged"]
    charged = SimpleNamespace(qc_params=dict(calculator.qc_params, charge=1))
    with pytest.raises(ValueError, match="different backend"):
        _load_stage(tmp_path, "ts", charged)
    geometry.write_text((DATA / "hnc.xyz").read_text())
    with pytest.raises(ValueError, match="Stale ts artifact"):
        _load_stage(tmp_path, "ts", calculator)


def test_xtb_model_is_physical_provenance_and_controls_stage_reuse(
    tmp_path, fake_path_backend,
):
    run_neb(
        ts_geometry=DATA / "guess.xyz", software="xtb", stage="ts",
        xtb_model="gxtb", output=tmp_path,
    )
    summary_path = tmp_path / "ts_summary.json"
    summary = json.loads(summary_path.read_text())
    assert summary["qc_params"]["xtb_model"] == "gxtb"
    assert "method" not in summary["qc_params"]
    assert "basis" not in summary["qc_params"]
    assert summary["backend_model"] == "g-xTB (--gxtb)"

    # The unchanged model permits a downstream stage to consume the TS stage.
    frequency = run_neb(
        software="xtb", stage="frequency", xtb_model="gxtb", output=tmp_path,
    )
    assert frequency["first_order_saddle_confirmed"]

    # A different Hamiltonian must fail before any downstream calculation.
    with pytest.raises(ValueError, match="different backend or electronic-structure settings"):
        run_neb(
            software="xtb", stage="frequency", xtb_model="gfn2", output=tmp_path,
        )

    gfn2_output = tmp_path / "gfn2"
    run_neb(
        ts_geometry=DATA / "guess.xyz", software="xtb", stage="ts",
        xtb_model="gfn2", output=gfn2_output,
    )
    gfn2_summary = json.loads((gfn2_output / "ts_summary.json").read_text())
    assert gfn2_summary["qc_params"]["xtb_model"] == "gfn2"
    assert "method" not in gfn2_summary["qc_params"]
    assert "basis" not in gfn2_summary["qc_params"]
    assert gfn2_summary["backend_model"] == "GFN2-xTB (--gfn 2)"


def test_xtb_model_is_ignored_for_non_xtb_provenance():
    calculator = SimpleNamespace(qc_params={
        "software": "gaussian", "method": "B97-D", "xtb_model": "gfn2", "nprocs": 2,
    })
    assert _physical_settings(calculator) == {"software": "gaussian", "method": "B97-D"}


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
        calculator.backend_energy_gradient_evaluations = (
            getattr(calculator, "backend_energy_gradient_evaluations", 0) + 2
        )
        if "frequency_exception" in failures:
            raise RuntimeError("fixture frequency failure after backend calls")
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
        summary["parameters"].pop(key, None)
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


def test_neb_candidate_selection_prefers_valid_climber_over_global_maximum():
    frames = [np.full((2, 3), float(index)) for index in range(5)]
    energies = [0.0, 0.1, 0.2, 2.0, 0.0]
    selection = _select_neb_ts_candidate(
        frames, energies, climbing_activated=True, climber_indices=[2],
    )
    assert selection["highest_energy_image_index"] == 3
    assert selection["interior_maximum"] is True
    assert selection["climbing_image_indices"] == [2]
    assert selection["selected_climbing_image_index"] == 2
    assert selection["ts_guess_image_index"] == 2
    assert selection["ts_guess_source"] == "climbing_image"


@pytest.mark.parametrize("bad_index", [0, 4, -1, 1.5, True, "2"])
def test_neb_candidate_selection_rejects_invalid_climber_indices(bad_index):
    frames = [np.full((2, 3), float(index)) for index in range(5)]
    selection = _select_neb_ts_candidate(
        frames, [0.0, 0.1, 0.2, 1.0, 0.0],
        climbing_activated=True, climber_indices=[bad_index],
    )
    assert selection["climbing_image_indices"] == []
    assert selection["selected_climbing_image_index"] is None
    assert selection["ts_guess_image_index"] == 3
    assert selection["ts_guess_source"] == "highest_energy_interior_image"


def test_neb_candidate_selection_rejects_nonfinite_climber_geometry_and_energy():
    frames = [np.zeros((2, 3)) for _ in range(5)]
    energies = [0.0, 0.1, 0.2, 1.0, 0.0]
    frames[2][0, 0] = np.nan
    selection = _select_neb_ts_candidate(
        frames, energies, climbing_activated=True, climber_indices=[2],
    )
    assert selection["climbing_image_indices"] == []
    assert selection["ts_guess_source"] == "highest_energy_interior_image"
    assert selection["ts_guess_image_index"] == 3

    frames[2][0, 0] = 0.0
    energies[2] = np.inf
    selection = _select_neb_ts_candidate(
        frames, energies, climbing_activated=True, climber_indices=[2],
    )
    assert selection["climbing_image_indices"] == []
    assert selection["highest_energy_image_index"] == 2
    assert selection["ts_guess_image_index"] is None


def _configure_fake_neb_chain(monkeypatch, *, activated, climbers=(), energies=None):
    import geometric.neb

    original_band = geometric.neb.ElasticBand
    original_optimize = geometric.neb.OptimizeChain
    captured = {}

    def make_band(*args, **kwargs):
        band = original_band(*args, **kwargs)
        if activated is not None:
            band.climbSet = activated
        if climbers is not None:
            band.climbers = list(climbers)
        if energies is not None:
            for structure, energy in zip(band.Structures, energies):
                structure.energy = float(energy)
        return band

    def optimize_chain(chain, engine, params):
        captured["ncimg"] = params.ncimg
        return original_optimize(chain, engine, params)

    monkeypatch.setattr(geometric.neb, "ElasticBand", make_band)
    monkeypatch.setattr(geometric.neb, "OptimizeChain", optimize_chain)
    return captured


@pytest.mark.parametrize("optimizer", ["geometric", "sella"])
def test_active_climber_is_exact_neb_artifact_and_ts_optimizer_input(
    tmp_path, fake_path_backend, monkeypatch, optimizer,
):
    import pyar.neb as neb

    _, _ = fake_path_backend
    # The global maximum is image 3; geomeTRIC reports image 2 as climbing.
    energies = [0.0] * 11
    energies[3] = 2.0
    energies[5] = 1.0
    captured = _configure_fake_neb_chain(
        monkeypatch, activated=True, climbers=[2], energies=energies,
    )
    ts_input = {}
    if optimizer == "sella":
        def fake_sella(symbols, coordinates, calculator, max_steps, fmax):
            ts_input["coordinates"] = np.asarray(coordinates).copy()
            return [np.asarray(coordinates)], [-2.0], True, "test-sella"

        monkeypatch.setattr(neb, "_optimize_sella", fake_sella)
    else:
        original_optimize = neb._optimize_geometry

        def record_ts_input(symbols, coordinates, calculator, output, label, max_cycles, **kwargs):
            if label == "ts":
                ts_input["coordinates"] = np.asarray(coordinates).copy()
            return original_optimize(symbols, coordinates, calculator, output, label, max_cycles, **kwargs)

        monkeypatch.setattr(neb, "_optimize_geometry", record_ts_input)

    result = run_neb(
        DATA / "hcn.xyz", DATA / "hnc.xyz", DATA / "guess.xyz", software="xtb",
        output=tmp_path, ts_optimizer=optimizer,
    )
    summary = json.loads((tmp_path / "neb_summary.json").read_text())
    symbols, path_frames = _read_xyz_trajectory(tmp_path / "neb_path.xyz")
    guess_symbols, guess_coordinates, comment = read_xyz(tmp_path / "ts_guess.xyz")
    assert captured["ncimg"] == 1
    assert summary["climbing_image_state_available"] is True
    assert summary["climbing_image_activated"] is True
    assert summary["climbing_image_indices"] == [2]
    assert summary["selected_climbing_image_index"] == 2
    assert summary["highest_energy_image_index"] == 3
    assert summary["ts_guess_image_index"] == 2
    assert summary["ts_guess_source"] == "climbing_image"
    assert summary["requested_climbing_images"] == 1
    assert summary["interior_maximum"] is True
    assert "first_order_saddle_confirmed" not in summary
    assert result["frequency"]["first_order_saddle_confirmed"] is True
    assert symbols == guess_symbols
    np.testing.assert_array_equal(guess_coordinates, path_frames[2])
    path_comment = Path(tmp_path / "neb_path.xyz").read_text().splitlines()[2 * (len(symbols) + 2) + 1]
    assert float(comment.split("energy_hartree=")[1]) == float(
        path_comment.split("energy_hartree=")[1]
    )
    np.testing.assert_array_equal(ts_input["coordinates"], guess_coordinates)


def test_inactive_climbing_state_uses_highest_energy_interior_fallback(
    tmp_path, fake_path_backend, monkeypatch,
):
    _configure_fake_neb_chain(monkeypatch, activated=False)
    run_neb(DATA / "hcn.xyz", DATA / "hnc.xyz", DATA / "guess.xyz",
            software="xtb", stage="relax", output=tmp_path)
    run_neb(ts_guess=DATA / "guess.xyz", software="xtb", stage="neb", output=tmp_path)
    summary = json.loads((tmp_path / "neb_summary.json").read_text())
    assert summary["climbing_image_state_available"] is True
    assert summary["climbing_image_activated"] is False
    assert summary["climbing_image_indices"] == []
    assert summary["selected_climbing_image_index"] is None
    assert summary["highest_energy_image_index"] == 5
    assert summary["ts_guess_image_index"] == 5
    assert summary["ts_guess_source"] == "highest_energy_interior_image"


@pytest.mark.parametrize("climbers", [[], [0], [10], [2.5], ["bad"]])
def test_active_but_unusable_climber_deliberately_falls_back(
    tmp_path, fake_path_backend, monkeypatch, climbers,
):
    _configure_fake_neb_chain(monkeypatch, activated=True, climbers=climbers)
    run_neb(DATA / "hcn.xyz", DATA / "hnc.xyz", software="xtb",
            stage="relax", output=tmp_path)
    run_neb(ts_guess=DATA / "guess.xyz", software="xtb", stage="neb", output=tmp_path)
    summary = json.loads((tmp_path / "neb_summary.json").read_text())
    assert summary["climbing_image_activated"] is True
    assert summary["climbing_image_indices"] == []
    assert summary["selected_climbing_image_index"] is None
    assert summary["ts_guess_image_index"] == summary["highest_energy_image_index"] == 5
    assert summary["ts_guess_source"] == "highest_energy_interior_image"


def test_multiple_valid_climbers_choose_highest_energy_climber_deterministically(
    tmp_path, fake_path_backend, monkeypatch,
):
    energies = [0.0] * 11
    energies[2], energies[3], energies[5] = 0.5, 0.75, 1.0
    _configure_fake_neb_chain(monkeypatch, activated=True, climbers=[3, 2], energies=energies)
    run_neb(DATA / "hcn.xyz", DATA / "hnc.xyz", software="xtb",
            stage="relax", output=tmp_path)
    run_neb(ts_guess=DATA / "guess.xyz", software="xtb", stage="neb", output=tmp_path)
    summary = json.loads((tmp_path / "neb_summary.json").read_text())
    assert summary["climbing_image_indices"] == [3, 2]
    assert summary["selected_climbing_image_index"] == 3
    assert summary["highest_energy_image_index"] == 5
    assert summary["ts_guess_image_index"] == 3


def test_missing_geometric_climbing_state_is_recorded_as_unavailable_and_falls_back(
    tmp_path, fake_path_backend, monkeypatch,
):
    _configure_fake_neb_chain(monkeypatch, activated=None, climbers=None)
    run_neb(DATA / "hcn.xyz", DATA / "hnc.xyz", software="xtb",
            stage="relax", output=tmp_path)
    run_neb(ts_guess=DATA / "guess.xyz", software="xtb", stage="neb", output=tmp_path)
    summary = json.loads((tmp_path / "neb_summary.json").read_text())
    assert summary["climbing_image_state_available"] is False
    assert summary["climbing_image_activated"] is None
    assert summary["climbing_image_indices"] == []
    assert summary["ts_guess_source"] == "highest_energy_interior_image"


def test_old_schema_two_neb_summary_without_climbing_fields_remains_reusable(
    tmp_path, fake_path_backend,
):
    run_neb(DATA / "hcn.xyz", DATA / "hnc.xyz", software="xtb",
            stage="relax", output=tmp_path)
    run_neb(ts_guess=DATA / "guess.xyz", software="xtb", stage="neb", output=tmp_path)
    summary_path = tmp_path / "neb_summary.json"
    summary = json.loads(summary_path.read_text())
    for key in (
        "climbing_image_state_available", "climbing_image_activated",
        "climbing_image_indices", "selected_climbing_image_index",
        "requested_climbing_images", "ts_guess_image_index", "ts_guess_source",
    ):
        summary.pop(key)
    summary_path.write_text(json.dumps(summary))
    assert run_neb(software="xtb", stage="ts", output=tmp_path)["ts_optimization_converged"]
    unchanged = json.loads(summary_path.read_text())
    assert all(key not in unchanged for key in (
        "climbing_image_state_available", "climbing_image_activated", "ts_guess_source",
    ))


def test_endpoint_energy_maximum_does_not_create_ts_guess(
    tmp_path, fake_path_backend, monkeypatch,
):
    energies = [2.0] + [0.0] * 10
    _configure_fake_neb_chain(monkeypatch, activated=False, energies=energies)
    run_neb(DATA / "hcn.xyz", DATA / "hnc.xyz", software="xtb",
            stage="relax", output=tmp_path)
    run_neb(ts_guess=DATA / "guess.xyz", software="xtb", stage="neb", output=tmp_path)
    summary = json.loads((tmp_path / "neb_summary.json").read_text())
    assert summary["highest_energy_image_index"] == 0
    assert summary["interior_maximum"] is False
    assert summary["ts_guess_image_index"] is None
    assert summary["ts_guess_source"] is None
    assert not (tmp_path / "ts_guess.xyz").exists()
    with pytest.raises(ValueError, match="converged NEB with an interior maximum"):
        run_neb(software="xtb", stage="ts", output=tmp_path)


def test_converged_neb_uses_valid_climber_even_when_global_maximum_is_endpoint(
    tmp_path, fake_path_backend, monkeypatch,
):
    energies = [2.0] + [0.0] * 10
    _configure_fake_neb_chain(
        monkeypatch, activated=True, climbers=[2], energies=energies,
    )
    result = run_neb(
        DATA / "hcn.xyz", DATA / "hnc.xyz", DATA / "guess.xyz",
        software="xtb", output=tmp_path,
    )
    summary = json.loads((tmp_path / "neb_summary.json").read_text())
    assert summary["converged"] is True
    assert summary["interior_maximum"] is False
    assert summary["highest_energy_image_index"] == 0
    assert summary["ts_guess_image_index"] == 2
    assert summary["ts_guess_source"] == "climbing_image"
    assert result["ts"]["ts_optimization_converged"] is True


def test_changing_climb_threshold_invalidates_neb_restart(tmp_path, fake_path_backend):
    run_neb(DATA / "hcn.xyz", DATA / "hnc.xyz", software="xtb",
            stage="relax", output=tmp_path)
    run_neb(ts_guess=DATA / "guess.xyz", software="xtb", stage="neb",
            images=3, climb=0.5, output=tmp_path)
    with pytest.raises(ValueError, match="Stage neb used different stage-specific parameters"):
        run_neb(software="xtb", stage="ts", images=3, climb=0.4, output=tmp_path)


@pytest.mark.parametrize("interpolation", ["linear", "idpp", "geodesic"])
def test_climbing_handoff_is_independent_of_interpolation(
    tmp_path, fake_path_backend, monkeypatch, interpolation,
):
    if interpolation == "geodesic" and find_spec("geodesic_interpolate") is None:
        pytest.skip("geodesic-interpolate is an optional dependency")
    kwargs = {"geodesic_max_iter": 3} if interpolation == "geodesic" else {}
    _configure_fake_neb_chain(monkeypatch, activated=True, climbers=[2])
    run_neb(DATA / "hcn.xyz", DATA / "hnc.xyz", software="xtb",
            stage="relax", output=tmp_path)
    run_neb(ts_guess=DATA / "guess.xyz", software="xtb", stage="neb",
            interpolation=interpolation, **kwargs, output=tmp_path)
    summary = json.loads((tmp_path / "neb_summary.json").read_text())
    assert summary["ts_guess_image_index"] == 2
    assert summary["ts_guess_source"] == "climbing_image"


@pytest.mark.parametrize("xtb_model", ["gxtb", "gfn2"])
def test_climbing_handoff_is_independent_of_xtb_model(
    tmp_path, fake_path_backend, monkeypatch, xtb_model,
):
    _configure_fake_neb_chain(monkeypatch, activated=True, climbers=[2])
    run_neb(DATA / "hcn.xyz", DATA / "hnc.xyz", software="xtb", xtb_model=xtb_model,
            stage="relax", output=tmp_path)
    run_neb(ts_guess=DATA / "guess.xyz", software="xtb", xtb_model=xtb_model,
            stage="neb", output=tmp_path)
    summary = json.loads((tmp_path / "neb_summary.json").read_text())
    assert summary["qc_params"]["xtb_model"] == xtb_model
    assert summary["ts_guess_image_index"] == 2


@pytest.mark.skipif(
    find_spec("geodesic_interpolate") is None,
    reason="geodesic-interpolate is an optional dependency",
)
def test_geodesic_workflow_keeps_geometric_neb_and_ts_handoff(
    tmp_path, fake_path_backend, monkeypatch,
):
    import pyar.neb as neb

    calls, _ = fake_path_backend
    original_optimize = neb._optimize_geometry
    ts_input = {}

    def record_ts_input(symbols, coordinates, calculator, output, label, max_cycles, **kwargs):
        if label == "ts":
            ts_input["coordinates"] = np.asarray(coordinates).copy()
        return original_optimize(symbols, coordinates, calculator, output, label, max_cycles, **kwargs)

    monkeypatch.setattr(neb, "_optimize_geometry", record_ts_input)
    result = run_neb(
        DATA / "hcn.xyz", DATA / "hnc.xyz", DATA / "guess.xyz", software="xtb",
        output=tmp_path, interpolation="geodesic", geodesic_max_iter=3,
    )

    assert "neb" in calls and "ts" in calls
    assert result["neb"]["geodesic_interpolate_version"]
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


def test_failed_stage_summary_retains_provider_calls_and_elapsed_time(tmp_path, fake_path_backend):
    _, failures = fake_path_backend
    failures.add("frequency_exception")
    with pytest.raises(RuntimeError, match="after backend calls"):
        run_neb(
            ts_geometry=DATA / "guess.xyz", software="xtb", stage="frequency",
            output=tmp_path,
        )
    summary = json.loads((tmp_path / "frequency_summary.json").read_text())
    assert summary["status"] == "failed"
    assert summary["backend_energy_gradient_evaluations"] == 2
    assert summary["wall_seconds"] >= 0


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
    assert result["hessian_source"] == "finite_difference_cartesian"
    assert result["hessian_evaluations"] is None
    assert result["hessian_wall_seconds"] >= 0
    assert not result["stationary"]
    assert not result["minimum_confirmed"]


def test_numerically_linear_geometry_retains_both_bending_modes(tmp_path):
    symbols, geometry, _ = read_xyz(DATA / "hcn.xyz")
    distorted = geometry.copy()
    distorted[0, 1] += 1e-6
    result = _frequency(symbols, distorted, HarmonicCalculator(geometry), tmp_path, "test", 20.)
    assert 0 < result["linear_geometry_correction_angstrom"] < 1e-5
    assert len(result["frequencies_cm-1"]) == 4
