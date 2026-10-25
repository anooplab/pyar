"""Manual-backed adapter boundaries and opt-in installed-backend calculations."""

import json
import os
from pathlib import Path

import numpy as np
import pytest

from pyar.core.molecule import Molecule
from pyar.scripts.singlepoint import build_parser
from pyar.singlepoint import native, request, service, worker
from pyar.singlepoint.io import BOHR_ANGSTROM, write_json


@pytest.fixture
def molecule():
    return Molecule(["H", "H"], [[0, 0, 0], [0.74, 0, 0]], name="h2")


def settings_for(tmp_path, backend):
    path = tmp_path / "h2.xyz"
    path.write_text("2\nH2\nH 0 0 0\nH 0.74 0 0\n")
    args = build_parser().parse_args(
        [str(path), "--backend", backend, "--output", str(tmp_path / "out")]
    )
    return request.resolve(args, "energy")[0]


@pytest.mark.parametrize("gradient", [False, True])
def test_xtb_adapter(tmp_path, molecule, monkeypatch, gradient):
    settings = settings_for(tmp_path, "xtb")

    def run(command, directory, settings):
        assert "--gfn" in command and "--uhf" in command
        assert ("--grad" in command) == gradient
        (directory / "gradient").write_text("0.1 0.2 0.3\n0.4 0.5 0.6\n")
        return "xtb version 6.7.1\nHamiltonian GFN2-xTB\n| TOTAL ENERGY -1.23 Eh |\nnormal termination of xtb"

    monkeypatch.setattr(native, "run_program", run)
    monkeypatch.setattr(native, "require_executable", lambda _: "/mock/xtb")
    energy, g, meta = native.xtb(
        settings, molecule, molecule.coordinates / BOHR_ANGSTROM, tmp_path, gradient
    )
    assert energy == -1.23
    assert meta["version"] == "6.7.1"
    if gradient:
        assert np.array_equal(g, [[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]])


@pytest.mark.parametrize("gradient", [False, True])
def test_orca_adapter(tmp_path, molecule, monkeypatch, gradient):
    settings = settings_for(tmp_path, "orca")

    def run(command, directory, settings):
        (directory / "input.engrad").write_text(
            "2\n-1.23\n0.1\n0.2\n0.3\n0.4\n0.5\n0.6\n"
        )
        return "Program Version 6.1.1\nFINAL SINGLE POINT ENERGY -1.23\nORCA TERMINATED NORMALLY"

    monkeypatch.setattr(native, "run_program", run)
    monkeypatch.setattr(native, "require_executable", lambda _: "/mock/orca")
    energy, g, meta = native.orca(
        settings, molecule, molecule.coordinates / BOHR_ANGSTROM, tmp_path, gradient
    )
    assert energy == -1.23 and meta["version"] == "6.1.1"
    assert (g is not None) == gradient


def test_orca_inconsistent_energy(tmp_path, molecule, monkeypatch):
    settings = settings_for(tmp_path, "orca")
    (tmp_path / "input.engrad").write_text("2\n-2\n0\n0\n0\n0\n0\n0\n")
    monkeypatch.setattr(native, "require_executable", lambda _: "/mock/orca")
    monkeypatch.setattr(
        native,
        "run_program",
        lambda *a: "FINAL SINGLE POINT ENERGY -1\nORCA TERMINATED NORMALLY",
    )
    with pytest.raises(ValueError, match="disagree"):
        native.orca(settings, molecule, molecule.coordinates, tmp_path, True)


@pytest.mark.parametrize("gradient", [False, True])
def test_gaussian_adapter(tmp_path, molecule, monkeypatch, gradient):
    settings = settings_for(tmp_path, "gaussian")

    def run(command, directory, settings, **kwargs):
        assert "PBE1PBE/def2SVP" in kwargs["stdin"]
        assert (" Force" in kwargs["stdin"]) == gradient
        return (
            "Gaussian 16 Revision C.01\n SCF Done: E(RHF) = -1.23 A.U.\n"
            "Forces (Hartrees/Bohr)\n1 1 0.1 0.2 0.3\n2 1 -0.1 -0.2 -0.3\n"
            "Normal termination of Gaussian"
        )

    monkeypatch.setattr(native, "run_program", run)
    monkeypatch.setattr(native, "require_executable", lambda _: "/mock/g16")
    energy, g, _ = native.gaussian(
        settings, molecule, molecule.coordinates, tmp_path, gradient
    )
    assert energy == -1.23
    if gradient:
        assert np.array_equal(g, [[-0.1, -0.2, -0.3], [0.1, 0.2, 0.3]])


def test_correlated_gaussian_rejected(tmp_path):
    path = tmp_path / "h2.xyz"
    path.write_text("2\nh2\nH 0 0 0\nH .74 0 0\n")
    args = build_parser().parse_args(
        [str(path), "--backend", "gaussian", "--method", "MP2"]
    )
    with pytest.raises(ValueError, match="HF and these SCF DFT"):
        request.resolve(args, "energy")


def test_python_worker_boundary(tmp_path, molecule, monkeypatch):
    settings = settings_for(tmp_path, "ani")

    def run(command, directory, settings):
        payload = json.loads((directory / "worker_request.json").read_text())
        assert payload["symbols"] == ["H", "H"]
        assert np.allclose(payload["coordinates_angstrom"], molecule.coordinates)
        write_json(
            directory / "worker_result.json",
            {
                "energy_hartree": -1.0,
                "gradient_hartree_per_bohr": [[0.0, 0.0, 0.0]] * 2,
                "metadata": {"version": "test"},
            },
        )

    monkeypatch.setattr(native, "run_program", run)
    energy, gradient, meta = native.python_backend(
        settings, molecule, molecule.coordinates / BOHR_ANGSTROM, tmp_path, True
    )
    assert energy == -1 and meta["version"] == "test"


def test_gxtb_preflight_rejects_missing_capability(monkeypatch):
    monkeypatch.setattr(request.shutil, "which", lambda _: "/mock/xtb")
    monkeypatch.setattr("pyar.backends.xtb.xtb_supports_gxtb", lambda _: False)
    with pytest.raises(ValueError, match="does not advertise"):
        request.preflight({"software": "xtb", "xtb_model": "gxtb"})


def test_model_state_constraints(molecule):
    molecule.multiplicity = 3
    with pytest.raises(ValueError, match="singlets"):
        request.validate_model_state("aimnet_2", molecule)
    molecule.multiplicity = 1
    molecule.charge = 2
    with pytest.raises(ValueError, match="neutral"):
        request.validate_model_state("ani", molecule)


@pytest.mark.external_backend
@pytest.mark.tier2
@pytest.mark.skipif(
    os.environ.get("PYAR_TEST_REAL_SINGLEPOINT") != "1",
    reason="opt-in real backend calculations",
)
@pytest.mark.parametrize("backend", ["xtb", "orca", "aimnet_2", "ani"])
def test_real_energy_gradient_and_minimum(tmp_path, molecule, backend):
    """Numerically locate an H2 minimum, then independently check its Hessian.

    This is a smoke/contract check, not a scientific accuracy benchmark.
    """
    from scipy.optimize import minimize_scalar

    settings = settings_for(tmp_path, backend)
    if backend == "orca":
        settings.update(method="HF", basis="STO-3G")
    request.preflight(settings)
    evaluate = service.make_evaluator(settings, molecule, tmp_path)
    if backend == "ani":
        molecule = optimized_ani_water()
    else:
        counter = iter(range(100))

        def objective(distance):
            coords = np.array([[0.0, 0.0, 0.0], [distance / BOHR_ANGSTROM, 0.0, 0.0]])
            return evaluate(coords, f"energy_{next(counter)}", gradient=False)[0]

        optimum = minimize_scalar(
            objective, bounds=(0.6, 1.0), method="bounded", options={"xatol": 2e-5}
        )
        molecule.coordinates[1, 0] = optimum.x
    # Rebind the evaluator to the actual molecular composition after test-only optimization.
    evaluate = service.make_evaluator(settings, molecule, tmp_path)
    coordinates = molecule.coordinates / BOHR_ANGSTROM
    energy, gradient = evaluate(coordinates)
    step = 0.01 if backend == "aimnet_2" else 0.002
    plus, minus = coordinates.copy(), coordinates.copy()
    plus[1, 0] += step
    minus[1, 0] -= step
    ep = evaluate(plus, "energy_plus", gradient=False)[0]
    em = evaluate(minus, "energy_minus", gradient=False)[0]
    assert gradient[1, 0] == pytest.approx((ep - em) / (2 * step), abs=3e-4)
    result = service.frequency_result(
        evaluate,
        molecule,
        coordinates,
        gradient,
        tmp_path,
        {
            "max_gradient": 4.5e-4,
            "rms_gradient": 3e-4,
            "step": step,
            "imaginary_threshold": 20.0,
        },
    )
    assert result["minimum_confirmed"], result
    expected_modes = 3 if backend == "ani" else 1
    assert len(result["frequencies_cm-1"]) == expected_modes
    assert min(result["frequencies_cm-1"]) > 0


def optimized_ani_water():
    from ase import Atoms
    from ase.optimize import BFGS

    atoms = Atoms(
        "OH2", positions=[[0.0, 0.0, 0.0], [0.957, 0.0, 0.0], [-0.24, 0.927, 0.0]]
    )
    atoms.calc = worker.ani_calculator({"settings": {"model": "ANI1x"}})
    optimizer = BFGS(atoms, logfile=None)
    assert optimizer.run(fmax=1e-4, steps=100)
    return Molecule(atoms.get_chemical_symbols(), atoms.get_positions(), name="water")


def test_worker_calculates_ase_units(tmp_path, monkeypatch):
    from ase.calculators.calculator import Calculator, all_changes
    from ase.units import Hartree, Bohr

    class Harmonic(Calculator):
        implemented_properties = ["energy", "forces"]

        def calculate(
            self, atoms=None, properties=("energy",), system_changes=all_changes
        ):
            super().calculate(atoms, properties, system_changes)
            self.results = {
                "energy": 2 * Hartree,
                "forces": np.ones((2, 3)) * Hartree / Bohr,
            }

    monkeypatch.chdir(tmp_path)
    monkeypatch.setitem(worker.CALCULATORS, "ani", lambda _: Harmonic())
    monkeypatch.setattr(worker, "calculator_metadata", lambda _: {"backend": "test"})
    write_json(
        "worker_request.json",
        {
            "settings": {"software": "ani"},
            "symbols": ["H", "H"],
            "coordinates_angstrom": [[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]],
            "gradient": True,
        },
    )
    worker.main()
    result = json.loads(Path("worker_result.json").read_text())
    assert result["energy_hartree"] == 2.0
    assert np.allclose(result["gradient_hartree_per_bohr"], -np.ones((2, 3)))


def test_model_factory_contracts(monkeypatch):
    from types import SimpleNamespace
    import sys

    fake = SimpleNamespace(double=lambda: SimpleNamespace(ase=lambda: "ani calculator"))
    monkeypatch.setitem(
        sys.modules,
        "torchani",
        SimpleNamespace(models=SimpleNamespace(ANI1x=lambda: fake)),
    )
    assert worker.ani_calculator({"settings": {"model": "ANI1x"}}) == "ani calculator"
    request.validate_ani_model("ANI1x")
    with pytest.raises(ValueError, match="does not provide"):
        request.validate_ani_model("ANI2x")


def test_metadata_preserves_model_hash(monkeypatch, tmp_path):
    path = tmp_path / "model.jpt"
    path.write_bytes(b"model bytes")
    monkeypatch.setattr("pyar.backends.aimnet2_assets.MODEL_PATH", str(path))
    monkeypatch.setattr(worker, "version", lambda _: "test-version")
    result = worker.calculator_metadata({"settings": {"software": "aimnet_2"}})
    assert result["model"] == "model.jpt"
    assert len(result["model_sha256"]) == 64
    assert result["torch_version"] == "test-version"
    result = worker.calculator_metadata(
        {"settings": {"software": "ani", "model": "ANI1x"}}
    )
    assert result["model"] == "ANI1x"
    assert result["torchani_version"] == "test-version"


def test_xtb_rejects_wrong_hamiltonian(tmp_path, molecule, monkeypatch):
    settings = settings_for(tmp_path, "xtb")
    monkeypatch.setattr(native, "require_executable", lambda _: "/mock/xtb")
    monkeypatch.setattr(
        native, "run_program", lambda *a: "GFN1-xTB\nnormal termination of xtb"
    )
    with pytest.raises(RuntimeError, match="Hamiltonian"):
        native.xtb(settings, molecule, molecule.coordinates, tmp_path, False)


def test_backend_capability_reporting(capsys):
    from pyar.scripts.cli_meta import backends

    backends(["ani", "--json"])
    entry = json.loads(capsys.readouterr().out)["backends"]["ani"]
    assert (
        entry["single_point_energy"]
        and entry["single_point_gradient"]
        and entry["numerical_frequency"]
    )
    assert not entry["energy_gradient"]  # Legacy optimizer capability is unchanged.


def test_aimnet_factory_without_loading_assets(monkeypatch):
    from types import SimpleNamespace
    import sys

    fake_module = SimpleNamespace(load_default_aimnet2_model=lambda: "model")
    monkeypatch.setitem(sys.modules, "pyar.backends.aimnet_2", fake_module)
    monkeypatch.setitem(
        sys.modules,
        "pyar.AIMNet2.calculators.aimnet2ase",
        SimpleNamespace(AIMNet2Calculator=lambda model, charge: (model, charge)),
    )
    assert worker.aimnet_calculator({"charge": 1}) == ("model", 1)


def test_extra_native_controls(tmp_path, molecule, monkeypatch):
    settings = settings_for(tmp_path, "orca")
    settings.update(dispersion="D4", basis=None)
    text = native.orca_input(settings, molecule, molecule.coordinates, False)
    assert " D4" in text
    settings = settings_for(tmp_path, "xtb")
    settings.update(solvent="water", solvent_model="alpb")

    def run(command, *args):
        assert command[-2:] == ["--alpb", "water"]
        return "GFN2-xTB\n| TOTAL ENERGY -1.0 Eh |\nnormal termination of xtb"

    monkeypatch.setattr(native, "run_program", run)
    monkeypatch.setattr(native, "require_executable", lambda _: "/mock/xtb")
    assert (
        native.xtb(settings, molecule, molecule.coordinates, tmp_path, False)[0] == -1.0
    )


def test_turbomole_missing_executable_is_actionable(monkeypatch):
    monkeypatch.setattr(
        request.shutil,
        "which",
        lambda name: None if name == "dscf" else "/mock/program",
    )
    with pytest.raises(ValueError, match="Missing executable 'dscf'"):
        request.preflight({"software": "turbomole"})


def test_preflight_assets_and_requirements(monkeypatch):
    monkeypatch.setattr(request, "module_available", lambda _: True)
    monkeypatch.setattr(request.shutil, "which", lambda _: "/mock/program")
    monkeypatch.setattr(
        "pyar.backends.aimnet2_assets.validate_aimnet2_runtime_assets", lambda: None
    )
    monkeypatch.setattr(request, "validate_ani_model", lambda _: None)
    monkeypatch.setattr("pyar.backends.xtb.xtb_supports_gxtb", lambda _: True)
    assert "torch" in request.preflight({"software": "aimnet_2"})
    assert "torchani" in request.preflight({"software": "ani", "model": "ANI1x"})
    assert "ase" in request.module_requirements("turbomole")
    assert "xtb" in request.preflight({"software": "xtb", "xtb_model": "gxtb"})


def test_worker_energy_only(tmp_path, monkeypatch):
    from ase.calculators.emt import EMT

    monkeypatch.chdir(tmp_path)
    monkeypatch.setitem(worker.CALCULATORS, "test", lambda _: EMT())
    monkeypatch.setattr(worker, "calculator_metadata", lambda _: {})
    write_json(
        "worker_request.json",
        {
            "settings": {"software": "test"},
            "symbols": ["H", "H"],
            "coordinates_angstrom": [[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]],
            "gradient": False,
        },
    )
    worker.main()
    assert (
        json.loads(Path("worker_result.json").read_text())["gradient_hartree_per_bohr"]
        is None
    )


@pytest.mark.external_backend
@pytest.mark.tier2
@pytest.mark.skipif(
    os.environ.get("PYAR_TEST_REAL_SINGLEPOINT") != "1",
    reason="opt-in real backend calculations",
)
@pytest.mark.parametrize(
    "backend,model",
    [("xtb", "alpb"), ("xtb", "gbsa"), ("orca", "CPCM"), ("orca", "SMD")],
)
def test_real_solvent_singlepoints(tmp_path, backend, model):
    molecule = Molecule(
        ["O", "H", "H"],
        [[0.0, 0.0, 0.0], [0.957, 0.0, 0.0], [-0.24, 0.927, 0.0]],
        name="water",
    )
    settings = settings_for(tmp_path, backend)
    request.preflight(settings)
    evaluate = service.make_evaluator(settings, molecule, tmp_path)
    coordinates = molecule.coordinates / BOHR_ANGSTROM
    gas_energy, gas_gradient = evaluate(coordinates, "gas")
    settings.update(solvent="water", solvent_model=model)
    solvent_energy, solvent_gradient = evaluate(coordinates, "solvent")
    assert abs(solvent_energy - gas_energy) > 1e-6
    assert gas_gradient.shape == solvent_gradient.shape == (3, 3)
    # Confirm solvent settings also apply to energy-only calls, not just gradients.
    energy_only, absent = evaluate(coordinates, "solvent_energy", gradient=False)
    assert energy_only == pytest.approx(solvent_energy, abs=2e-7)
    assert absent is None


@pytest.mark.external_backend
@pytest.mark.tier2
@pytest.mark.skipif(
    os.environ.get("PYAR_TEST_REAL_SINGLEPOINT") != "1",
    reason="opt-in real backend calculations",
)
@pytest.mark.parametrize("backend", ["xtb", "orca"])
def test_real_charged_doublet(tmp_path, molecule, backend):
    molecule.charge, molecule.multiplicity = 1, 2
    settings = settings_for(tmp_path, backend)
    request.preflight(settings)
    evaluate = service.make_evaluator(settings, molecule, tmp_path)
    energy, gradient = evaluate(molecule.coordinates / BOHR_ANGSTROM)
    assert np.isfinite(energy)
    assert gradient.shape == (2, 3)
