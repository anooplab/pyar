"""Backend construction contracts; stubs do not validate native chemistry."""

import builtins
from copy import deepcopy
from types import ModuleType, SimpleNamespace
from unittest.mock import Mock
import sys

import pytest

from pyar import backends, optimiser
from pyar.backend_errors import BackendExecutionError


# Independent public routes, including historical aliases and the OBabel signature.
ROUTES = [
    ("mlatom_aiqm1", "aiqm1_mlatom", "AIQM1"),
    ("gaussian", "gaussian", "Gaussian"),
    ("g16", "gaussian", "Gaussian"),
    ("orca", "orca", "Orca"),
    ("orca16", "orca", "Orca"),
    ("orca-aiqm1", "orca_aiqm1", "OrcaAIQM1"),
    ("psi4", "psi4", "Psi4"),
    ("xtb", "xtb", "Xtb"),
    ("xtb_turbo", "xtb_turbo", "XtbTurbo"),
    ("turbomole", "turbomole", "Turbomole"),
    ("mopac", "mopac", "Mopac"),
    ("aimnet_2", "aimnet_2", "Aimnet2"),
    ("aiqm1_mlatom", "aiqm1_mlatom", "AIQM1"),
    ("xtb-aimnet2", "xtb_aimnet2", "XtbAimnet2"),
    ("xtb_aimnet2", "xtb_aimnet2", "XtbAimnet2"),
    ("xtb-aiqm1", "xtb_aiqm1", "XtbAIQM1"),
    ("xtb_aiqm1", "xtb_aiqm1", "XtbAIQM1"),
    ("obabel", "babel", "OBabel"),
]


@pytest.fixture
def constructors(monkeypatch):
    result = {}
    for _, module_name, class_name in ROUTES + [
        ("", "geometric", "Geometric"),
        ("", "mlatom_aiqm1", "MlatomAiqm1"),
    ]:
        if module_name in result:
            continue
        module = ModuleType(f"pyar.backends.{module_name}")
        instance = SimpleNamespace(set_gamma=Mock(), backend=module_name)
        constructor = Mock(return_value=instance)
        setattr(module, class_name, constructor)
        monkeypatch.setitem(sys.modules, module.__name__, module)
        monkeypatch.setattr(backends, module_name, module, raising=False)
        result[module_name] = constructor
    return result


def assert_only_constructor(constructors, selected):
    for name, constructor in constructors.items():
        assert constructor.call_count == (1 if name == selected else 0), name


@pytest.mark.parametrize("software,module_name,class_name", ROUTES)
def test_native_route_preserves_input_identity_and_settings(
    constructors, software, module_name, class_name
):
    molecule = SimpleNamespace(charge=1, multiplicity=1, name="charged-input")
    settings = {"software": software, "method": "owner-method", "basis": "owner-basis"}
    before = deepcopy(settings)

    result = optimiser.build_geometry(molecule, settings)

    assert result is constructors[module_name].return_value
    expected_args = (molecule,) if module_name == "babel" else (molecule, settings)
    actual_args = constructors[module_name].call_args.args
    assert all(
        actual is expected for actual, expected in zip(actual_args, expected_args)
    )
    assert len(actual_args) == len(expected_args)
    assert constructors[module_name].call_args.kwargs == {}
    assert_only_constructor(constructors, module_name)
    assert settings == before
    assert (molecule.charge, molecule.multiplicity) == (1, 1)


@pytest.mark.parametrize(
    "gamma,selected", [(None, "xtb_turbo"), (0.0, "xtb"), (1.0, "xtb_turbo")]
)
def test_xtb_turbo_zero_gamma_is_the_native_xtb_route(constructors, gamma, selected):
    molecule = object()
    settings = {"software": "xtb_turbo", "gamma": gamma}
    assert (
        optimiser.build_geometry(molecule, settings)
        is constructors[selected].return_value
    )
    constructors[selected].assert_called_once_with(molecule, settings)
    assert_only_constructor(constructors, selected)


@pytest.mark.parametrize("gamma", [None, 0.0, 100.0])
def test_orca_aiqm1_gamma_is_forwarded_only_when_supplied(constructors, gamma):
    settings = {"software": "orca-aiqm1", "gamma": gamma}
    result = optimiser.build_geometry(object(), settings)
    if gamma is None:
        result.set_gamma.assert_not_called()
    else:
        result.set_gamma.assert_called_once_with(gamma)
    assert_only_constructor(constructors, "orca_aiqm1")


@pytest.mark.parametrize("software", ["xtb", "orca", "orca16", "gaussian", "g16"])
def test_geometric_precedes_native_dispatch_with_unchanged_request(
    constructors, software
):
    molecule = object()
    settings = {"software": software, "geometry_optimizer": "geometric", "gamma": 100.0}
    before = deepcopy(settings)
    assert (
        optimiser.build_geometry(molecule, settings)
        is constructors["geometric"].return_value
    )
    constructors["geometric"].assert_called_once_with(molecule, settings)
    assert_only_constructor(constructors, "geometric")
    assert settings == before


@pytest.mark.parametrize(
    "software,optimizer",
    [("unsupported", "native"), ("psi4", "geometric"), ("unsupported", "geometric")],
)
def test_unsupported_route_never_constructs_or_falls_back(
    constructors, software, optimizer
):
    with pytest.raises(ValueError):
        optimiser.build_geometry(
            object(), {"software": software, "geometry_optimizer": optimizer}
        )
    assert_only_constructor(constructors, None)


@pytest.mark.parametrize(
    "error",
    [
        BackendExecutionError("backend failed"),
        ValueError("bad settings"),
        ModuleNotFoundError("missing extra"),
    ],
)
def test_constructor_failure_is_propagated_without_fallback(constructors, error):
    constructors["xtb"].side_effect = error
    with pytest.raises(type(error)) as caught:
        optimiser.build_geometry(object(), {"software": "xtb"})
    assert caught.value is error
    assert_only_constructor(constructors, "xtb")


def test_optional_import_failure_is_not_replaced_by_another_backend(
    monkeypatch, constructors
):
    original_import = builtins.__import__
    error = ModuleNotFoundError("No module named 'torch'", name="torch")

    def import_without_aimnet(name, globals=None, locals=None, fromlist=(), level=0):
        if name == "pyar.backends" and "aimnet_2" in fromlist:
            raise error
        return original_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr(builtins, "__import__", import_without_aimnet)
    with pytest.raises(ModuleNotFoundError) as caught:
        optimiser.build_geometry(object(), {"software": "aimnet_2"})
    assert caught.value is error
    assert_only_constructor(constructors, None)


def test_missing_backend_setting_does_not_guess_a_default(constructors):
    with pytest.raises(KeyError, match="software"):
        optimiser.build_geometry(object(), {})
    assert_only_constructor(constructors, None)
