"""Direct read-only path preflight contracts; no calculators or QC are run."""

import builtins
from copy import deepcopy
from types import SimpleNamespace
from unittest.mock import Mock

import numpy as np
import pytest

from pyar.core.molecule import Molecule
from pyar import neb, optimization_request
from pyar.workflows import reaction_characterization as workflow


@pytest.fixture
def preflight_context(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    reactants = [
        Molecule(["He"], [[0.0, 0.0, 0.0]], charge=0, multiplicity=1),
        Molecule(["H"], [[3.0, 0.0, 0.0]], charge=0, multiplicity=2, scftype="uhf"),
    ]
    settings = dict(
        software="xtb",
        method=None,
        basis=None,
        nprocs=2,
        scf_cycles=100,
        opt_cycles=80,
        opt_threshold="tight",
        xtb_model="gfn2",
        gamma=100.0,
        bias_controller="adaptive",
        bias_potential="afir",
    )
    request = SimpleNamespace(reactants=reactants, qc_params=settings)
    backend = Mock(return_value=["xtb", "geometric", "xtb"])
    monkeypatch.setattr(optimization_request, "preflight", backend)
    execution = Mock(side_effect=AssertionError("preflight must not execute stages"))
    monkeypatch.setattr(workflow, "run_path_stages", execution)
    return SimpleNamespace(
        request=request, backend=backend, execution=execution, root=tmp_path
    )


def assert_unchanged(context, settings, coordinates):
    assert context.request.qc_params == settings
    for molecule, expected in zip(context.request.reactants, coordinates):
        np.testing.assert_array_equal(molecule.coordinates, expected)
    assert list(context.root.iterdir()) == []
    context.execution.assert_not_called()


@pytest.mark.parametrize(
    "through",
    ["neb", "ts", "frequency", "irc", "endpoints", "endpoint-frequency", "all"],
)
def test_path_preflight_is_unbiased_and_preserves_combined_input_state(
    preflight_context, through
):
    context = preflight_context
    original = deepcopy(context.request.qc_params)
    coordinates = [m.coordinates.copy() for m in context.request.reactants]
    assert workflow.preflight_path_request(context.request, through, {}) == [
        "geometric",
        "xtb",
    ]
    settings, molecules = context.backend.call_args.args
    assert settings == dict(
        software="xtb",
        method=None,
        basis=None,
        nprocs=2,
        scf_cycles=100,
        opt_cycles=80,
        opt_threshold="tight",
        geometry_optimizer="geometric",
        gamma=0.0,
        opt_target="minimum",
        xtb_model="gfn2",
    )
    assert len(molecules) == 1
    merged = molecules[0]
    assert merged.atoms_list == ["He", "H"]
    assert (merged.charge, merged.multiplicity, merged.scftype) == (0, 2, "uhf")
    assert merged.fragments == [[0], [1]]
    np.testing.assert_array_equal(
        merged.coordinates, [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]]
    )
    assert through in context.backend.call_args.kwargs["check_example"]
    assert_unchanged(context, original, coordinates)
    assert [
        (m.charge, m.multiplicity, m.scftype) for m in context.request.reactants
    ] == [(0, 1, "rhf"), (0, 2, "uhf")]


def test_charged_merged_state_is_forwarded_without_changing_fragments(
    preflight_context,
):
    context = preflight_context
    hydrogen = context.request.reactants[1]
    hydrogen.charge, hydrogen.multiplicity, hydrogen.scftype = -1, 1, "rhf"
    context.request.qc_params.update(charge=-1, multiplicity=1, scftype="rhf")
    original = deepcopy(context.request.qc_params)
    coordinates = [m.coordinates.copy() for m in context.request.reactants]
    workflow.preflight_path_request(context.request, "neb", {})
    merged = context.backend.call_args.args[1][0]
    assert (merged.charge, merged.multiplicity, merged.scftype) == (-1, 1, "rhf")
    assert [(m.charge, m.multiplicity) for m in context.request.reactants] == [
        (0, 1),
        (-1, 1),
    ]
    assert_unchanged(context, original, coordinates)


def test_discovery_only_does_not_require_path_backend_or_extras(preflight_context):
    assert workflow.preflight_path_request(preflight_context.request, "react", {}) == []
    preflight_context.backend.assert_not_called()
    assert list(preflight_context.root.iterdir()) == []


@pytest.mark.parametrize(
    "through,options",
    [
        ("unsupported", {}),
        ("neb", {"images": 4}),
        ("neb", {"max_cycles": 0}),
        ("neb", {"unknown": True}),
        ("neb", {"pathway_ts_source": "unsupported"}),
        ("neb", {"reactant_reference_mode": "unsupported"}),
    ],
)
def test_invalid_path_request_is_rejected_before_dependency_checks(
    preflight_context, through, options
):
    with pytest.raises(ValueError):
        workflow.preflight_path_request(preflight_context.request, through, options)
    preflight_context.backend.assert_not_called()
    preflight_context.execution.assert_not_called()
    assert list(preflight_context.root.iterdir()) == []


@pytest.mark.parametrize(
    "through,needs_sella", [("neb", False), ("ts", True), ("all", True)]
)
@pytest.mark.parametrize("available", [False, True])
def test_sella_requirement_applies_only_to_stages_that_use_it(
    preflight_context, monkeypatch, through, needs_sella, available
):
    original_import = builtins.__import__
    calls = []

    def import_sella(name, globals=None, locals=None, fromlist=(), level=0):
        if name == "sella":
            calls.append(name)
            if not available:
                raise ImportError("synthetic missing Sella")
            return SimpleNamespace()
        return original_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr(builtins, "__import__", import_sella)
    if needs_sella and not available:
        with pytest.raises(ImportError, match="pyar-chem\\[sella\\]") as caught:
            workflow.preflight_path_request(
                preflight_context.request, through, {"ts_optimizer": "sella"}
            )
        assert isinstance(caught.value.__cause__, ImportError)
    else:
        result = workflow.preflight_path_request(
            preflight_context.request, through, {"ts_optimizer": "sella"}
        )
        assert ("sella" in result) is needs_sella
    assert len(calls) == int(needs_sella)
    assert list(preflight_context.root.iterdir()) == []


@pytest.mark.parametrize(
    "interpolation,expected",
    [("linear", []), ("idpp", ["ase"]), ("geodesic", ["ase", "geodesic"])],
)
def test_interpolation_checks_only_selected_optional_apis(
    preflight_context, monkeypatch, interpolation, expected
):
    calls = []
    monkeypatch.setattr(neb, "_ase_idpp_api", lambda: calls.append("ase"))
    monkeypatch.setattr(neb, "_geodesic_api", lambda: calls.append("geodesic"))
    options = {"interpolation": interpolation}
    original = dict(options)
    result = workflow.preflight_path_request(preflight_context.request, "neb", options)
    assert calls == expected
    assert ("ASE IDPP" in result) is (interpolation != "linear")
    assert ("geodesic-interpolate" in result) is (interpolation == "geodesic")
    assert options == original
    assert list(preflight_context.root.iterdir()) == []


@pytest.mark.parametrize(
    "error",
    [
        ValueError("invalid route"),
        ImportError("missing backend"),
        OSError("inaccessible executable"),
    ],
)
def test_backend_preflight_failure_propagates_before_optional_api_checks(
    preflight_context, monkeypatch, error
):
    preflight_context.backend.side_effect = error
    api = Mock()
    monkeypatch.setattr(neb, "_ase_idpp_api", api)
    with pytest.raises(type(error)) as caught:
        workflow.preflight_path_request(
            preflight_context.request, "neb", {"interpolation": "idpp"}
        )
    assert caught.value is error
    api.assert_not_called()
    assert list(preflight_context.root.iterdir()) == []


@pytest.mark.parametrize("api", ["_ase_idpp_api", "_geodesic_api"])
@pytest.mark.parametrize(
    "error",
    [ImportError("missing optional API"), RuntimeError("incompatible optional API")],
)
def test_selected_optional_api_failure_is_actionable_and_read_only(
    preflight_context, monkeypatch, api, error
):
    monkeypatch.setattr(neb, "_ase_idpp_api", Mock())
    monkeypatch.setattr(neb, "_geodesic_api", Mock())
    getattr(neb, api).side_effect = error
    with pytest.raises(type(error)) as caught:
        workflow.preflight_path_request(
            preflight_context.request, "neb", {"interpolation": "geodesic"}
        )
    assert caught.value is error
    preflight_context.execution.assert_not_called()
    assert list(preflight_context.root.iterdir()) == []
