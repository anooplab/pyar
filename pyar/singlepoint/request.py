"""Resolve fixed-geometry requests without inheriting optimization defaults."""

from pathlib import Path
from types import SimpleNamespace
import re
import shutil

import numpy as np

from pyar.backend_capabilities import get_backend_capabilities, normalize_backend_name
from pyar.dependencies import (
    missing_executable_message,
    missing_python_message,
    module_available,
)
from pyar.optimization_request import load_inputs
from pyar.singlepoint.native import EVALUATORS


DEFAULTS = {
    "xtb": {"xtb_model": "gfn2"},
    "orca": {"method": "PBE", "basis": "def2-SVP"},
    "gaussian": {"method": "PBE1PBE", "basis": "def2-SVP"},
    "turbomole": {"method": "pbe", "basis": "def2-SVP"},
    "aimnet_2": {},
    "ani": {"model": "ANI1x"},
}
SOLVENT_MODELS = {
    "xtb": ("alpb", "gbsa"),
    "orca": ("CPCM", "SMD"),
    "gaussian": ("PCM", "SMD"),
}


def positive(value, name):
    if not np.isfinite(value) or value <= 0:
        raise ValueError(f"--{name.replace('_', '-')} must be finite and positive")


def token(value, name):
    if value is not None and not re.fullmatch(r"[A-Za-z0-9_+*(),.\-]+", value):
        raise ValueError(
            f"--{name} must be a single backend keyword without whitespace or control characters"
        )


def theory_settings(args, backend):
    values = dict(DEFAULTS[backend])
    for name in ("method", "basis", "xtb_model", "model"):
        value = getattr(args, name)
        if value is not None:
            if name not in values:
                raise ValueError(
                    f"Backend {backend!r} does not support --{name.replace('_', '-')}"
                )
            token(value, name)
            reject_workflow_keyword(value)
            values[name] = value
    if backend == "orca":
        resolve_orca_basis(args, values)
    return values


def reject_workflow_keyword(value):
    forbidden = (
        r"(opt|freq|engrad|^sp$|^scan$|^irc$|^md$|goat|neb|cpcm|smd|pal\d|noautostart)"
    )
    if re.search(forbidden, value, re.I):
        raise ValueError(
            "Theory options must not contain workflow or execution keywords"
        )


def resolve_orca_basis(args, values):
    from pyar.backends.orca_methods import orca_method

    method, xtb = orca_method(values["method"])
    if method == "g-xTB":
        raise ValueError("Use --backend xtb --xtb-model gxtb for external g-xTB")
    builtin = method.lower() in {"r2scan-3c", "b97-3c", "pbeh-3c", "hf-3c"}
    if xtb or builtin:
        if args.basis is not None:
            raise ValueError("This ORCA method defines its own basis; omit --basis")
        values["basis"] = None
    values["method"] = method


def solvent_settings(args, backend):
    if not args.solvent:
        if args.solvent_model:
            raise ValueError("--solvent-model requires --solvent")
        return {"solvent": None, "solvent_model": None}
    token(args.solvent, "solvent")
    model = resolve_solvent_model(backend, args.solvent_model)
    return {"solvent": args.solvent, "solvent_model": model}


def resolve_solvent_model(backend, supplied):
    choices = SOLVENT_MODELS.get(backend, ())
    if not choices:
        raise ValueError(
            f"Backend {backend!r} does not support solvent settings in this workflow"
        )
    model = supplied or choices[0]
    normalized = {item.lower(): item for item in choices}
    if model.lower() not in normalized:
        raise ValueError(f"Choose a solvent model from {choices} for {backend}")
    return normalized[model.lower()]


def validate_backend_options(args, backend, settings):
    if settings.get("xtb_model") not in {None, "gfn0", "gfn1", "gfn2", "gfnff", "gxtb"}:
        raise ValueError("Choose --xtb-model gfn0, gfn1, gfn2, gfnff or gxtb")
    if backend == "ani" and settings["model"] not in {"ANI1x", "ANI1ccx", "ANI2x"}:
        raise ValueError("Choose --model ANI1x, ANI1ccx or ANI2x")
    validate_method(backend, settings)
    validate_dispersion(args, backend)
    validate_resources(args, backend)
    validate_electronic_controls(args, backend)


def validate_method(backend, settings):
    # The Gaussian reader extracts SCF Done energies. Correlated energies must
    # never be mislabeled as MP2/CC energies by accepting the SCF reference.
    gaussian = {
        "hf",
        "pbe1pbe",
        "pbepbe",
        "bp86",
        "b3lyp",
        "cam-b3lyp",
        "wb97xd",
        "m06",
        "m062x",
        "m06l",
        "tpsstpss",
        "blyp",
        "bhandhlyp",
    }
    if backend == "gaussian" and settings["method"].lower() not in gaussian:
        raise ValueError(
            "Gaussian standalone calculations currently support HF and these SCF DFT methods: "
            + ", ".join(sorted(gaussian))
        )


def validate_dispersion(args, backend):
    if args.dispersion == "none":
        return
    if backend not in {"orca", "gaussian"}:
        raise ValueError(
            "Explicit --dispersion is supported only for ORCA and Gaussian"
        )
    if backend == "gaussian" and args.dispersion == "D4":
        raise ValueError("This Gaussian adapter supports D3BJ or D3ZERO, not D4")


def validate_resources(args, backend):
    if args.memory_mb is not None and backend not in {"orca", "gaussian"}:
        raise ValueError("--memory-mb is supported only for ORCA and Gaussian")
    if backend == "turbomole" and args.nprocs != 1:
        raise ValueError(
            "The standalone Turbomole adapter currently supports --nprocs 1 only"
        )


def validate_electronic_controls(args, backend):
    if args.accuracy is not None and backend != "xtb":
        raise ValueError("--accuracy is supported only for xTB")
    if backend in {"aimnet_2", "ani"} and args.scf_cycles is not None:
        raise ValueError("ML potentials do not accept --scf-cycles")


RESOURCE_DEFAULTS = {
    "xtb": {"accuracy": 0.1, "scf_cycles": 250},
    "orca": {"memory_mb": 2000, "scf_cycles": 250},
    "gaussian": {"memory_mb": 2000, "scf_cycles": 250},
    "turbomole": {"scf_cycles": 250},
    "ani": {},
    "aimnet_2": {},
}


def common_settings(args, backend):
    values = {
        "software": backend,
        "nprocs": args.nprocs,
        "model": None,
        "memory_mb": None,
        "scf_cycles": None,
        "accuracy": None,
        "dispersion": args.dispersion,
    }
    for name, default in RESOURCE_DEFAULTS[backend].items():
        supplied = getattr(args, name)
        values[name] = default if supplied is None else supplied
    return values


def validate_numbers(args, workflow):
    for name in ("nprocs", "memory_mb", "scf_cycles", "accuracy"):
        value = getattr(args, name)
        if value is not None:
            positive(value, name)
    if workflow == "frequency":
        for name in ("step", "max_gradient", "rms_gradient", "imaginary_threshold"):
            positive(getattr(args, name), name)


def resolve(args, workflow):
    backend = normalize_backend_name(args.backend.lower())
    if backend not in EVALUATORS:
        raise ValueError(
            f"Unsupported fixed-geometry backend {backend!r}; choose {', '.join(EVALUATORS)}"
        )
    validate_numbers(args, workflow)
    settings = {
        **common_settings(args, backend),
        **theory_settings(args, backend),
        **solvent_settings(args, backend),
    }
    validate_backend_options(args, backend, settings)
    molecules = load_inputs(SimpleNamespace(**vars(args), scftype=None))
    validate_states(backend, molecules)
    root = Path(args.output).expanduser().resolve()
    if root.exists():
        raise ValueError(
            f"Output already exists: {root}; choose a new --output directory"
        )
    return settings, molecules, root


def validate_states(backend, molecules):
    for molecule in molecules:
        validate_name(molecule.name)
        validate_model_state(backend, molecule)


def validate_name(name):
    if not re.fullmatch(r"[\w.\-]+", name) or name in {".", ".."}:
        raise ValueError(f"Unsafe calculation basename {name!r}; rename the input")


def validate_model_state(backend, molecule):
    if backend in {"ani", "aimnet_2"} and molecule.multiplicity != 1:
        raise ValueError(f"{backend} adapter currently supports singlets only")
    if backend == "ani" and molecule.charge != 0:
        raise ValueError("ANI models support neutral molecules only")


def preflight(settings):
    backend = settings["software"]
    requirements = executable_requirements(backend) + module_requirements(backend)
    if backend == "aimnet_2":
        from pyar.backends.aimnet2_assets import validate_aimnet2_runtime_assets

        validate_aimnet2_runtime_assets()
    if backend == "ani":
        validate_ani_model(settings["model"])
    if settings.get("xtb_model") == "gxtb":
        from pyar.backends.xtb import xtb_supports_gxtb

        if not xtb_supports_gxtb(shutil.which("xtb")):
            raise ValueError("Installed xTB does not advertise --gxtb support")
    return requirements


def executable_requirements(backend):
    requirements = sorted(get_backend_capabilities(backend).single_point_executables)
    for executable in requirements:
        if not shutil.which(executable):
            raise ValueError(missing_executable_message(executable))
    return requirements


def module_requirements(backend):
    requirements = sorted(get_backend_capabilities(backend).required_python_modules)
    if backend == "turbomole":
        requirements.append("ase")
    for module in requirements:
        if not module_available(module):
            raise ValueError(missing_python_message(module))
    return requirements


def validate_ani_model(model):
    import torchani

    if not hasattr(torchani.models, model):
        raise ValueError(
            f"Installed TorchANI does not provide {model}; use ANI1x or ANI1ccx"
        )
