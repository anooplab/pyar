"""Input fixtures and exact serialization fingerprints from the original request."""

from copy import deepcopy
import hashlib
import json

from pyar.core.molecule import Molecule

CASES = {
    "default": {},
    "sites": {
        "site": ["1", "2"],
        "site_groups": [(" hydrogens ", ["1", "2"])],
        "contact_distance": 3.2,
    },
    "energy_shell": {
        "placement_score": "energy",
        "backend_parameters": {"software": "xtb"},
    },
    "explicit_wall": {
        "wall": "outer-sphere",
        "wall_radius": 5.0,
        "wall_force_constant": 0.5,
        "backend_parameters": {"software": "xtb", "geometry_optimizer": "geometric"},
    },
    "automatic_shell": {
        "wall": "auto",
        "backend_parameters": {"software": "xtb", "geometry_optimizer": "geometric"},
    },
    "cluster": {
        "objective": "cluster",
        "placement_score": "energy",
        "backend_parameters": {"software": "xtb"},
    },
    "automatic_cluster": {
        "wall": "auto",
        "objective": "cluster",
        "placement_score": "energy",
        "backend_parameters": {"software": "xtb", "geometry_optimizer": "geometric"},
    },
    "resampling": {
        "objective": "cluster",
        "placement_score": "energy",
        "resampling_moves": 2,
        "resampling_evaluation_budget": 20,
        "resampling_time_limit": 30.0,
        "resampling_seed": 12,
        "backend_parameters": {"software": "xtb", "geometry_optimizer": "geometric"},
    },
    "backend_alias": {"backend_parameters": {"software": "orca16"}},
}


def water(name):
    return Molecule(
        ["O", "H", "H"],
        [[0, 0, 0], [0.9572, 0, 0], [-0.239, 0.927, 0]],
        name=name,
        charge=0,
        multiplicity=1,
    )


def request_values(case):
    return dict(
        solute=water("solute"),
        solvent=water("solvent"),
        count=2,
        **deepcopy(CASES[case]),
    )


def fingerprint(request):
    state = request.to_state_dict()
    return dict(
        canonical_sha256=hashlib.sha256(
            json.dumps(state, sort_keys=True, separators=(",", ":")).encode()
        ).hexdigest(),
        ordered_sha256=hashlib.sha256(
            json.dumps(state, separators=(",", ":")).encode()
        ).hexdigest(),
        site=None if request.site is None else list(request.site),
        site_groups=[[name, list(indices)] for name, indices in request.site_groups],
    )
