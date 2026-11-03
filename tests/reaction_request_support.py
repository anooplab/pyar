"""Stable request evidence; no chemistry backend or invented reference energies."""

import hashlib
import json

CASES = {
    "adaptive": {},
    "none_defaults": {
        "orientations": None,
        "proximity_factor": None,
        "softmin_beta": None,
    },
    "fixed": {"bias_controller": "fixed", "bias_min": 10},
    "scheduled": {
        "bias_controller": "scheduled",
        "bias_min": 10,
        "bias_scheduled_alpha": 0.002,
    },
    "charged_softmin": {
        "charge": [-1, 0],
        "multiplicity": [2, 2],
        "site": [1, 0],
        "bias_potential": "softmin",
        "softmin_beta": 2.0,
    },
    "orca_dft": {"backend": "orca16", "method": "BP86", "basis": "def2-SVP"},
    "orca_xtb": {"backend": "orca", "method": "XTB2"},
    # Preserve existing validation semantics; bool settings are not scientifically certified.
    "existing_boolean_acceptance": {
        "orientations": True,
        "nprocs": True,
        "release_retry_limit": True,
    },
}


def snapshot(request, root):
    restart = request.restart_request
    digest = hashlib.sha256(
        json.dumps(restart, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    result = dict(
        qc_params=request.qc_params,
        bias_min=request.bias_min,
        bias_max=request.bias_max,
        orientations=request.orientations,
        site=request.site,
        proximity_factor=request.proximity_factor,
        restart_request=restart,
        restart_request_sha256=digest,
        reactants=[
            dict(
                name=m.name,
                atoms=list(m.atoms_list),
                coordinates=m.coordinates.tolist(),
                charge=m.charge,
                multiplicity=m.multiplicity,
                scftype=m.scftype,
            )
            for m in request.reactants
        ],
    )

    def normalize(value):
        if isinstance(value, dict):
            return {key: normalize(item) for key, item in value.items()}
        if isinstance(value, list):
            return [normalize(item) for item in value]
        if isinstance(value, str):
            return value.replace(str(root), "<ROOT>")
        return value

    return normalize(result)
