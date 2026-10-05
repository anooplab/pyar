"""Validated fixed-seed, repeated-addend growth requests."""
from dataclasses import dataclass
import numpy as np

from pyar.aggregation.request import AggregateRequest, _molecule_signature
from pyar.backend_capabilities import (normalize_backend_name, validate_backend_capability,
                                       BACKEND_CAPABILITIES, get_backend_capabilities, unsupported_qc_options)


@dataclass(frozen=True)
class GrowRequest:
    seed: object
    monomer: object
    count: int
    number_of_orientations: int = 8
    backend_parameters: dict | None = None
    maximum_number_of_seeds: int = 12
    site: tuple[int, int] | None = None
    connectivity_policy: str = "auto"
    selection_feature: str = "auto"
    selection_algorithm: str = "auto"
    selection_distance: str = "euclidean"
    selection_system_type: str = "auto"

    def __post_init__(self):
        for label, value in (("count", self.count), ("orientations", self.number_of_orientations),
                             ("maximum_number_of_seeds", self.maximum_number_of_seeds)):
            if isinstance(value, bool) or not isinstance(value, int) or value < 1:
                raise ValueError(f"{label} must be a positive integer")
        for molecule in (self.seed, self.monomer):
            if not hasattr(molecule, "atoms_list"):
                raise ValueError("Grow requires exactly one seed molecule and one monomer molecule")
            if not molecule.atoms_list or not np.isfinite(molecule.coordinates).all():
                raise ValueError("Grow requires valid seed and monomer geometries")
            electrons = sum(molecule.atomic_number) - molecule.charge
            if (electrons < 0 or molecule.multiplicity < 1
                    or (electrons - (molecule.multiplicity - 1)) % 2
                    or molecule.multiplicity - 1 > electrons):
                raise ValueError(f"Invalid charge/multiplicity for {molecule.name}")
        if self.site is not None:
            if len(self.site) != 2:
                raise ValueError("Site requires one local atom index per input")
            for index, molecule in zip(self.site, (self.seed, self.monomer)):
                if not isinstance(index, int) or not 0 <= index < len(molecule):
                    raise ValueError(f"Site index {index} is out of range for {molecule.name}")
        qc = dict(self.backend_parameters or {})
        if qc.get("software"):
            from pyar.data.defualt_parameters import values
            for key in ("nprocs", "opt_cycles", "scf_cycles", "opt_threshold", "scf_threshold", "geometry_optimizer"):
                qc.setdefault(key, values[key])
        if not qc.get("software") and any(qc.get(k) is not None for k in ("method", "basis", "model", "xtb_model")):
            raise ValueError("Backend options require software")
        if qc.get("software"):
            qc["software"] = normalize_backend_name(qc["software"])
            if qc["software"] not in BACKEND_CAPABILITIES:
                raise ValueError(f"Unknown backend: {qc['software']}")
            caps = get_backend_capabilities(qc["software"])
            for molecule in (self.seed, self.monomer):
                if (molecule.charge and not caps.supports_charge) or (molecule.multiplicity != 1 and not caps.supports_multiplicity):
                    raise ValueError("Backend does not support the requested electronic state")
            unsupported = unsupported_qc_options(qc["software"], {k for k in ("method", "basis", "model") if qc.get(k) is not None})
            if unsupported:
                raise ValueError(f"Unsupported backend options: {', '.join(sorted(unsupported))}")
            if qc["software"] == "orca":
                from pyar.backends.orca_methods import orca_method_keywords, orca_external_method_block
                orca_method_keywords(qc)
                orca_external_method_block(qc)
            if qc["software"] == "gaussian" and not all(qc.get(k) for k in ("method", "basis")):
                raise ValueError("Gaussian requires method and basis")
            if qc.get("xtb_model") is not None:
                from pyar.backends.xtb_utils import canonical_xtb_model
                if qc["software"] != "xtb":
                    raise ValueError("xtb_model is only valid for standalone xtb")
                qc["xtb_model"] = canonical_xtb_model(qc["xtb_model"])
            if qc.get("gxtb_wrapper") is not None:
                import os
                from pathlib import Path
                wrapper = Path(qc["gxtb_wrapper"]).expanduser().resolve()
                if not wrapper.is_file() or not os.access(wrapper, os.X_OK):
                    raise ValueError("gxtb_wrapper must be an executable file")
                qc["gxtb_wrapper"] = str(wrapper)
            optimizer = qc.get("geometry_optimizer", "native")
            if optimizer not in {"native", "geometric"}:
                raise ValueError("Geometry optimizer must be native or geometric")
            validate_backend_capability(qc["software"],
                                        ["native_optimization" if optimizer == "native" else "energy_gradient"],
                                        context="grow")
        if qc.get("opt_target", "minimum") not in {"min", "minimum"} or qc.get("gamma") not in (None, 0):
            raise ValueError("Grow supports unbiased minimum optimization only")
        if qc.get("opt_target") == "min":
            qc["opt_target"] = "minimum"
        for option in ("nprocs", "opt_cycles", "scf_cycles"):
            if option in qc and int(qc[option]) < 1:
                raise ValueError(f"{option} must be positive")
        from pyar.molecule_merge import combine_multiplicity
        charge, multiplicity = self.seed.charge, self.seed.multiplicity
        electrons = sum(self.seed.atomic_number)
        for _ in range(self.count):
            charge += self.monomer.charge
            multiplicity = combine_multiplicity(multiplicity, self.monomer.multiplicity)
            electrons += sum(self.monomer.atomic_number)
            if electrons - charge < multiplicity - 1 or (electrons - charge - multiplicity + 1) % 2:
                raise ValueError("The existing fragment spin-combination rule yields an invalid growth state")
            if qc.get("software"):
                if (charge and not caps.supports_charge) or (multiplicity != 1 and not caps.supports_multiplicity):
                    raise ValueError("Backend does not support a combined growth electronic state")
        normalized = AggregateRequest.from_options(
            [self.seed, self.monomer], [1, self.count], self.number_of_orientations,
            qc, self.maximum_number_of_seeds, 0, 1, self.site, self.connectivity_policy,
            self.selection_feature, self.selection_algorithm, self.selection_distance,
            self.selection_system_type)
        object.__setattr__(self, "backend_parameters", qc)
        for key in ("connectivity_policy", "selection_feature", "selection_algorithm",
                    "selection_distance", "selection_system_type"):
            object.__setattr__(self, key, getattr(normalized, key))

    def to_state_dict(self):
        from pyar.growth.service import sampling_configuration
        from pyar.selection.policy import SELECTION_POLICY_VERSION
        return {"seed": _molecule_signature(self.seed), "monomer": _molecule_signature(self.monomer),
                "input_growth_metadata": [
                    {key: getattr(m, key) for key in ("connectivity_policy_hint", "growth_kind", "source_kind") if hasattr(m, key)}
                    for m in (self.seed, self.monomer)],
                "count": self.count, "orientations": self.number_of_orientations,
                "backend_parameters": self.backend_parameters,
                "maximum_number_of_seeds": self.maximum_number_of_seeds,
                "site": None if self.site is None else list(self.site),
                "connectivity_policy": self.connectivity_policy,
                "selection_feature": self.selection_feature,
                "selection_algorithm": self.selection_algorithm,
                "selection_distance": self.selection_distance,
                "selection_system_type": self.selection_system_type,
                "selection_policy_version": SELECTION_POLICY_VERSION,
                "sampling": sampling_configuration(number_of_orientations=self.number_of_orientations,
                                                    use_angles=len(self.monomer) > 1)}
