"""Typed request model for conformer-search workflows."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping


class ConformerRequestError(RuntimeError):
    """Raised when a conformer-search request is invalid."""


def _normalize_torsion_mode(torsion_mode):
    """Return the canonical torsion mode or raise a request error."""
    normalized = str(torsion_mode).lower().replace("-", "_")
    if normalized == "torsion_kick":
        return "random"
    if normalized != "random":
        raise ConformerRequestError("--torsion-mode must be 'random'")
    return normalized


@dataclass(frozen=True)
class ConformerRequest:
    """Validated options for one conformer-search run."""

    input: str
    input_format: str
    num_conformers: int
    top_n: int
    backend_top_n: int | None
    num_seeds: int
    diversity_fraction: float
    compactness_fraction: float
    rms_threshold: float
    generation_dedup_rms: float
    use_random_coords: bool
    torsion_kicks: bool
    torsion_rounds: int
    torsion_mode: str
    torsion_kicks_per_conformer: int
    torsion_max_bonds: int
    torsion_dedup_rms: float
    force_field: str
    seed: int
    seed_values: tuple[int, ...]
    num_threads: int
    max_iterations: int
    charge: int | None
    multiplicity: int
    scftype: str
    backend_parameters: Mapping[str, Any] = field(default_factory=dict)

    @classmethod
    def from_options(
        cls,
        input_spec,
        *,
        input_format,
        num_conformers,
        top_n,
        backend_top_n,
        num_seeds,
        diversity_fraction,
        compactness_fraction,
        rms_threshold,
        use_random_coords,
        torsion_kicks,
        torsion_rounds,
        torsion_mode,
        torsion_kicks_per_conformer,
        torsion_max_bonds,
        torsion_dedup_rms,
        force_field,
        seed,
        num_threads,
        max_iterations,
        charge,
        multiplicity,
        scftype,
        qc_params,
    ):
        """Build a normalized request from public workflow arguments."""
        num_conformers = int(num_conformers)
        top_n = int(top_n)
        num_seeds = int(num_seeds)
        rms_threshold = float(rms_threshold)
        torsion_rounds = int(torsion_rounds)
        torsion_kicks_per_conformer = int(torsion_kicks_per_conformer)
        torsion_max_bonds = int(torsion_max_bonds)
        torsion_dedup_rms = float(torsion_dedup_rms)
        diversity_fraction = float(diversity_fraction)
        compactness_fraction = float(compactness_fraction)

        if num_conformers < 1:
            raise ConformerRequestError("--num-conformers must be at least 1")
        if top_n < 1:
            raise ConformerRequestError("--top-n must be at least 1")
        if num_seeds < 1:
            raise ConformerRequestError("--num-seeds must be at least 1")
        if rms_threshold < 0.0:
            raise ConformerRequestError("--rms-threshold must be non-negative")
        if torsion_rounds < 0:
            raise ConformerRequestError("--torsion-rounds must be non-negative")
        if torsion_kicks_per_conformer < 0:
            raise ConformerRequestError("--torsion-kicks-per-conformer must be non-negative")
        if torsion_max_bonds < 1:
            raise ConformerRequestError("--torsion-max-bonds must be at least 1")
        if torsion_dedup_rms < 0.0:
            raise ConformerRequestError("--torsion-dedup-rms must be non-negative")
        if backend_top_n is not None and int(backend_top_n) < 1:
            raise ConformerRequestError("--backend-top-n must be at least 1")
        if not 0.0 <= diversity_fraction <= 1.0:
            raise ConformerRequestError("--diversity-fraction must be between 0 and 1")
        if not 0.0 <= compactness_fraction <= 1.0:
            raise ConformerRequestError("--compactness-fraction must be between 0 and 1")

        seed = int(seed)
        seed_values = tuple(seed + index for index in range(num_seeds))
        generation_dedup_rms = max(torsion_dedup_rms, rms_threshold, 0.5)
        return cls(
            input=str(input_spec),
            input_format=input_format,
            num_conformers=num_conformers,
            top_n=top_n,
            backend_top_n=None if backend_top_n is None else int(backend_top_n),
            num_seeds=num_seeds,
            diversity_fraction=diversity_fraction,
            compactness_fraction=compactness_fraction,
            rms_threshold=rms_threshold,
            generation_dedup_rms=float(generation_dedup_rms),
            use_random_coords=bool(use_random_coords),
            torsion_kicks=bool(torsion_kicks),
            torsion_rounds=torsion_rounds,
            torsion_mode=_normalize_torsion_mode(torsion_mode),
            torsion_kicks_per_conformer=torsion_kicks_per_conformer,
            torsion_max_bonds=torsion_max_bonds,
            torsion_dedup_rms=torsion_dedup_rms,
            force_field=force_field,
            seed=seed,
            seed_values=seed_values,
            num_threads=int(num_threads),
            max_iterations=int(max_iterations),
            charge=charge,
            multiplicity=int(multiplicity),
            scftype=scftype,
            backend_parameters=dict(qc_params or {}),
        )

    @property
    def backend_requested(self):
        """Return whether the request asks for backend refinement."""
        return bool(self.backend_parameters.get("software"))

    def to_state_dict(self):
        """Return the JSON-serializable representation stored in state files."""
        return {
            "input": self.input,
            "input_format": self.input_format,
            "num_conformers": self.num_conformers,
            "top_n": self.top_n,
            "backend_top_n": self.backend_top_n,
            "num_seeds": self.num_seeds,
            "diversity_fraction": self.diversity_fraction,
            "compactness_fraction": self.compactness_fraction,
            "rms_threshold": self.rms_threshold,
            "generation_dedup_rms": self.generation_dedup_rms,
            "use_random_coords": self.use_random_coords,
            "torsion_kicks": self.torsion_kicks,
            "torsion_rounds": self.torsion_rounds,
            "torsion_mode": self.torsion_mode,
            "torsion_kicks_per_conformer": self.torsion_kicks_per_conformer,
            "torsion_max_bonds": self.torsion_max_bonds,
            "torsion_dedup_rms": self.torsion_dedup_rms,
            "force_field": self.force_field,
            "seed": self.seed,
            "seed_values": list(self.seed_values),
            "num_threads": self.num_threads,
            "max_iterations": self.max_iterations,
            "charge": self.charge,
            "multiplicity": self.multiplicity,
            "scftype": self.scftype,
            "backend_parameters": dict(self.backend_parameters),
        }
