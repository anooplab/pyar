"""Adapter for the external RGD1-TSopt-GFN2 dataset."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import random
from collections import defaultdict
from pathlib import Path

from pyar.benchmarks.ts_optimizer import (
    TSOptimizerBenchmarkError,
    load_ts_optimizer_benchmark,
)
from pyar.neb import read_xyz


DATASET_NAME = "RGD1-TSopt-GFN2"
DATASET_VERSION = "1.0"
DATASET_DOI = "10.5281/zenodo.20489312"
DATASET_URL = "https://zenodo.org/records/20489312"
ARCHIVE_MD5 = "f07008a7cc6f0e3441acb439a9ab6717"
ARCHIVE_SHA256 = "f97cea45e78c6a93fd5439f0bd6c6f99ce3db1ceba3491f0dd4b2cdc876394c0"
TIERS = ("easy", "med", "hard")
TIER_ALPHA_ANGSTROM = {"easy": 0.06, "med": 0.11, "hard": 0.15}
ATOM_BINS = ((4, 8), (9, 12), (13, 16), (17, 24))
MANIFEST_FIELDS = {
    "id", "tier", "alpha_react_A", "beta_nonreact_A", "total_rmsd_A",
    "react_overlap", "n_atoms", "charge",
}


def _sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_dataset(dataset_dir):
    root = Path(dataset_dir).expanduser().resolve()
    manifest_path = root / "manifest.tsv"
    settings_path = root / "settings.json"
    if not manifest_path.is_file() or not settings_path.is_file():
        raise TSOptimizerBenchmarkError(
            "dataset directory must contain manifest.tsv and settings.json"
        )
    try:
        settings = json.loads(settings_path.read_text(encoding="utf-8"))
        with manifest_path.open(newline="", encoding="utf-8") as stream:
            reader = csv.DictReader(stream, delimiter="\t")
            if not reader.fieldnames or set(reader.fieldnames) != MANIFEST_FIELDS:
                raise TSOptimizerBenchmarkError(
                    "RGD1 manifest.tsv has an unexpected header"
                )
            rows = list(reader)
    except (OSError, json.JSONDecodeError, csv.Error) as exc:
        raise TSOptimizerBenchmarkError(f"cannot read RGD1 dataset metadata: {exc}") from exc

    if not isinstance(settings, dict) or settings.get("method") != "GFN2-xTB":
        raise TSOptimizerBenchmarkError("dataset settings.json must describe GFN2-xTB")
    if (settings.get("charge") != 0 or settings.get("unpaired_electrons") != 0
            or settings.get("implementation") != "tblite 0.6.0"
            or settings.get("electronic_temperature_K") != 300
            or settings.get("scf_accuracy") != 0.01
            or settings.get("solvation") != "none (gas phase)"
            or not isinstance(settings.get("units"), dict)
            or settings["units"].get("xyz") != "Angstrom"
            or settings["units"].get("energy") != "Hartree"):
        raise TSOptimizerBenchmarkError(
            "settings.json differs from the published GFN2-xTB, neutral-singlet reference protocol"
        )
    grouped = defaultdict(dict)
    atom_counts = {}
    charges = {}
    for line_number, row in enumerate(rows, start=2):
        case_id, tier = row.get("id", ""), row.get("tier", "")
        if not case_id or any(
                c not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_.-"
                for c in case_id):
            raise TSOptimizerBenchmarkError(
                f"invalid RGD1 case identifier at manifest line {line_number}"
            )
        if tier not in TIERS:
            raise TSOptimizerBenchmarkError(
                f"unsupported RGD1 tier {tier!r} at manifest line {line_number}"
            )
        if tier in grouped[case_id]:
            raise TSOptimizerBenchmarkError(f"duplicate RGD1 row: {case_id}/{tier}")
        try:
            n_atoms = int(row["n_atoms"])
            charge = int(row["charge"])
            alpha = float(row["alpha_react_A"])
            beta = float(row["beta_nonreact_A"])
            total_rmsd = float(row["total_rmsd_A"])
            overlap = float(row["react_overlap"])
        except (TypeError, ValueError, OverflowError) as exc:
            raise TSOptimizerBenchmarkError(
                f"invalid numeric value in RGD1 manifest line {line_number}"
            ) from exc
        if (n_atoms < 1 or charge != 0 or not all(
                math.isfinite(value) for value in (alpha, beta, total_rmsd, overlap))):
            raise TSOptimizerBenchmarkError(
                f"invalid geometry metadata in RGD1 manifest line {line_number}"
            )
        if (not math.isclose(alpha, TIER_ALPHA_ANGSTROM[tier], abs_tol=1e-9)
                or not math.isclose(beta, 0.12, abs_tol=1e-9)):
            raise TSOptimizerBenchmarkError(
                f"unexpected displacement parameters for {case_id}/{tier}"
            )
        if case_id in atom_counts and (
                atom_counts[case_id] != n_atoms or charges[case_id] != charge):
            raise TSOptimizerBenchmarkError(f"inconsistent metadata for RGD1 case {case_id}")
        atom_counts[case_id], charges[case_id] = n_atoms, charge
        grouped[case_id][tier] = {
            **row,
            "n_atoms": n_atoms,
            "charge": charge,
            "alpha_react_A": alpha,
            "beta_nonreact_A": beta,
            "total_rmsd_A": total_rmsd,
            "react_overlap": overlap,
        }
    if not grouped:
        raise TSOptimizerBenchmarkError("RGD1 manifest.tsv contains no cases")
    for case_id, tiers in grouped.items():
        if set(tiers) != set(TIERS):
            raise TSOptimizerBenchmarkError(
                f"RGD1 case {case_id} must contain each of easy, med, and hard"
            )
        case_root = root / "cases" / case_id
        required_files = [case_root / name for name in (
            "ts_ref.xyz", "ts_ref.energy", "reactant.xyz", "product.xyz",
            "ts_ref.hessian", "provenance.txt",
            *(f"start_{tier}.xyz" for tier in TIERS),
        )]
        missing = [str(path) for path in required_files if not path.is_file()]
        if missing:
            raise TSOptimizerBenchmarkError(
                f"RGD1 case {case_id} is missing required files: " + ", ".join(missing)
            )
        symbols, _, _ = read_xyz(case_root / "ts_ref.xyz")
        if len(symbols) != atom_counts[case_id]:
            raise TSOptimizerBenchmarkError(
                f"RGD1 case {case_id} atom count does not match manifest.tsv"
            )
        if not set(symbols).issubset({"C", "H", "N", "O"}):
            raise TSOptimizerBenchmarkError(
                f"RGD1 case {case_id} contains elements outside the published C/H/N/O scope"
            )
        for geometry in required_files:
            if geometry.suffix != ".xyz":
                continue
            other_symbols, _, _ = read_xyz(geometry)
            if other_symbols != symbols:
                raise TSOptimizerBenchmarkError(
                    f"RGD1 case {case_id} {geometry.name} has a different ordered atom list"
                )
        try:
            reference_energy = float((case_root / "ts_ref.energy").read_text().split()[0])
        except (OSError, IndexError, ValueError) as exc:
            raise TSOptimizerBenchmarkError(
                f"RGD1 case {case_id} has an invalid ts_ref.energy"
            ) from exc
        if not math.isfinite(reference_energy):
            raise TSOptimizerBenchmarkError(
                f"RGD1 case {case_id} has a non-finite reference energy"
            )
        for row in tiers.values():
            row["reference_energy_hartree"] = reference_energy
    return root, settings, grouped, atom_counts


def _stratified_reactions(atom_counts, count, seed):
    ids_by_bin = {bounds: [] for bounds in ATOM_BINS}
    for case_id, n_atoms in atom_counts.items():
        for low, high in ATOM_BINS:
            if low <= n_atoms <= high:
                ids_by_bin[(low, high)].append(case_id)
                break
        else:
            raise TSOptimizerBenchmarkError(
                f"RGD1 case {case_id} has unsupported atom count {n_atoms}"
            )
    count = min(count, len(atom_counts))
    quotas = {bounds: count // len(ATOM_BINS) for bounds in ATOM_BINS}
    for bounds in ATOM_BINS[:count % len(ATOM_BINS)]:
        quotas[bounds] += 1
    rng = random.Random(seed)
    selected = []
    remaining = count
    for bounds in ATOM_BINS:
        pool = sorted(ids_by_bin[bounds])
        quota = min(quotas[bounds], len(pool))
        selected.extend(rng.sample(pool, quota))
        remaining -= quota
    if remaining:
        unselected = sorted(set(atom_counts) - set(selected))
        selected.extend(rng.sample(unselected, remaining))
    return selected


def prepare_rgd1_tsopt_manifest(
    dataset_dir, output, *, reactions=15, seed=20261002, tiers=TIERS,
    ts_fmax=0.02, ts_max_cycles=200,
):
    """Create a deterministic, atom-count-stratified RGD1 benchmark manifest."""
    root, dataset_settings, grouped, atom_counts = _read_dataset(dataset_dir)
    if isinstance(reactions, bool) or not isinstance(reactions, int) or reactions < 1:
        raise TSOptimizerBenchmarkError("reactions must be a positive integer")
    if isinstance(seed, bool) or not isinstance(seed, int):
        raise TSOptimizerBenchmarkError("seed must be an integer")
    selected_tiers = tuple(tiers)
    if (not selected_tiers or len(set(selected_tiers)) != len(selected_tiers)
            or any(tier not in TIERS for tier in selected_tiers)):
        raise TSOptimizerBenchmarkError(
            "tiers must be a unique non-empty subset of easy, med, hard"
        )
    if isinstance(ts_max_cycles, bool) or not isinstance(ts_max_cycles, int) or ts_max_cycles < 1:
        raise TSOptimizerBenchmarkError("ts_max_cycles must be a positive integer")
    try:
        ts_fmax = float(ts_fmax)
    except (TypeError, ValueError, OverflowError):
        raise TSOptimizerBenchmarkError("ts_fmax must be positive and finite") from None
    if not math.isfinite(ts_fmax) or ts_fmax <= 0:
        raise TSOptimizerBenchmarkError("ts_fmax must be positive and finite")

    selected_ids = _stratified_reactions(atom_counts, reactions, seed)
    manifest_path = Path(output).expanduser().resolve()
    if manifest_path.exists():
        raise TSOptimizerBenchmarkError(
            f"output already exists: {manifest_path}; choose a new path to preserve prior data"
        )
    cases = []
    for reaction_id in selected_ids:
        case_root = root / "cases" / reaction_id
        for tier in selected_tiers:
            row = grouped[reaction_id][tier]
            cases.append({
                "id": f"{reaction_id}_{tier}",
                "difficulty": tier,
                "ts_guess": Path(os.path.relpath(
                    case_root / f"start_{tier}.xyz", manifest_path.parent,
                )).as_posix(),
                "reference_ts": Path(os.path.relpath(
                    case_root / "ts_ref.xyz", manifest_path.parent,
                )).as_posix(),
                "reactant": Path(os.path.relpath(
                    case_root / "reactant.xyz", manifest_path.parent,
                )).as_posix(),
                "product": Path(os.path.relpath(
                    case_root / "product.xyz", manifest_path.parent,
                )).as_posix(),
                "charge": 0,
                "multiplicity": 1,
                "reference_energy_hartree": row["reference_energy_hartree"],
                "perturbation": {
                    "method": "GFN2-xTB mass-weighted Hessian normal modes",
                    "tier": tier,
                    "alpha_react_A": row["alpha_react_A"],
                    "beta_nonreact_A": row["beta_nonreact_A"],
                    "total_rmsd_A": row["total_rmsd_A"],
                    "react_overlap": row["react_overlap"],
                },
            })
    pilot = {
        "kind": "atom_count_stratified_pilot",
        "seed": seed,
        "selected_reactions": selected_ids,
        "reaction_count": len(selected_ids),
        "tiers": list(selected_tiers),
        "experimental_unit": "reaction_id + difficulty tier + starting geometry",
    }
    manifest = {
        "name": f"RGD1-TSopt-GFN2-pilot-{len(selected_ids)}-reactions",
        "source": {
            "name": DATASET_NAME,
            "version": DATASET_VERSION,
            "license": "Copyright; the Zenodo record states no redistribution license",
            "doi": DATASET_DOI,
            "url": DATASET_URL,
            "zenodo_archive_md5": ARCHIVE_MD5,
            "archive_sha256_reference": ARCHIVE_SHA256,
            "dataset_reaction_count": len(grouped),
            "dataset_starting_geometry_count": sum(len(tiers) for tiers in grouped.values()),
            "dataset_manifest_sha256": _sha256(root / "manifest.tsv"),
            "settings_sha256": _sha256(root / "settings.json"),
            "reference_implementation": dataset_settings["implementation"],
            "reference_settings": dataset_settings,
            "pilot_selection": pilot,
        },
        "qc_model": {"software": "xtb", "xtb_model": "gfn2", "nprocs": 1},
        "settings": {
            "ts_fmax": ts_fmax,
            "ts_max_cycles": ts_max_cycles,
            "imaginary_frequency_threshold": 20.0,
            "irc_max_cycles": 200,
            "endpoint_max_cycles": 300,
        },
        "cases": cases,
    }
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        with manifest_path.open("x", encoding="utf-8") as stream:
            json.dump(manifest, stream, indent=2, sort_keys=True)
            stream.write("\n")
        load_ts_optimizer_benchmark(manifest_path)
    except Exception:
        manifest_path.unlink(missing_ok=True)
        raise
    return manifest_path
