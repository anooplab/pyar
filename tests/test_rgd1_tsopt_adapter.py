import csv
import json
from pathlib import Path

import pytest

from pyar.benchmarks import rgd1_tsopt as adapter
from pyar.benchmarks.ts_optimizer import load_ts_optimizer_benchmark
from pyar.scripts.ts_optimizer_benchmark import argument_parse


def _dataset(tmp_path, atom_counts=(4, 9, 13, 17)):
    root = tmp_path / "RGD1-TSopt-GFN2"
    cases_root = root / "cases"
    cases_root.mkdir(parents=True)
    fields = [
        "id", "tier", "alpha_react_A", "beta_nonreact_A", "total_rmsd_A",
        "react_overlap", "n_atoms", "charge",
    ]
    with (root / "manifest.tsv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for reaction_index, atom_count in enumerate(atom_counts):
            case_id = f"MR_sample_{reaction_index}"
            case_root = cases_root / case_id
            case_root.mkdir(parents=True)
            symbols = ["C"] + ["H"] * (atom_count - 1)
            geometry = "\n".join(
                [str(atom_count), case_id]
                + [f"{symbol} {index * 0.4:.5f} 0.0 0.0" for index, symbol in enumerate(symbols)]
            ) + "\n"
            (case_root / "ts_ref.xyz").write_text(geometry)
            (case_root / "ts_ref.hessian").write_text("synthetic Hessian fixture\n")
            (case_root / "provenance.txt").write_text(f"source=RGD1 reaction {case_id}\n")
            (case_root / "reactant.xyz").write_text(geometry)
            (case_root / "product.xyz").write_text(geometry)
            (case_root / "ts_ref.energy").write_text("-1.2345\n")
            for tier, alpha in zip(adapter.TIERS, (0.06, 0.11, 0.15)):
                (case_root / f"start_{tier}.xyz").write_text(geometry)
                writer.writerow({
                    "id": case_id,
                    "tier": tier,
                    "alpha_react_A": alpha,
                    "beta_nonreact_A": 0.12,
                    "total_rmsd_A": 0.14,
                    "react_overlap": 0.447,
                    "n_atoms": atom_count,
                    "charge": 0,
                })
    (root / "settings.json").write_text(json.dumps({
        "method": "GFN2-xTB", "implementation": "tblite 0.6.0", "charge": 0,
        "unpaired_electrons": 0, "electronic_temperature_K": 300,
        "scf_accuracy": 0.01, "solvation": "none (gas phase)",
        "units": {"xyz": "Angstrom", "energy": "Hartree"},
    }))
    (root / "README.md").write_text("synthetic test dataset")
    return root


def test_rgd1_adapter_creates_deterministic_stratified_pilot_manifest(tmp_path):
    dataset = _dataset(tmp_path, atom_counts=(4, 4, 9, 9, 13, 13, 17, 17))
    first = adapter.prepare_rgd1_tsopt_manifest(
        dataset, dataset / "pilot.json", reactions=4, seed=73,
    )
    raw = json.loads(first.read_text())
    spec = load_ts_optimizer_benchmark(first)

    assert len(spec.cases) == 12
    assert [case.difficulty for case in spec.cases].count("easy") == 4
    assert [case.difficulty for case in spec.cases].count("med") == 4
    assert [case.difficulty for case in spec.cases].count("hard") == 4
    assert spec.qc_model == {"software": "xtb", "xtb_model": "gfn2", "nprocs": 1}
    assert all(case.charge == 0 and case.multiplicity == 1 for case in spec.cases)
    assert all(Path(case.ts_guess).is_file() for case in spec.cases)
    assert all(case.reference_energy_hartree == pytest.approx(-1.2345) for case in spec.cases)
    assert raw["source"]["pilot_selection"]["reaction_count"] == 4
    assert raw["source"]["dataset_reaction_count"] == 8
    assert raw["source"]["dataset_starting_geometry_count"] == 24
    assert raw["source"]["license"].startswith("Copyright")
    assert raw["source"]["zenodo_archive_md5"] == adapter.ARCHIVE_MD5

    second = adapter.prepare_rgd1_tsopt_manifest(
        dataset, dataset / "pilot_second.json", reactions=4, seed=73,
    )
    second_raw = json.loads(second.read_text())
    assert [case["id"] for case in raw["cases"]] == [
        case["id"] for case in second_raw["cases"]
    ]
    selected_ids = raw["source"]["pilot_selection"]["selected_reactions"]
    dataset_rows = list(csv.DictReader((dataset / "manifest.tsv").open(), delimiter="\t"))
    selected_counts = {
        case_id: int(next(row["n_atoms"] for row in dataset_rows if row["id"] == case_id))
        for case_id in selected_ids
    }
    assert {next(low for low, high in adapter.ATOM_BINS if low <= n <= high)
            for n in selected_counts.values()} == {4, 9, 13, 17}


def test_rgd1_adapter_preserves_tier_metadata_and_refuses_overwrite(tmp_path):
    dataset = _dataset(tmp_path)
    path = adapter.prepare_rgd1_tsopt_manifest(
        dataset, dataset / "pilot.json", reactions=1, seed=11, tiers=("hard",),
        ts_fmax=0.02, ts_max_cycles=80,
    )
    spec = load_ts_optimizer_benchmark(path)
    assert len(spec.cases) == 1
    assert spec.cases[0].difficulty == "hard"
    assert spec.cases[0].perturbation == {
        "method": "GFN2-xTB mass-weighted Hessian normal modes",
        "tier": "hard",
        "alpha_react_A": 0.15,
        "beta_nonreact_A": 0.12,
        "total_rmsd_A": 0.14,
        "react_overlap": 0.447,
    }
    assert spec.settings["ts_fmax"] == 0.02
    assert spec.settings["ts_max_cycles"] == 80
    with pytest.raises(adapter.TSOptimizerBenchmarkError, match="output already exists"):
        adapter.prepare_rgd1_tsopt_manifest(dataset, path, reactions=1)


def test_rgd1_adapter_rejects_missing_tier_and_wrong_atom_order(tmp_path):
    dataset = _dataset(tmp_path)
    manifest_path = dataset / "manifest.tsv"
    rows = list(csv.DictReader(manifest_path.open(), delimiter="\t"))
    manifest_path.write_text("\t".join(rows[0]) + "\n" + "\t".join(
        rows[0][key] for key in rows[0]
    ) + "\n")
    with pytest.raises(adapter.TSOptimizerBenchmarkError, match="each of easy, med, and hard"):
        adapter.prepare_rgd1_tsopt_manifest(dataset, tmp_path / "bad.json", reactions=1)

    dataset = _dataset(tmp_path / "other")
    case_root = dataset / "cases" / "MR_sample_0"
    (case_root / "product.xyz").write_text("4\nwrong order\nH 0 0 0\nC 1 0 0\nH 2 0 0\nH 3 0 0\n")
    with pytest.raises(adapter.TSOptimizerBenchmarkError, match="different ordered atom list"):
        adapter.prepare_rgd1_tsopt_manifest(dataset, tmp_path / "bad_order.json", reactions=1)


def test_rgd1_prepare_cli_options():
    args = argument_parse([
        "prepare-rgd1", "--dataset-dir", "dataset", "--output", "pilot.json",
        "--reactions", "12", "--seed", "4", "--tiers", "easy", "hard",
        "--ts-fmax", "0.02", "--ts-max-cycles", "100",
    ])
    assert args.command == "prepare-rgd1"
    assert args.reactions == 12 and args.seed == 4
    assert args.tiers == ["easy", "hard"]
    assert args.ts_fmax == pytest.approx(0.02)
    assert args.ts_max_cycles == 100
