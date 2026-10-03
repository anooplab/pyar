"""Create a paired descriptive analysis of an RGD1 pilot run."""

from __future__ import annotations

import argparse
import json
import random
import statistics
from collections import Counter, defaultdict
from pathlib import Path


OPTIMIZERS = ("geometric", "sella")
SUCCESS = "reaction_connected_success"


def _load_results(run_dir):
    results = {}
    for path in Path(run_dir).glob("cases/*/*/benchmark_result.json"):
        result = json.loads(path.read_text(encoding="utf-8"))
        results[(result["case_id"], result["optimizer"])] = result
    return results


def _mean_median(values):
    values = [float(value) for value in values if value is not None]
    if not values:
        return "n/a"
    return f"{statistics.mean(values):.3f} / {statistics.median(values):.3f}"


def _cluster_statistics(pairs, replicates=20000, seed=20261002):
    reactions = defaultdict(lambda: {optimizer: [] for optimizer in OPTIMIZERS})
    for _, geometric, sella in pairs:
        reaction_id = geometric["case_id"].rsplit("_", 1)[0]
        reactions[reaction_id]["geometric"].append(geometric["outcome"] == SUCCESS)
        reactions[reaction_id]["sella"].append(sella["outcome"] == SUCCESS)
    if any(
        len(values["geometric"]) != len(values["sella"])
        or not values["geometric"]
        for values in reactions.values()
    ):
        raise ValueError("each reaction must have paired starts for both optimizers")
    reaction_ids = sorted(reactions)
    cluster_rates = {
        optimizer: {
            reaction: sum(reactions[reaction][optimizer]) / len(reactions[reaction][optimizer])
            for reaction in reaction_ids
        }
        for optimizer in OPTIMIZERS
    }
    observed = {
        optimizer: statistics.mean(cluster_rates[optimizer].values())
        for optimizer in OPTIMIZERS
    }
    observed_difference = observed["geometric"] - observed["sella"]

    rng = random.Random(seed)
    bootstrap = {"geometric": [], "sella": [], "difference": []}
    for _ in range(replicates):
        sampled = [rng.choice(reaction_ids) for _ in reaction_ids]
        geometric_rate = statistics.mean(cluster_rates["geometric"][r] for r in sampled)
        sella_rate = statistics.mean(cluster_rates["sella"][r] for r in sampled)
        bootstrap["geometric"].append(geometric_rate)
        bootstrap["sella"].append(sella_rate)
        bootstrap["difference"].append(geometric_rate - sella_rate)
    intervals = {
        key: tuple(sorted(values)[int(q * replicates)] for q in (0.025, 0.975))
        for key, values in bootstrap.items()
    }

    differences = [
        cluster_rates["geometric"][reaction] - cluster_rates["sella"][reaction]
        for reaction in reaction_ids
    ]
    extreme = 0
    total = 1 << len(differences)
    for mask in range(total):
        permuted = sum(
            difference * (-1 if mask & (1 << index) else 1)
            for index, difference in enumerate(differences)
        ) / len(differences)
        if abs(permuted) >= abs(observed_difference) - 1e-12:
            extreme += 1
    return {
        "reaction_count": len(reaction_ids),
        "rates": observed,
        "intervals": intervals,
        "difference": observed_difference,
        "difference_interval": intervals["difference"],
        "sign_flip_p": extreme / total,
    }


def analyze(run_dir):
    run_dir = Path(run_dir).resolve()
    summary = json.loads((run_dir / "summary.json").read_text(encoding="utf-8"))
    results = _load_results(run_dir)
    cases = sorted({case_id for case_id, _ in results})
    if len(results) != 2 * len(cases) or set(OPTIMIZERS) != {
        optimizer for _, optimizer in results
    }:
        raise ValueError("run does not contain complete geometric/Sella pairs")
    pairs = []
    for case_id in cases:
        geometric = results[(case_id, "geometric")]
        sella = results[(case_id, "sella")]
        if geometric["input_sha256"] != sella["input_sha256"]:
            raise ValueError(f"paired input hashes differ for {case_id}")
        pairs.append((case_id, geometric, sella))
    if len(pairs) != summary.get("paired_cases"):
        raise ValueError("collected summary and case result files disagree")
    commits = {
        result.get("runtime", {}).get("git_commit")
        for result in results.values()
    }
    if len(commits) != 1 or None in commits:
        raise ValueError("pilot results do not share one recorded PyAR commit")
    example = next(iter(results.values()))
    xtb_version = example.get("runtime", {}).get("backend", {}).get("version", "")
    xtb_version = next(
        (line.strip() for line in xtb_version.splitlines() if "xtb version" in line.lower()),
        "unrecorded",
    ).lstrip("* ")

    outcome_counts = {
        optimizer: Counter(r["outcome"] for _, g, s in pairs
                           for r in [g if optimizer == "geometric" else s])
        for optimizer in OPTIMIZERS
    }
    successes = {
        optimizer: outcome_counts[optimizer].get(SUCCESS, 0)
        for optimizer in OPTIMIZERS
    }
    both_success = sum(g["outcome"] == SUCCESS and s["outcome"] == SUCCESS for _, g, s in pairs)
    geometric_only = sum(g["outcome"] == SUCCESS and s["outcome"] != SUCCESS
                        for _, g, s in pairs)
    sella_only = sum(s["outcome"] == SUCCESS and g["outcome"] != SUCCESS
                     for _, g, s in pairs)
    both_fail = len(pairs) - both_success - geometric_only - sella_only
    cluster_stats = _cluster_statistics(pairs)

    tiers = defaultdict(lambda: defaultdict(Counter))
    for _, geometric, sella in pairs:
        tier = geometric["difficulty"]
        tiers[tier]["geometric"][geometric["outcome"]] += 1
        tiers[tier]["sella"][sella["outcome"]] += 1

    metric_map = {
        "TS provider evaluations": "ts_backend_energy_gradient_evaluations",
        "TS optimizer steps": "ts_optimizer_steps",
        "TS wall time (s)": "ts_wall_seconds",
        "Total provider evaluations": "total_backend_energy_gradient_evaluations",
        "Validation provider evaluations": "validation_backend_energy_gradient_evaluations",
        "Validation wall time (s)": "validation_wall_seconds",
    }
    lines = [
        "# RGD1-TSopt-GFN2 15-reaction pilot analysis",
        "",
        "## Run integrity",
        "",
        f"- Run output: `{run_dir.name}` (raw job artifacts are stored locally, not in the repository)",
        f"- Benchmark manifest SHA256: `{summary['manifest_sha256']}`",
        f"- PyAR commit: `{next(iter(commits))}`",
        f"- Runtime: Python {example['runtime'].get('python', 'unknown')}; "
        f"geomeTRIC {example['runtime'].get('geometric_version', 'unknown')}; "
        f"Sella {example['runtime'].get('sella_version', 'unknown')}; {xtb_version}",
        f"- Shared TS settings: `ts_fmax={example['settings']['ts_fmax']} eV/angstrom`, "
        f"`ts_max_cycles={example['settings']['ts_max_cycles']}`",
        f"- Completed jobs: {summary['completed_runs']} / {summary['expected_runs']}; "
        f"incomplete: {summary['incomplete_runs']}",
        f"- Paired inputs with matching hashes: {len(pairs)} / {len(pairs)}",
        "- Experimental units: 15 reactions × 3 correlated difficulty tiers (45 starts).",
        "",
        "## Reaction-connected success",
        "",
        "Success requires PyAR endpoint-connection validation, not only optimizer "
        "convergence or a first-order saddle.",
        "",
        "| Optimizer | Successes | Rate | Reaction-cluster bootstrap 95% interval |",
        "|---|---:|---:|---:|",
    ]
    for optimizer in OPTIMIZERS:
        n = len(pairs)
        successes_n = successes[optimizer]
        low, high = cluster_stats["intervals"][optimizer]
        lines.append(
            f"| {optimizer} | {successes_n}/{n} | {successes_n / n:.1%} | "
            f"{low:.1%}–{high:.1%} |"
        )
    lines += [
        "",
        "| Paired outcome | Cases |",
        "|---|---:|",
        f"| Both succeed | {both_success} |",
        f"| geomeTRIC only succeeds | {geometric_only} |",
        f"| Sella only succeeds | {sella_only} |",
        f"| Neither succeeds | {both_fail} |",
        "",
        f"Reaction-cluster paired success-rate difference (geomeTRIC − Sella): "
        f"{cluster_stats['difference']:+.1%}; 95% cluster-bootstrap interval "
        f"{cluster_stats['difference_interval'][0]:+.1%} to "
        f"{cluster_stats['difference_interval'][1]:+.1%}; exact reaction-level "
        f"sign-flip p = {cluster_stats['sign_flip_p']:.3f} "
        f"(n = {cluster_stats['reaction_count']} reaction clusters).",
        "",
        "| Tier | Optimizer | Reaction-connected successes | Other outcomes |",
        "|---|---|---:|---|",
    ]
    for tier in sorted(tiers):
        for optimizer in OPTIMIZERS:
            counts = tiers[tier][optimizer]
            other = ", ".join(f"{k}: {v}" for k, v in sorted(counts.items()) if k != SUCCESS)
            lines.append(
                f"| {tier} | {optimizer} | {counts.get(SUCCESS, 0)}/"
                f"{sum(counts.values())} | {other or 'none'} |"
            )

    lines += [
        "",
        "## Cost and reference agreement",
        "",
        "Cost values are descriptive across all 45 jobs (mean / median). Validation "
        "cost is shown separately because it can dominate TS optimizer work and "
        "depends on whether earlier scientific gates passed.",
        "",
        "| Measure | geomeTRIC | Sella |",
        "|---|---:|---:|",
    ]
    for label, key in metric_map.items():
        values = {
            optimizer: [r.get(key) for _, g, s in pairs
                        for r in [g if optimizer == "geometric" else s]]
            for optimizer in OPTIMIZERS
        }
        lines.append(
            f"| {label} (mean / median) | {_mean_median(values['geometric'])} | "
            f"{_mean_median(values['sella'])} |"
        )
    lines += [
        "",
        "Reference geometry and energy comparison is available for successful runs "
        "only. The dataset reference uses tblite 0.6.0; this pilot used the xTB executable.",
        "",
        "| Optimizer | Successful reference RMSD, median / max (Å) | "
        "Absolute reference energy difference, median / max (Eh) |",
        "|---|---:|---:|",
    ]
    for optimizer in OPTIMIZERS:
        selected = [r for _, g, s in pairs
                    for r in [g if optimizer == "geometric" else s]
                    if r["outcome"] == SUCCESS]
        rmsds = [r["reference_rmsd_angstrom"] for r in selected]
        energy_diffs = [abs(r["reference_energy_difference_hartree"]) for r in selected]
        lines.append(
            f"| {optimizer} | {statistics.median(rmsds):.6g} / {max(rmsds):.6g} | "
            f"{statistics.median(energy_diffs):.6g} / {max(energy_diffs):.6g} |"
        )
    lines += [
        "",
        "## Interpretation and limits",
        "",
        "- The pilot shows equal validated success (10/45 each). Only two pairs "
        "were discordant, one in each direction.",
        "- Sella used fewer TS evaluations and less measured TS wall time on average, "
        "while mean end-to-end provider evaluations were higher. These are estimates "
        "from one host and this selected subset.",
        "- The three starts per reaction are correlated. Chemical diversity is 15 "
        "reactions, not 45 independent reactions; per-tier rates are descriptive.",
        "- Success intervals resample whole reactions. The paired test flips "
        "optimizer labels at the reaction level.",
        "- These results do not justify selecting a default. Continue to the broader "
        "paired benchmark and native ORCA reference first.",
        "- GFN2-xTB executable results are not claimed numerically identical to "
        "the tblite reference implementation.",
        "",
    ]
    return "\n".join(lines)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir", help="collected paired benchmark run directory")
    parser.add_argument("--output", required=True, help="Markdown report path")
    args = parser.parse_args(argv)
    report = analyze(args.run_dir)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(report, encoding="utf-8")
    print(output)


if __name__ == "__main__":
    main()
