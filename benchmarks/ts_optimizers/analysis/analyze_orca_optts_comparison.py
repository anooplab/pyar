"""Analyze paired PyAR and ORCA OptTS results from the RGD1 pilot."""

from __future__ import annotations

import argparse
import json
import random
import re
import statistics
from collections import Counter, defaultdict
from pathlib import Path


METHODS = ("orca_optts", "geometric", "sella")
SUCCESS = "reaction_connected_success"
DISPLAY = {"orca_optts": "ORCA OptTS", "geometric": "geomeTRIC", "sella": "Sella"}


def _results(directory, method):
    pattern = (
        "cases/*/orca_optts/benchmark_result.json" if method == "orca_optts"
        else f"cases/*/{method}/benchmark_result.json"
    )
    loaded = {}
    for path in Path(directory).glob(pattern):
        result = json.loads(path.read_text(encoding="utf-8"))
        if result.get("optimizer") != method or result["case_id"] in loaded:
            raise ValueError(f"duplicate case or incorrect optimizer in {path}")
        loaded[result["case_id"]] = result
    return loaded


def _summary(values):
    values = [float(value) for value in values if value is not None]
    if not values:
        return "n/a"
    return f"{statistics.mean(values):.3f} / {statistics.median(values):.3f}"


def _orca_used_redundant_internals(result):
    if result.get("settings", {}).get("orca_coordinate_system") == "redundant_internal":
        return True
    run_directory = result.get("run_directory")
    if not run_directory:
        return False
    input_path = Path(run_directory) / "orca_optts.inp"
    return input_path.is_file() and "CoordSys redundant" in input_path.read_text(encoding="utf-8")


def _cluster_comparison(
    results, method, *, reference="orca_optts", replicates=20000, seed=20261003,
):
    groups = defaultdict(lambda: defaultdict(list))
    for case_id in sorted(results[reference]):
        group = case_id.rsplit("_", 1)[0]
        for name in (reference, method):
            groups[group][name].append(results[name][case_id]["outcome"] == SUCCESS)
    if not groups or any(
        set(values) != {reference, method}
        or len(values[reference]) != len(values[method])
        for values in groups.values()
    ):
        raise ValueError("paired reaction groups are incomplete")
    rates = {
        name: {group: sum(values[name]) / len(values[name]) for group, values in groups.items()}
        for name in (reference, method)
    }
    reaction_ids = sorted(groups)
    observed = statistics.mean(
        rates[reference][group] - rates[method][group]
        for group in reaction_ids
    )
    rng = random.Random(seed)
    boot = []
    for _ in range(replicates):
        sampled = [rng.choice(reaction_ids) for _ in reaction_ids]
        boot.append(statistics.mean(
            rates[reference][group] - rates[method][group]
            for group in sampled
        ))
    boot.sort()
    interval = (boot[int(0.025 * replicates)], boot[int(0.975 * replicates)])
    differences = [
        rates[reference][group] - rates[method][group]
        for group in reaction_ids
    ]
    # Zero differences never change the permutation statistic. Enumerate small
    # sets exactly, but bound work for larger benchmark collections.
    nonzero = [difference for difference in differences if difference != 0]
    exact = len(nonzero) <= 16
    samples = 1 << len(nonzero) if exact else replicates
    extreme = 0
    for iteration in range(samples):
        mask = iteration if exact else rng.getrandbits(len(nonzero))
        shuffled = sum(
            difference * (-1 if mask & (1 << index) else 1)
            for index, difference in enumerate(nonzero)
        ) / len(differences)
        if abs(shuffled) >= abs(observed) - 1e-12:
            extreme += 1
    return {
        "reaction_count": len(reaction_ids),
        "observed_difference": observed,
        "interval": interval,
        "sign_flip_p": extreme / samples if exact else (extreme + 1) / (samples + 1),
        "sign_flip_method": "exact" if exact else "Monte Carlo (seeded)",
        "clusters_with_success": {
            name: sum(rate > 0 for rate in values.values())
            for name, values in rates.items()
        },
    }


def _provenance_gaps(results):
    """Reject recorded contradictions; never fill in missing historical evidence."""
    from pyar.benchmarks.ts_optimizer import VALIDATION_SETTINGS

    gaps = set()
    for case_id in results["orca_optts"]:
        rows = [cases[case_id] for cases in results.values()]
        fields = {
            "charge": [row.get("case_qc_settings", {}).get("charge") for row in rows],
            "multiplicity": [row.get("case_qc_settings", {}).get("multiplicity") for row in rows],
        }
        models = []
        for row in rows:
            qc = row.get("qc_model", {})
            model = qc.get("xtb_model")
            if model is None and qc.get("method") == "GFN2-xTB via ORCA xtb interface":
                model = "gfn2"
            models.append(model)
        fields["Hamiltonian"] = models
        for key in VALIDATION_SETTINGS:
            fields[key] = [row.get("settings", {}).get(key) for row in rows]
        for key in ("reference_ts", "reactant", "product"):
            fields[f"{key} hash"] = [row.get("input_hashes", {}).get(key) for row in rows]
        fields["xTB executable"] = [
            row.get("qc_model", {}).get("xtb_executable")
            or row.get("runtime", {}).get("backend", {}).get("path") for row in rows
        ]
        for label, values in fields.items():
            recorded = [value for value in values if value is not None]
            if len(set(recorded)) > 1:
                raise ValueError(f"{label} mismatch for {case_id}")
            if len(recorded) < len(values):
                gaps.add(label)
    return sorted(gaps)


def analyze(
    orca_directory, paired_directory, *, geometric_directory=None,
    sella_directory=None, sella_cartesian_directory=None,
):
    results = {"orca_optts": _results(orca_directory, "orca_optts")}
    results["geometric"] = _results(geometric_directory or paired_directory, "geometric")
    results["sella"] = _results(sella_directory or paired_directory, "sella")
    if sella_cartesian_directory is not None:
        results["sella_cartesian"] = _results(sella_cartesian_directory, "sella")
    case_ids = set(results["orca_optts"])
    if not case_ids or any(set(method_results) != case_ids for method_results in results.values()):
        raise ValueError("ORCA, geomeTRIC, and Sella must contain the same cases")
    for case_id in case_ids:
        hashes = {method_results[case_id]["input_sha256"] for method_results in results.values()}
        if len(hashes) != 1 or not all(isinstance(value, str) and value for value in hashes):
            raise ValueError(f"input geometry hash mismatch for {case_id}")
    difficulties = {
        case_id: results["orca_optts"][case_id]["difficulty"] for case_id in case_ids
    }
    if any(
        results[name][case_id]["difficulty"] != difficulties[case_id]
        for name in results for case_id in case_ids
    ):
        raise ValueError("difficulty labels do not match across optimizers")
    _provenance_gaps(results)
    return results


def render(results):
    gaps = _provenance_gaps(results)
    ids = sorted(results["orca_optts"])
    example = results["orca_optts"][ids[0]]
    runtime = example.get("runtime", {})
    runtime_backend = runtime.get("backend", {})
    orca_version = runtime.get("orca_version", "not recorded")
    raw_xtb_version = runtime_backend.get("version", "not recorded")
    match = re.search(r"xtb version\s+([^\s]+)", raw_xtb_version, re.I)
    xtb_version = match.group(1) if match else raw_xtb_version.replace("\n", " ")
    pyar_commit = runtime.get("git_commit", "not recorded")
    orca_results = results["orca_optts"].values()
    matched_force_threshold = all(
        result.get("settings", {}).get("reference_pyAR_ts_fmax_applied_to_orca") is True
        for result in orca_results
    )
    sella_internal = all(
        result.get("settings", {}).get("sella_internal_coordinates") is True
        for result in results["sella"].values()
    )
    orca_internal = all(
        _orca_used_redundant_internals(result)
        for result in results["orca_optts"].values()
    )
    force_targets = {
        result.get("settings", {}).get("reference_pyAR_ts_fmax_ev_per_angstrom")
        for result in results["orca_optts"].values()
    } | {
        result.get("settings", {}).get("ts_fmax")
        for method in ("geometric", "sella")
        for result in results[method].values()
    }
    cycle_limits = {
        result.get("settings", {}).get("orca_optts_max_cycles")
        for result in results["orca_optts"].values()
    } | {
        result.get("settings", {}).get("ts_max_cycles")
        for method in ("geometric", "sella")
        for result in results[method].values()
    }
    common_force_target = next(iter(force_targets)) if len(force_targets) == 1 else None
    common_cycle_limit = next(iter(cycle_limits)) if len(cycle_limits) == 1 else None
    if matched_force_threshold and orca_internal and sella_internal and common_force_target is not None and common_cycle_limit is not None:
        protocol_lines = [
            "- All three optimizers used internal coordinates: ORCA redundant internal coordinates, geomeTRIC delocalized internals, and Sella internal coordinates.",
            f"- All used a {common_cycle_limit}-cycle limit and a common maximum-gradient target of `{float(common_force_target):g} eV/angstrom`. ORCA's maximum and RMS gradient cutoffs were converted from that target; its native energy and displacement criteria remained enabled, so the full convergence definitions are not identical.",
            "- ORCA used `InHess XTB2` and Bofill updates. PyAR used its existing Hessian setup for geomeTRIC and Sella.",
            "- Wall-time comparisons describe the supplied runs; contemporaneous execution and equal machine load have not been verified.",
        ]
    else:
        protocol_lines = [
            "- A common gradient threshold and internal-coordinate protocol could not be verified from all supplied records. Consult each run's settings and optimizer output before attributing differences to the algorithm.",
        ]
    lines = [
        "# ORCA OptTS comparison on the RGD1-TSopt-GFN2 pilot",
        "",
        "## Integrity and protocol",
        "",
        f"- Paired starts: {len(ids)}; reaction clusters: {len({i.rsplit('_', 1)[0] for i in ids})}.",
        "- All three methods used the exact same starting geometry for every case (SHA-256 checked).",
        f"- The first ORCA record reports ORCA {orca_version} and xTB `{xtb_version}`.",
        f"- PyAR revision recorded by the optimizer runs: `{pyar_commit}`.",
        "- Provenance incomplete; unverified fields: " + ", ".join(gaps) + "." if gaps else "- Recorded Hamiltonian, spin, validation settings, endpoint/reference hashes and executable paths match across methods.",
        *protocol_lines,
        "- Reaction success requires independent frequency, IRC, and endpoint validation. Missing provenance prevents verification that these checks used an identical protocol.",
        "- ORCA does not expose provider calls in the same way as the PyAR energy-gradient adapter. ORCA optimizer timing-record counts are descriptive only and are not equated to PyAR provider evaluations.",
        "",
        "## Reaction-connected outcomes",
        "",
        "| Method | Validated starts | Rate | Reaction groups with at least one validated start |",
        "|---|---:|---:|---:|",
    ]
    for name in METHODS:
        successes = sum(results[name][case_id]["outcome"] == SUCCESS for case_id in ids)
        group_count = (
            sum(any(results[name][case_id]["outcome"] == SUCCESS for case_id in ids
                    if case_id.rsplit("_", 1)[0] == group)
                for group in {case_id.rsplit("_", 1)[0] for case_id in ids})
        )
        group_total = len({case_id.rsplit("_", 1)[0] for case_id in ids})
        lines.append(f"| {DISPLAY[name]} | {successes}/{len(ids)} | {100 * successes / len(ids):.1f}% | {group_count}/{group_total} |")
    lines.extend([
        "",
        "| Paired validated-start result | geomeTRIC | Sella |",
        "|---|---:|---:|",
    ])
    for method in METHODS[1:]:
        both = sum(
            results["orca_optts"][case_id]["outcome"] == SUCCESS
            and results[method][case_id]["outcome"] == SUCCESS for case_id in ids
        )
        orca_only = sum(
            results["orca_optts"][case_id]["outcome"] == SUCCESS
            and results[method][case_id]["outcome"] != SUCCESS for case_id in ids
        )
        pyar_only = sum(
            results["orca_optts"][case_id]["outcome"] != SUCCESS
            and results[method][case_id]["outcome"] == SUCCESS for case_id in ids
        )
        neither = len(ids) - both - orca_only - pyar_only
        lines.append(f"| Both / ORCA only / PyAR only / neither ({DISPLAY[method]}) | {both} / {orca_only} / {pyar_only} / {neither} | — |")
    lines.extend([
        "",
        "The paired reaction-cluster bootstrap compares each method's success fraction within each reaction, then resamples the reaction groups present in these results. A cluster-level sign-flip test is also reported.",
        "",
        "| ORCA minus PyAR | Mean paired difference | Reaction-cluster bootstrap 95% interval | Sign-flip p (method) |",
        "|---|---:|---:|---:|",
    ])
    for method in METHODS[1:]:
        stats = _cluster_comparison(results, method)
        low, high = stats["interval"]
        lines.append(
            f"| {DISPLAY[method]} | {stats['observed_difference']:+.1%} | "
            f"{low:+.1%} to {high:+.1%} | {stats['sign_flip_p']:.3f} ({stats['sign_flip_method']}) |"
        )
    direct = _cluster_comparison(results, "sella", reference="geometric")
    low, high = direct["interval"]
    lines.extend([
        "",
        "Direct paired optimizer comparison (geomeTRIC minus Sella):",
        f"mean reaction-cluster difference {direct['observed_difference']:+.1%}; bootstrap 95% interval {low:+.1%} to {high:+.1%}; {direct['sign_flip_method']} sign-flip p={direct['sign_flip_p']:.3f}.",
    ])
    lines.extend([
        "",
        "## Outcomes and cost",
        "",
        "| Method | Outcome counts | Median optimizer cycles | Median optimizer wall time (s) | Median ORCA gradient-log events / PyAR provider evaluations |",
        "|---|---|---:|---:|---:|",
    ])
    for name in METHODS:
        counts = Counter(results[name][case_id]["outcome"] for case_id in ids)
        cycles = _summary(results[name][case_id].get("ts_optimizer_steps") for case_id in ids).split(" / ")[-1]
        wall = _summary(results[name][case_id].get("ts_wall_seconds") for case_id in ids).split(" / ")[-1]
        if name == "orca_optts":
            events = _summary(results[name][case_id].get("ts_orca_gradient_events") for case_id in ids).split(" / ")[-1]
        else:
            events = _summary(results[name][case_id].get("ts_backend_energy_gradient_evaluations") for case_id in ids).split(" / ")[-1]
        lines.append(f"| {DISPLAY[name]} | `{dict(sorted(counts.items()))}` | {cycles} | {wall} | {events} |")
    if "sella_cartesian" in results:
        both = sum(
            results["sella_cartesian"][case_id]["outcome"] == SUCCESS
            and results["sella"][case_id]["outcome"] == SUCCESS
            for case_id in ids
        )
        cartesian_only = sum(
            results["sella_cartesian"][case_id]["outcome"] == SUCCESS
            and results["sella"][case_id]["outcome"] != SUCCESS
            for case_id in ids
        )
        internal_only = sum(
            results["sella_cartesian"][case_id]["outcome"] != SUCCESS
            and results["sella"][case_id]["outcome"] == SUCCESS
            for case_id in ids
        )
        neither = len(ids) - both - cartesian_only - internal_only
        changed = sum(
            results["sella_cartesian"][case_id]["outcome"]
            != results["sella"][case_id]["outcome"]
            for case_id in ids
        )
        lines.extend([
            "",
            "## Sella coordinate-system comparison",
            "",
            f"The same {len(ids)} starting geometries were also run with Sella in Cartesian coordinates. Validated successes were {both + cartesian_only}/{len(ids)} for Cartesian and {both + internal_only}/{len(ids)} for internal coordinates; both succeeded on {both}, Cartesian alone on {cartesian_only}, internal alone on {internal_only}, and neither on {neither}. The full outcome label changed for {changed}/{len(ids)} starts. This shows that coordinate choice can affect outcomes, but this pilot does not establish a universal setting.",
        ])
    success_counts = {
        name: sum(results[name][case_id]["outcome"] == SUCCESS for case_id in ids)
        for name in METHODS
    }
    full_agreement = sum(
        len({results[name][case_id]["outcome"] == SUCCESS for name in METHODS}) == 1
        for case_id in ids
    )
    group_success_counts = {
        name: sum(
            any(results[name][case_id]["outcome"] == SUCCESS for case_id in ids
                if case_id.rsplit("_", 1)[0] == group)
            for group in {case_id.rsplit("_", 1)[0] for case_id in ids}
        )
        for name in METHODS
    }
    orca_converged = sum(
        results["orca_optts"][case_id].get("ts_result", {}).get("optimizer_converged") is True
        for case_id in ids
    )
    success_summary = ", ".join(
        f"{DISPLAY[name]} {success_counts[name]}/{len(ids)}"
        for name in METHODS
    )
    group_summary = ", ".join(
        f"{DISPLAY[name]} {group_success_counts[name]}/{len({case_id.rsplit('_', 1)[0] for case_id in ids})}"
        for name in METHODS
    )
    lines.extend([
        "",
        "## Conclusion",
        "",
        f"Validated reaction-connected outcomes were: {success_summary}. The corresponding counts of reaction groups with at least one validated start were: {group_summary}. All three methods agreed on success or failure for {full_agreement}/{len(ids)} paired starts. The paired reaction-cluster estimates and their intervals above describe this pilot; they do not establish that the methods are equivalent or that one is generally better.",
        "",
        f"ORCA's optimizer converged for {orca_converged}/{len(ids)} starts, while full validation determined whether each saddle connected the intended endpoints. An optimizer convergence message alone is not sufficient; frequency, IRC, and endpoint checks remain necessary.",
        "",
        "Optimizer wall times are descriptive for these recorded runs. ORCA's gradient-log event count is not equivalent to PyAR provider evaluations, so the count columns must not be read as directly comparable electronic-structure call counts.",
        "",
        "This pilot alone does not establish a general default optimizer; broader validation is still needed.",
        "",
    ])
    if matched_force_threshold and orca_internal and sella_internal:
        geo_success = success_counts["geometric"]
        sella_success = success_counts["sella"]
        sella_only = sum(
            results["geometric"][case_id]["outcome"] != SUCCESS
            and results["sella"][case_id]["outcome"] == SUCCESS
            for case_id in ids
        )
        recovered_text = (
            "It added one validated start"
            if sella_only == 1 else f"It added {sella_only} validated starts"
        )
        lines.extend([
            "## Provisional default and recovery recommendation",
            "",
            f"Retain the current default pending broader validation. Internal-coordinate Sella validated {sella_success}/{len(ids)} starts versus {geo_success}/{len(ids)} for geomeTRIC. Reaction-group coverage is reported above.",
            f"Keep internal-coordinate Sella available as an explicit retry after geomeTRIC fails end-to-end validation. {recovered_text} among the {len(ids) - geo_success} geomeTRIC failures here. Retry from the original TS guess, then run frequency, IRC, and endpoint checks again. The gain is small, so automatic retries should wait for broader validation.",
            "The comparison does not isolate the effect of P-RFO from Hessian initialization, updates, or energy-gradient call overhead. ORCA timing counters are not directly comparable to PyAR provider calls.",
            "",
        ])
    return "\n".join(lines)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("orca_run", help="directory containing ORCA cases/*/orca_optts results")
    parser.add_argument("paired_run", help="directory containing paired geomeTRIC/Sella results")
    parser.add_argument("--output", help="write Markdown report to this file (default: stdout)")
    parser.add_argument("--geometric-run", help="separate run directory containing geomeTRIC results")
    parser.add_argument("--sella-run", help="separate run directory containing Sella results")
    parser.add_argument("--sella-cartesian-run", help="optional Cartesian-coordinate Sella run for coordinate comparison")
    args = parser.parse_args(argv)
    report = render(analyze(
        args.orca_run, args.paired_run,
        geometric_directory=args.geometric_run,
        sella_directory=args.sella_run,
        sella_cartesian_directory=args.sella_cartesian_run,
    ))
    if args.output:
        Path(args.output).write_text(report, encoding="utf-8")
    else:
        print(report)


if __name__ == "__main__":
    main()
