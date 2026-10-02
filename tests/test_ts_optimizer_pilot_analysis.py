from benchmarks.ts_optimizers.analysis.analyze_rgd1_pilot import _cluster_statistics


def test_pilot_analysis_resamples_and_compares_at_reaction_level():
    pairs = []
    for reaction_id, geometric_success, sella_success in (
        ("reaction_a", (True, True, False), (True, False, False)),
        ("reaction_b", (False, False, False), (False, False, True)),
    ):
        for tier, geometric, sella in zip(
            ("easy", "med", "hard"), geometric_success, sella_success,
        ):
            pairs.append((
                f"{reaction_id}_{tier}",
                {"case_id": f"{reaction_id}_{tier}", "outcome":
                 "reaction_connected_success" if geometric else "other"},
                {"case_id": f"{reaction_id}_{tier}", "outcome":
                 "reaction_connected_success" if sella else "other"},
            ))

    stats = _cluster_statistics(pairs, replicates=1000, seed=7)

    assert stats["reaction_count"] == 2
    assert stats["rates"] == {"geometric": 1 / 3, "sella": 1 / 3}
    assert stats["difference"] == 0
    assert stats["sign_flip_p"] == 1
