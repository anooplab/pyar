"""Tests for conservative accepted-geometry release evidence."""

import numpy as np
import pytest

from pyar.bonding_analysis import BondingAnalysis
from pyar.release import BondOrderTrajectory, ReleaseTracker


SYMBOLS = ["C", "C"]
FRAGMENTS = [[0], [1]]


def observe(tracker, distance, alpha_critical=0.01):
    return tracker.observe_accepted(
        SYMBOLS,
        np.array([[0.0, 0.0, 0.0], [distance, 0.0, 0.0]]),
        FRAGMENTS,
        alpha=0.1,
        alpha_critical=alpha_critical,
        alpha_target=0.1,
        segment_index=tracker.accepted_step,
    )


def test_transient_contact_does_not_release():
    tracker = ReleaseTracker()
    observe(tracker, 2.0)
    evidence = observe(tracker, 3.0)

    assert evidence.state == "DRIVING"
    assert evidence.persistence_counter == 0


def test_persistent_stabilized_low_resistance_contact_becomes_candidate():
    tracker = ReleaseTracker()
    observe(tracker, 1.4)
    observe(tracker, 1.38)
    evidence = observe(tracker, 1.37)

    assert evidence.state == "CANDIDATE"
    assert evidence.persistence_counter == 3
    assert evidence.forming_pairs == ((0, 1),)
    assert evidence.normalized_distances[0] < 1.0


def test_hysteresis_prevents_topology_flicker():
    tracker = ReleaseTracker()
    observe(tracker, 1.40)
    observe(tracker, 2.7)
    evidence = observe(tracker, 1.40)

    assert evidence.state == "DRIVING"
    assert evidence.persistence_counter == 1


def test_persistent_contact_triggers_probe_even_at_high_resistance():
    tracker = ReleaseTracker()
    observe(tracker, 1.4, alpha_critical=0.2)
    observe(tracker, 1.38, alpha_critical=0.2)
    evidence = observe(tracker, 1.37, alpha_critical=0.2)
    assert evidence.state == "CANDIDATE"
    assert evidence.reason == "persistent_stabilized_contact_probe_high_resistance"
    assert evidence.alpha_critical == 0.2
    assert tracker.low_resistance_count == 0


@pytest.mark.parametrize("radius_ratio", [1.0, 1.05, 1.10])
def test_contact_must_be_strictly_inside_covalent_radius_sum(radius_ratio):
    from pyar.data import new_atomic_data

    radii_sum = (
        new_atomic_data.covalent_radius["C"]
        + new_atomic_data.covalent_radius["N"]
    )
    tracker = ReleaseTracker()
    for _ in range(4):
        evidence = tracker.observe_accepted(
            ["C", "N"],
            np.array([[0.0, 0.0, 0.0], [radii_sum * radius_ratio, 0.0, 0.0]]),
            [[0], [1]],
            alpha=0.1,
            alpha_critical=0.01,
            alpha_target=0.1,
            segment_index=tracker.accepted_step,
        )
    assert evidence.state == "DRIVING"
    assert evidence.forming_pairs == ()


def test_marginal_contact_inside_radius_sum_is_not_release_ready():
    from pyar.data import new_atomic_data

    radii_sum = new_atomic_data.covalent_radius["C"] * 2
    tracker = ReleaseTracker()
    for _ in range(5):
        evidence = tracker.observe_accepted(
            SYMBOLS,
            np.array([[0.0, 0.0, 0.0], [0.96 * radii_sum, 0.0, 0.0]]),
            FRAGMENTS,
            alpha_critical=0.2,
        )
    assert evidence.state == "DRIVING"
    assert evidence.forming_pairs == ()


def test_release_distance_fraction_is_validated_and_restart_compatible():
    tracker = ReleaseTracker(distance_fraction=0.94)
    restored = ReleaseTracker(distance_fraction=0.94)
    restored.load_state_dict(tracker.state_dict())
    assert restored.distance_fraction == 0.94
    with pytest.raises(ValueError, match="release tolerances"):
        ReleaseTracker(distance_fraction=1.01)


def test_restart_preserves_release_evidence():
    tracker = ReleaseTracker()
    observe(tracker, 1.4)
    observe(tracker, 1.38)
    restored = ReleaseTracker()
    restored.load_state_dict(tracker.state_dict())

    evidence = observe(restored, 1.37)
    assert evidence.state == "CANDIDATE"
    assert evidence.accepted_step == 3


def bond_analysis(value, scheme="wiberg_ao"):
    return BondingAnalysis(
        scheme=scheme,
        bond_orders={(0, 1): value},
        metadata={"backend": "xtb"},
    )


def test_bond_order_growth_and_plateau_can_support_release():
    tracker = ReleaseTracker()
    observe(tracker, 1.4)
    tracker.observe_accepted(SYMBOLS, np.array([[0, 0, 0], [1.38, 0, 0]]), FRAGMENTS,
                             alpha_critical=0.01, bonding_analysis=bond_analysis(0.20))
    evidence = tracker.observe_accepted(
        SYMBOLS, np.array([[0, 0, 0], [1.37, 0, 0]]), FRAGMENTS,
        alpha_critical=0.01, bonding_analysis=bond_analysis(0.22),
    )

    assert evidence.state == "CANDIDATE"
    assert evidence.bond_order_available is True
    assert evidence.bond_order_scheme == "wiberg_ao"


def test_transient_bond_order_growth_does_not_support_release():
    tracker = ReleaseTracker()
    for distance, order in ((1.4, 0.10), (1.38, 0.30), (1.37, 0.10)):
        evidence = tracker.observe_accepted(
            SYMBOLS,
            np.array([[0, 0, 0], [distance, 0, 0]]),
            FRAGMENTS,
            alpha_critical=0.01,
            bonding_analysis=bond_analysis(order),
        )
    assert evidence.state == "DRIVING"


def test_negligible_bond_order_blocks_an_electronic_release():
    tracker = ReleaseTracker()
    for distance in (1.4, 1.38, 1.37):
        evidence = tracker.observe_accepted(
            SYMBOLS,
            np.array([[0, 0, 0], [distance, 0, 0]]),
            FRAGMENTS,
            alpha_critical=0.01,
            bonding_analysis=bond_analysis(0.0),
        )
    assert evidence.state == "DRIVING"


def test_unavailable_bond_order_uses_geometry_fallback():
    tracker = ReleaseTracker()
    unavailable = BondingAnalysis(None, available=False, reason="not_supported")
    for distance in (1.4, 1.38, 1.37):
        evidence = tracker.observe_accepted(
            SYMBOLS,
            np.array([[0, 0, 0], [distance, 0, 0]]),
            FRAGMENTS,
            alpha_critical=0.01,
            bonding_analysis=unavailable,
        )
    assert evidence.state == "CANDIDATE"
    assert evidence.bond_order_available is False


def test_omitted_early_bond_order_allows_geometry_probe_after_restart():
    tracker = ReleaseTracker()
    tracker.observe_accepted(
        SYMBOLS, np.array([[0, 0, 0], [2.5, 0, 0]]), FRAGMENTS,
        bonding_analysis=BondingAnalysis("mayer", {}),
    )
    restored = ReleaseTracker()
    restored.load_state_dict(tracker.state_dict())
    for distance, order in ((1.40, 0.5), (1.38, 0.8), (1.37, 0.8)):
        evidence = restored.observe_accepted(
            SYMBOLS, np.array([[0, 0, 0], [distance, 0, 0]]), FRAGMENTS,
            bonding_analysis=bond_analysis(order, "mayer"), alpha_critical=0.01,
        )
    assert evidence.state == "CANDIDATE"
    assert evidence.bond_order_delta["(0, 1)"] is None
    assert "incomplete_bond_order_geometry_fallback" in evidence.reason
    assert not restored.bond_order_tracker.last_diagnostics["qualifies"]


def test_missing_pair_does_not_override_complete_negative_bond_evidence():
    trajectory = BondOrderTrajectory()
    pairs = {(0, 2), (1, 3)}
    for _ in range(3):
        diagnostics = trajectory.observe(
            ["C"] * 4, pairs, pairs, BondingAnalysis("mayer", {(0, 2): 0.0}),
            [[0, 1], [2, 3]],
        )
    assert diagnostics["geometry_fallback_pairs"] == ["(1, 3)"]
    assert not diagnostics["permits_geometry_probe"]


def test_bond_order_scheme_change_after_unavailable_step_resets_history():
    trajectory = BondOrderTrajectory()
    pairs = {(0, 1)}
    trajectory.observe(SYMBOLS, pairs, pairs, bond_analysis(0.1, "mayer"))
    trajectory.observe(SYMBOLS, pairs, pairs, None)
    diagnostics = trajectory.observe(SYMBOLS, pairs, pairs, bond_analysis(0.8, "wiberg_ao"))
    assert diagnostics["delta_from_initial"]["(0, 1)"] is None
    assert trajectory.history[(0, 1)] == [0.8]


def test_comparable_bond_orders_are_recorded_without_a_universal_cutoff():
    trajectory = BondOrderTrajectory()
    diagnostics = trajectory.observe(
        ["C", "C", "C", "C"],
        {(0, 1), (2, 3)},
        {(2, 3)},
        BondingAnalysis(
            scheme="mayer",
            bond_orders={(0, 1): 1.2, (2, 3): 0.4},
        ),
    )

    assert diagnostics["comparable_bond_orders"]["(2, 3)"] == [1.2]


def test_two_forming_bonds_keep_separate_trajectories():
    trajectory = BondOrderTrajectory()
    pairs = {(0, 2), (1, 3)}
    first = BondingAnalysis("mayer", {(0, 2): 0.20, (1, 3): 0.10})
    second = BondingAnalysis("mayer", {(0, 2): 0.22, (1, 3): 0.12})
    trajectory.observe(["C"] * 4, pairs, pairs, first)
    diagnostics = trajectory.observe(["C"] * 4, pairs, pairs, second)

    assert diagnostics["delta_from_initial"]["(0, 2)"] == pytest.approx(0.02)
    assert diagnostics["delta_from_initial"]["(1, 3)"] == pytest.approx(0.02)
    assert set(trajectory.history) == pairs
