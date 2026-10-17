import importlib.util
from pathlib import Path


PILOT_PATH = (Path(__file__).parents[1] / "benchmarks" / "microsolvation"
              / "ms4_perturbation_pilot.py")
SPEC = importlib.util.spec_from_file_location("ms4_perturbation_pilot", PILOT_PATH)
PILOT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PILOT)


def test_matched_budget_keeps_only_whole_completed_optimizations():
    attempts = [
        {"attempt": 0, "status": "completed", "backend_calls_before": 0,
         "backend_calls_after": 12},
        {"attempt": 1, "status": "failed_optimization", "backend_calls_before": 12,
         "backend_calls_after": 20},
        {"attempt": 2, "status": "completed", "backend_calls_before": 20,
         "backend_calls_after": 31},
        {"attempt": 3, "status": "completed", "backend_calls_before": 31,
         "backend_calls_after": 45},
    ]

    kept = PILOT.matched_completed_prefix(attempts, call_budget=31)

    # The failed optimization consumes budget, and the attempt crossing the
    # common cap is excluded as a whole rather than counted as a partial result.
    assert [attempt["attempt"] for attempt in kept] == [0, 2]


def test_common_budget_ends_at_last_complete_result_in_every_arm():
    baseline = {"attempts": [
        {"status": "completed", "backend_calls_after": 12},
        {"status": "completed", "backend_calls_after": 24},
    ]}
    perturb = {"attempts": [
        {"status": "completed", "backend_calls_after": 17},
        {"status": "failed_optimization", "backend_calls_after": 31},
    ]}

    assert PILOT.common_completed_budget({"baseline": baseline, "perturb": perturb}, 30) == 17
