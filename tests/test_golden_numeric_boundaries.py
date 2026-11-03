"""Analytical coverage of numerical services controlled by workflow goldens."""

from types import SimpleNamespace

import numpy as np
import pytest

from pyar.selection import deduplication


@pytest.mark.scientific
@pytest.mark.tier1
@pytest.mark.parametrize(
    "stretch, expected_threshold", [(0.0, 0.05), (0.08, 0.06), (0.2, 0.15), (0.4, 0.15)]
)
def test_adaptive_threshold_uses_real_aligned_dimer_rmsd_and_existing_bounds(
    monkeypatch, stretch, expected_threshold
):
    # Each endpoint moves by `stretch` in opposite directions: after optimal
    # alignment, the RMSD is exactly `stretch`. Only the fingerprint prefilter
    # is controlled; the actual permutation/Kabsch/SVD computation still runs.
    first = SimpleNamespace(
        atoms_list=["He", "He"], coordinates=np.array([[-1.0, 0, 0], [1.0, 0, 0]])
    )
    second = SimpleNamespace(
        atoms_list=["He", "He"],
        coordinates=np.array([[-1.0 - stretch, 0, 0], [1.0 + stretch, 0, 0]]),
    )
    monkeypatch.setattr(deduplication, "_structure_is_similar", lambda *_: True)

    value = deduplication._adaptive_duplicate_rmsd_threshold([first, second])

    # Independent analytical expectations for 0.75*RMSD, floor .05 and cap .15;
    # tolerance covers floating-point linear algebra, not chemical calibration.
    assert value == pytest.approx(expected_threshold, rel=0, abs=1e-12)


def test_adaptive_threshold_singleton_preserves_existing_default():
    assert deduplication._adaptive_duplicate_rmsd_threshold([object()]) == 0.10
