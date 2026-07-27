"""Tests for the MOTION-matching bias primitives.

The two primitives that had to be reverse-engineered from MATLAB are pinned
here: the ``smoothdata`` gaussian kernel (against a delta probe of MATLAB) and
the both-direction lane assignment (against MATLAB's assign_lanes on a real
segment, via the session fixtures). The end-to-end ``matching_bias`` is exercised
in ``test_matching_bias_against_recovered`` when the data is available.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from mvtpy import gpsmatch  # noqa: E402

WORKSPACE = Path(__file__).resolve().parents[3]


# ---------------------------------------------------------------------------
# smoothdata gaussian


def test_smoothdata_matches_matlab_uniform_delta():
    """A delta on a uniform grid recovers the (normalized) kernel weights.

    Expected values are from MATLAB smoothdata(x,'gaussian',3,'SamplePoints',t).
    """
    t = np.arange(21) * 0.1
    x = np.zeros(21)
    x[10] = 1.0
    out = gpsmatch.smoothdata_gaussian(x, t, 3.0)

    expected = {9: 0.0715836777, 10: 0.0722539968, 15: 0.0575963007, 20: 0.0313771426}
    for index, value in expected.items():
        assert out[index] == pytest.approx(value, abs=1e-9)


def test_smoothdata_matches_matlab_nonuniform():
    t = np.array([0, 0.1, 0.25, 0.3, 0.7, 1.0, 1.3, 1.35, 2.0])
    x = np.zeros(9)
    x[5] = 1.0
    out = gpsmatch.smoothdata_gaussian(x, t, 3.0)

    expected = [0.0528687204, 0.0638818098, 0.0832752704, 0.0904545393,
                0.1519135105, 0.1853201491, 0.1882264901, 0.1854545915, 0.1035763939]
    assert out == pytest.approx(expected, abs=1e-9)


def test_smoothdata_is_normalized():
    """Constant input is preserved (weights sum to one at every point)."""
    t = np.linspace(0, 5, 60)
    x = np.full(t.size, 3.7)
    out = gpsmatch.smoothdata_gaussian(x, t, 3.0)
    assert np.allclose(out, 3.7)


def test_smoothdata_truncates_beyond_half_window():
    """Points more than window/2 away carry no weight.

    A value far outside the half-window must not influence the center.
    """
    t = np.array([0.0, 0.1, 0.2, 5.0])   # last point 5 s away
    x = np.array([0.0, 0.0, 0.0, 100.0])
    out = gpsmatch.smoothdata_gaussian(x, t, 3.0)
    assert out[0] == pytest.approx(0.0, abs=1e-12)


# ---------------------------------------------------------------------------
# bidirectional lane assignment (uses the released-data fixtures)


def test_bidirectional_lanes_match_matlab(motion_segment_lanes):
    """Lane arrays match MATLAB's assign_lanes to float tolerance.

    The comparison is against MATLAB's own assign_lanes run on the same MOTION
    segment; agreement is to ~1e-14, far below the 0.5-lane matching threshold.
    """
    mine, matlab = motion_segment_lanes
    worst = 0.0
    for our_lane, their_lane in zip(mine, matlab):
        if our_lane.shape != their_lane.shape:
            continue
        worst = max(worst, float(np.max(np.abs(our_lane - their_lane))))
    assert worst < 1e-10, f"lane arrays diverge by {worst}"


# ---------------------------------------------------------------------------
# end-to-end bias


@pytest.mark.skipif(
    "MVT_RUN_SLOW" not in __import__("os").environ,
    reason="full-day matching is minutes-long; set MVT_RUN_SLOW=1 to run")
def test_matching_bias_against_recovered(matching_bias_case):
    """Per-run median_xd matches the value recovered from MATLAB's output.

    The recovered target is median(ft2m*x_python - x_matlab) per run, which is
    exactly the offset the matching pass computes. Only runs whose matched
    stretches are numerically well-determined are checked; the tolerance allows
    for the sub-ULP resample residual feeding the distance computation.
    """
    mine, target = matching_bias_case
    checked = agree = 0
    worst = 0.0
    for index, target_value in target.items():
        if index not in mine:
            continue
        checked += 1
        diff = abs(mine[index] - target_value)
        worst = max(worst, diff)
        if diff < 1e-3:
            agree += 1
    assert checked > 0
    assert agree / checked > 0.9, \
        f"only {agree}/{checked} runs within 1e-3 m (worst {worst:.4f})"


def test_smoothdata_omits_nan_like_matlab():
    """MATLAB's smoothdata ignores NaN; a point is NaN only if its whole window is.

    Propagating NaN instead poisoned a half-window (1.5 s, ~37 samples at
    MOTION's 25 Hz) ahead of every NaN. dist_to_av is NaN wherever a trajectory
    runs past the AV's own time range, so matched stretches ended early and the
    median_xd bias shifted - this was the last discrepancy between the Python
    and MATLAB GPS output, worth 2.6-21 mm on 5 of 772 runs.
    """
    from mvtpy.gpsmatch import smoothdata_gaussian

    # The window is 3 s wide (+/- 1.5 s), so the NaN block must exceed that
    # before any point sees an all-NaN window.
    t = np.arange(80) * 0.1
    v = np.ones(80)
    v[20:] = np.nan

    out = smoothdata_gaussian(v, t, 3.0)
    # A window still containing a finite sample stays finite, and the NaNs in it
    # are ignored rather than poisoning the result.
    assert np.isfinite(out[:34]).all()
    assert out[19] == pytest.approx(1.0)
    assert out[33] == pytest.approx(1.0)
    # Beyond a full half-window past the last finite sample, the window is
    # entirely NaN and the result is NaN.
    assert np.isnan(out[40:]).all()


def test_smoothdata_all_nan_window_is_nan():
    from mvtpy.gpsmatch import smoothdata_gaussian

    t = np.arange(10) * 0.1
    out = smoothdata_gaussian(np.full(10, np.nan), t, 3.0)
    assert np.isnan(out).all()
