"""Parity tests for the full data set (eastbound + reference trajectories).

The `full` stage was ported after `slim` and exercises code paths `slim` never
touches - eastbound lane assignment above all. These tests pin the behaviour
that differs, so a future edit to the shared lane code cannot silently break
the direction `slim` does not use.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from mvtpy import build, full, lanes  # noqa: E402


# ---------------------------------------------------------------------------
# eastbound lane handling
#
# MATLAB's assign_lanes does two things for eastbound that are easy to miss,
# because slim filters eastbound out before it ever reaches them:
#
#   yEast = -yEast                      (before the outlier bounds)
#   y_corr = -(Se*(-y - dl(x)) + Ce)    (negate, correct, negate back)
#
# Folding those signs into the scale/offset constants gives a different answer,
# because the driving line is subtracted from the flipped y. Porting it that way
# produced corrected positions off by up to 0.54 m, which then changed where
# lane-change clipping cut the trajectory and shifted every later record.


def _straight_record(direction, y_value, n=200):
    return {
        "direction": direction,
        "timestamp": np.arange(n, dtype=float) / 25.0,
        "x_position": np.linspace(1000.0, 2000.0, n),
        "y_position": np.full(n, y_value, dtype=float),
    }


def test_eastbound_driving_line_flips_y_before_filtering():
    """Raw eastbound y is negative; unflipped it falls outside the bounds."""
    options = lanes.LaneIdentificationOptions()
    east = [_straight_record(1, -1.5 * options.lane_width) for _ in range(40)]

    # With the flip, these samples are in range and produce a driving line.
    line = lanes.estimate_driving_line(east, direction=1)
    assert np.isfinite(line.shift).all()

    # Westbound records with the same magnitude but positive y are what the
    # westbound estimator expects, and must be unaffected by the eastbound path.
    west = [_straight_record(-1, 1.5 * options.lane_width) for _ in range(40)]
    west_line = lanes.estimate_driving_line(west, direction=-1)
    assert np.isfinite(west_line.shift).all()


def test_eastbound_correction_is_not_a_sign_folded_westbound():
    """The eastbound expression is structurally different, not just re-signed."""
    options = lanes.LaneIdentificationOptions()
    line = lanes.DrivingLine(x_cells=np.array([0.0, 5000.0]),
                             shift=np.array([0.7, 0.7]))
    record = _straight_record(1, -20.0, n=50)

    y_corr, _ = lanes.assign_lanes(record, line)

    y = np.asarray(record["y_position"], dtype=float)
    expected = -(options.scale_east * (-y - line(np.asarray(record["x_position"]))
                                       ) + options.offset_east)
    assert np.allclose(y_corr, expected, rtol=0, atol=0)

    # The naive form - same expression as westbound with eastbound constants -
    # gives a materially different answer, which is the bug this guards.
    naive = options.scale_east * (y - line(np.asarray(record["x_position"]))) \
        + options.offset_east
    assert not np.allclose(y_corr, naive, rtol=0, atol=1e-6)


def test_westbound_correction_is_unchanged_by_the_eastbound_fix():
    """slim's parity depends on this expression being exactly as it was."""
    options = lanes.LaneIdentificationOptions()
    line = lanes.DrivingLine(x_cells=np.array([0.0, 5000.0]),
                             shift=np.array([0.3, 0.3]))
    record = _straight_record(-1, 20.0, n=50)

    y_corr, _ = lanes.assign_lanes(record, line)
    x = np.asarray(record["x_position"], dtype=float)
    y = np.asarray(record["y_position"], dtype=float)
    expected = options.scale_west * (y - line(x)) + options.offset_west
    assert np.allclose(y_corr, expected, rtol=0, atol=0)


# ---------------------------------------------------------------------------
# reference trajectory


def test_reference_trajectory_matches_the_real_drive_end_to_end():
    """The two-phase reference covers the same distance in the same time."""
    t = np.linspace(0.0, 60.0, 601)
    v = np.full(t.size, 25.0)
    x = 25.0 * t

    a1, a2, v_ref, a_ref, x_ref = full.reference_trajectory(x, v, t)

    assert x_ref[0] == pytest.approx(x[0], abs=1e-9)
    assert x_ref[-1] == pytest.approx(x[-1], abs=1e-6)
    assert v_ref[0] == pytest.approx(v[0], abs=1e-9)
    # A constant-speed drive needs no acceleration in either phase.
    assert a1 == pytest.approx(0.0, abs=1e-9)
    assert a2 == pytest.approx(0.0, abs=1e-9)
    assert a_ref.shape == t.shape


def test_reference_trajectory_switches_form_when_speed_would_go_negative():
    """A near-stop drive takes the three-phase branch, and never reverses.

    This is the `any(vRef < 0)` branch. It is a discontinuity: a trajectory an
    ULP either side of the test takes a different formula and every reference
    field changes together. The test pins that the branch exists and produces a
    non-negative profile, not that any particular trajectory takes it.
    """
    t = np.linspace(0.0, 100.0, 1001)
    # Starts fast, ends nearly stopped, covering little ground: the two-phase
    # form implies a negative speed partway through.
    x = np.concatenate([np.linspace(0.0, 50.0, 500), np.full(501, 50.0)])
    v = np.concatenate([np.full(500, 30.0), np.full(501, 0.05)])

    _, _, v_ref, _, _ = full.reference_trajectory(x, v, t)
    assert (v_ref >= 0).all()


# ---------------------------------------------------------------------------
# stage wiring


def test_full_is_a_stage_but_never_a_default_target():
    """`full` must be buildable on request and skipped by a plain build.

    Nothing downstream reads it and it is ~1.7x the size of slim, so it is
    opt-in in both implementations. plan() also has to keep an explicitly
    requested non-default target: filtering through DEFAULT_TARGETS dropped it
    silently, which looked exactly like "full refuses to build".
    """
    assert "full" in build.STAGES
    assert "full" not in build.DEFAULT_TARGETS
    assert "full.py" in build.STAGE_SOURCES["full"]


def test_full_field_order_matches_the_matlab_struct():
    """45 fields: init_data_struct's 41, then the four upstream fields."""
    assert len(full.FIELD_ORDER) == 45
    assert full.FIELD_ORDER[0] == "trajectory_id"
    assert full.FIELD_ORDER[5] == "direction"
    assert full.FIELD_ORDER[-4:] == (
        "upstream_engaged_av_id",
        "distance_to_upstream_engaged_av_meters",
        "upstream_av_id",
        "distance_to_upstream_av_meters",
    )
    # The four fuel totals that the deterministic quadrature governs.
    for name in ("total_fuel_consumed_grams",
                 "total_fuel_consumed_flat_road_grams",
                 "total_reference_fuel_consumed_grams",
                 "total_reference_fuel_consumed_flat_road_grams"):
        assert name in full.FIELD_ORDER
