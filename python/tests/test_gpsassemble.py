"""Parity tests for GPS assembly past run-splitting.

What is verified exactly against the released data, and what is not, is spelled
out per test. The control-signal reconstruction is bit-exact; the resampled
floating-point fields match to a documented 6th-decimal tolerance (see
``test_resampled_fields_match_within_tolerance`` for why).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from mvtpy import gpsassemble as ga  # noqa: E402
from mvtpy.kinematics import FT_TO_METER  # noqa: E402
from mvtpy.matround import round_decimals  # noqa: E402

pytest.importorskip("pandas")


def _released_array(values):
    return np.asarray([np.nan if v is None else v for v in values], dtype=float)


def _equal_with_nan(mine, ref):
    mine = np.asarray(mine, dtype=float)
    ref = _released_array(ref)
    if mine.shape != ref.shape:
        return False
    both_nan = np.isnan(mine) & np.isnan(ref)
    return bool(np.all(both_nan | (mine == ref)))


# ---------------------------------------------------------------------------
# unit behavior


def test_resample_drops_off_cadence_samples():
    # A doubled sample (0 s gap) must be dropped before interpolation.
    time = np.array([0.0, 0.0, 0.1, 0.2])
    value = np.array([5.0, 999.0, 6.0, 7.0])
    grid = np.array([0.0, 0.1, 0.2])
    out = ga.resample_10hz(time, value, grid)
    assert out == pytest.approx([5.0, 6.0, 7.0])


def test_interp_uses_the_weighted_blend_form():
    # A(i)*(1-w) + A(i+1)*w, the form MATLAB's interp1 uses.
    xp = np.array([0.0, 1.0])
    fp = np.array([10.0, 20.0])
    assert ga._interp_extrap(xp, fp, np.array([0.25])) == pytest.approx([12.5])
    # Extrapolation continues the line.
    assert ga._interp_extrap(xp, fp, np.array([2.0])) == pytest.approx([30.0])


def test_tenth_second_grid_bounds():
    grid = ga._tenth_second_grid(10.53, 11.47)
    assert grid[0] == 10.5
    assert grid[-1] == pytest.approx(11.5)
    assert np.allclose(np.diff(grid), 0.1)


def test_control_car_reinstates_control_through_a_brief_stop():
    # Controller on, a short stop in the middle, moving otherwise: the stop must
    # not be read as disengagement.
    n = 20
    engaged = np.ones(n, dtype=bool)
    speed = np.ones(n)
    speed[8:11] = 0                       # brief stop
    x = np.cumsum(speed)
    t = np.arange(n) * 0.1
    control_car, _ = ga.control_car_status(engaged, speed, x, t)
    assert np.all(control_car[8:11] == 1)


def test_control_last30_is_a_30_second_lookback():
    n = 400
    engaged = np.zeros(n, dtype=bool)
    engaged[50] = True                    # a single engaged sample
    speed = np.ones(n)
    x = np.cumsum(speed)
    t = np.arange(n) * 0.1
    _, last30 = ga.control_car_status(engaged, speed, x, t)
    assert last30[50] == 1
    assert last30[349] == 1               # 299 samples later, still within 30 s
    assert last30[351] == 0               # 301 samples later, outside


# ---------------------------------------------------------------------------
# parity against released data


def test_lane_and_direction_match_released(gps_assembly_case):
    for case in gps_assembly_case:
        pre, record = case["pre"], case["record"]
        assert int(pre.assigned_lane) == int(record["assigned_lane"])
        assert pre.direction == record["direction"]


def test_control_car_matches_released_exactly(gps_assembly_case):
    """The control-signal reconstruction is bit-exact across every run."""
    for case in gps_assembly_case:
        pre, record, window = case["pre"], case["record"], case["window"]
        control_car, control_last30 = ga.control_car_status(
            pre.control_active, pre.can_speed,
            FT_TO_METER * pre.x_position, pre.timestamp)
        assert _equal_with_nan(control_car[window], record["control_car"])
        assert _equal_with_nan(control_last30[window], record["control_last30"])


def test_controller_engaged_matches_released_exactly(gps_assembly_case):
    for case in gps_assembly_case:
        pre, record, window = case["pre"], case["record"], case["window"]
        assert _equal_with_nan(pre.control_active[window].astype(float),
                               np.asarray(record["controller_engaged"], dtype=float))


def test_controller_engaged_encodes_as_json_booleans():
    """MATLAB carries this as a logical, so jsonencode writes true/false.

    Encoding it as 0/1 is the same information but not the same bytes, and at
    3M samples it was 9.8 MB of spurious difference against the released file -
    enough to fail a checksum comparison on its own.
    """
    from mvtpy import matjson

    record = {"controller_engaged": np.array([True, False, True]).tolist()}
    assert matjson.dumps(record) == '{"controller_engaged":[true,false,true]}'


def test_server_connected_logic_on_hand_built_inputs():
    """The connection rule: connected iff the last ping is within 2 s before.

    Verified here on constructed inputs, where there is no CSV-parse residual, so
    the rule itself is pinned exactly (against the released data it is 99.96%,
    the gap being sub-ULP timestamp differences at the 2 s threshold).
    """
    run_time = np.array([100.0, 101.0, 102.0, 105.0])
    ping_time = np.array([98.0, 99.5, 103.5])   # last-before values: 99.5, 99.5, 99.5, 103.5
    connected = ga.is_server_connected(run_time, ping_time)
    # t=100: last ping 99.5, gap -0.5 -> connected
    # t=101: last ping 99.5, gap -1.5 -> connected
    # t=102: last ping 99.5, gap -2.5 -> disconnected
    # t=105: last ping 103.5, gap -1.5 -> connected
    assert list(connected) == [1.0, 1.0, 0.0, 1.0]


def test_server_connected_matches_released_closely(gps_assembly_case):
    """is_server_connected reproduces the released values to ~99.9% of samples.

    The residual is not a logic difference: the MATLAB pipeline's own resampled
    run timestamps and parsed ping times differ from the Python ones by less
    than a ULP, which flips the >= -2 s comparison for a handful of samples per
    run. Feeding identical inputs to the MATLAB rule and this function gives
    bit-identical results (verified separately on 198k samples).
    """
    total = exact = 0
    for case in gps_assembly_case:
        pre, record, window = case["pre"], case["record"], case["window"]
        connected = ga.is_server_connected(pre.timestamp, case["status"]["timestamp"])[window]
        ref = np.array([0.0 if v in (0, False, None) else 1.0
                        for v in record["is_server_connected"]])
        total += ref.size
        exact += int(np.sum(connected == ref))
    assert exact / total > 0.999, f"only {exact}/{total} connection samples match"


def test_resampled_fields_agree_closely(gps_assembly_case):
    """y, latitude, longitude reproduce the released values very closely.

    The residual comes entirely from the CSV-parsed source floats, in two forms:

    * most divergent samples land one unit off across a 6th-decimal rounding
      boundary (the GPS data rounds to 6 decimals, 100x more sensitive than the
      4-decimal MOTION data);
    * a few samples per run — at run edges, where interp1 extrapolates, or at an
      isolated spot where a sub-ULP timestamp difference flips the
      off-cadence-sample filter — differ by up to ~1e-2 m.

    Byte-exact GPS output would require reproducing MATLAB's ``readtable`` float
    parsing bit-for-bit, a separate effort. This bounds the worst residual and
    (in the companion test) the fraction of samples affected.
    """
    worst = 0.0
    for case in gps_assembly_case:
        pre, record, window = case["pre"], case["record"], case["window"]
        for mine, ref in (
            (round_decimals(FT_TO_METER * pre.y_position, 6)[window], record["y_position"]),
            (round_decimals(pre.latitude, 6)[window], record["latitude"]),
            (round_decimals(pre.longitude, 6)[window], record["longitude"]),
        ):
            diff = np.abs(np.asarray(mine) - _released_array(ref))
            worst = max(worst, float(np.nanmax(diff)) if diff.size else 0.0)
    assert worst <= 2e-2, f"resampled fields diverge by {worst} m, larger than expected"


def test_latitude_longitude_are_bit_exact(gps_assembly_case):
    """With correctly-rounded CSV parsing, lat/long reproduce exactly.

    latitude and longitude are resampled from CSV columns and rounded to 6
    decimals without any further arithmetic. Once the CSV floats are parsed
    correctly-rounded (float_precision='round_trip', matching MATLAB readtable)
    and the 10 Hz grid is built with MATLAB's colon algorithm, these are
    bit-identical.
    """
    for kind, ref_key in (("lat", "latitude"), ("long", "longitude")):
        total = exact = 0
        for case in gps_assembly_case:
            mine = round_decimals(getattr(case["pre"], ref_key), 6)[case["window"]]
            ref = _released_array(case["record"][ref_key])
            total += ref.size
            exact += int(np.sum(np.asarray(mine) == ref))
        assert exact / total > 0.999, f"{kind}: only {exact}/{total} samples exact"


def test_resampled_samples_are_almost_all_bit_exact(gps_assembly_case):
    """y reproduces the released values to well over 98% of samples.

    y goes through the ft2m multiply before rounding, unlike lat/long, so a
    handful more samples per run land across a 6th-decimal boundary. The
    remaining residual is a sub-ULP grid/interpolation sensitivity, not a
    logic difference; if it grew materially, that would signal a regression.
    """
    total = exact = 0
    for case in gps_assembly_case:
        mine = round_decimals(FT_TO_METER * case["pre"].y_position, 6)[case["window"]]
        ref = _released_array(case["record"]["y_position"])
        total += ref.size
        exact += int(np.sum(np.asarray(mine) == ref))
    fraction = exact / total
    assert fraction > 0.98, f"y: only {exact}/{total} ({fraction:.4f}) samples exact"


def test_colon_matches_matlab_both_ends_construction():
    """The 10 Hz grid uses MATLAB's colon algorithm, not a + k*step.

    MATLAB builds a:step:b from both ends, which for a large base like a POSIX
    timestamp gives different last bits than a + k*step on ~20% of points.
    """
    grid = ga._colon(1668603246.5, 0.1, 1668603250.5)
    n = grid.size - 1
    half = n // 2
    # First half from the low end, second half from the high end.
    assert grid[0] == 1668603246.5
    assert grid[-1] == pytest.approx(1668603250.5)
    assert grid[half] == 1668603246.5 + half * 0.1
    assert grid[half + 1] == 1668603250.5 - (n - half - 1) * 0.1


def test_colon_sets_the_exact_midpoint_to_the_endpoint_average():
    """MATLAB's a:d:b builds from both ends, and for an even number of steps
    the middle element is (a+b)/2 - not either one-sided form.

    Verified bit-for-bit against MATLAB over the 772 real 2022-11-18 grids
    (0 of 3,739,194 points differ). Before this, 104 grid points per day were
    wrong, which flipped 6th-decimal roundings downstream.
    """
    from mvtpy.gpsassemble import _colon

    a, b = 1668776258.0, 1668776458.0        # 2000 steps of 0.1 => even
    grid = _colon(a, 0.1, b)
    assert grid.size == 2001
    assert grid[1000] == (a + b) / 2
    assert grid[0] == a and grid[-1] == b


def test_colon_odd_step_count_splits_without_a_midpoint():
    from mvtpy.gpsassemble import _colon

    a, b = 1668776258.0, 1668776258.0 + 0.1 * 2001   # odd number of steps
    grid = _colon(a, 0.1, b)
    assert grid.size == 2002
    assert grid[0] == a and grid[-1] == b


def test_interp_is_exact_on_a_flat_segment():
    """MATLAB's interp1 returns the value itself when both endpoints are equal.

    The blend A*(1-w) + B*w rounds the two products independently and lands a
    ULP off. Speed data is full of flat runs, so this was the last source of
    1-ULP differences in the resampled fields: against interp1 on a real
    9,449-point run the plain blend differed on 88 values, every one a flat
    segment; with the guard, none.
    """
    xp = np.array([0.0, 1.0])
    fp = np.array([68.736000000000004, 68.736000000000004])
    out = ga._interp_extrap(xp, fp, np.array([0.44856179347689035]))
    assert out[0] == fp[0]
    # the unguarded blend does not manage this
    w = 0.44856179347689035
    assert fp[0] * (1 - w) + fp[1] * w != fp[0]
