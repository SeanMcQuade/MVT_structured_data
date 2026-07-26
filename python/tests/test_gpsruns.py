"""Parity tests for control-vehicle run parsing.

``parse_gps_data`` decides how each vehicle's day is split into testbed runs.
The released GPS file holds those runs after further preprocessing (10 Hz
resampling and trimming, not yet ported), so the check here is structural: the
same number of runs, per vehicle, each released run falling inside the parsed
run it came from.
"""

from __future__ import annotations

import sys
from collections import Counter
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from mvtpy import avdist  # noqa: E402
from mvtpy.gpsruns import GpsRunOptions, day_time_limits, parse_gps_data  # noqa: E402

WORKSPACE = Path(__file__).resolve().parents[3]
CARS_GPS = WORKSPACE / "data" / "cars" / "cars_gps"
RELEASED_GPS = WORKSPACE / "results" / "gps" / "CIRCLES_GPS_10Hz_2022-11-16.json"

pytest.importorskip("pandas", reason="pandas is needed to read the vehicle CSVs")


@pytest.fixture(scope="module")
def parsed_runs():
    if not CARS_GPS.is_dir():
        pytest.skip("raw vehicle GPS data not available")
    return parse_gps_data(CARS_GPS, 16)


@pytest.fixture(scope="module")
def released_runs():
    if not RELEASED_GPS.is_file():
        pytest.skip("released GPS results not available")
    return avdist.load_av_runs(RELEASED_GPS)


def test_day_limits_are_local_time():
    lower, upper = day_time_limits(16)
    assert upper - lower == 15 * 3600          # 03:00 to 18:00
    assert lower == pytest.approx(1668589200)  # 2022-11-16 03:00 America/Chicago


def test_run_totals_match_released(parsed_runs, released_runs):
    assert len(parsed_runs) == len(released_runs)


def test_per_vehicle_run_counts_match_released(parsed_runs, released_runs):
    mine = Counter(run.vin for run in parsed_runs)
    theirs = Counter(int(run.av_id) for run in released_runs)
    assert mine == theirs


def test_every_released_run_falls_inside_a_parsed_run(parsed_runs, released_runs):
    """The released runs are trimmed versions of these, so they must be nested."""
    by_vehicle = {}
    for run in parsed_runs:
        by_vehicle.setdefault(run.vin, []).append(run)

    orphans = []
    for released in released_runs:
        candidates = by_vehicle.get(int(released.av_id), [])
        inside = any(run.starting_time - 0.5 <= released.first_timestamp
                     and released.last_timestamp <= run.ending_time + 0.5
                     for run in candidates)
        if not inside:
            orphans.append(released.av_id)
    assert not orphans, f"{len(orphans)} released runs not inside a parsed run"


def test_runs_are_directional_sweeps(parsed_runs):
    """Each run crosses the testbed, so direction agrees with net displacement."""
    for run in parsed_runs:
        assert run.direction in (1.0, -1.0)
        assert run.direction == (1.0 if run.x_position[-1] > run.x_position[0] else -1.0)


def test_runs_satisfy_the_filter_thresholds(parsed_runs):
    options = GpsRunOptions()
    for run in parsed_runs:
        assert run.ending_time - run.starting_time > options.min_run_time
        span = abs(run.x_position[-1] - run.x_position[0])
        assert span > options.min_run_length_km * options.km_to_ft


def test_run_numbers_increment_per_vehicle(parsed_runs):
    seen = {}
    for run in parsed_runs:
        expected = seen.get(run.vin, 0) + 1
        assert run.run_num == expected
        seen[run.vin] = expected


def _synthetic_vehicle_day(options: GpsRunOptions, exit_index: int = 800):
    """One westbound sweep that leaves the testbed at `exit_index`.

    Sized to clear the run filters (>60 s, >1.2 km). Sample `exit_index` is the
    first one outside the testbed; MATLAB keeps it as the run's final sample.
    """
    import numpy as np
    import pandas as pd

    n = exit_index + 100
    x = np.empty(n)
    x[:exit_index] = np.linspace(20000.0, options.min_rcs_x + 1.0, exit_index)
    x[exit_index:] = options.min_rcs_x - 100.0      # outside from here on
    y = np.linspace(100.0, 50.0, n)                 # moving toward the centre
    return pd.DataFrame({
        "Systime": 1668772800.0 + np.arange(n) * 0.1,
        "rcs_x": x, "rcs_y": y,
        "Long": np.full(n, -86.6), "Lat": np.full(n, 36.0),
        "state_x": x, "state_y": y,
        "can_speed": np.full(n, 25.0),
        "Status": ["A"] * n,
        "control_active": pd.array(["True"] * n, dtype="string"),
    })


def test_run_keeps_the_first_sample_outside_the_testbed():
    """MATLAB slices vehTable(1:runEnd,:) inclusively.

    Slicing to `end` instead dropped one sample from the tail of every run,
    which shortened 214 of the 772 released records for 2022-11-18.
    """
    from mvtpy.gpsruns import _split_runs

    options = GpsRunOptions()
    exit_index = 800
    table = _synthetic_vehicle_day(options, exit_index)
    runs = _split_runs(table, vehicle_id=1, options=options)

    assert runs, "synthetic sweep should produce a run"
    assert runs[0].timestamp.size == exit_index + 1
    # The kept sample is the one outside the testbed.
    assert runs[0].x_position[-1] < options.min_rcs_x
