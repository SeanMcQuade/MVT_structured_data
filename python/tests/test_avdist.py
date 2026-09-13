"""Parity tests for distance-to-control-vehicle computation.

Runs the port on segments derived from raw data and requires the four distance
fields and their vehicle ids to equal what MATLAB wrote into the released data,
including the empty-versus-null distinction.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from mvtpy import avdist  # noqa: E402
from mvtpy.matround import round_decimals  # noqa: E402

#: released field -> (key in the helper's output, sign applied by the caller)
DISTANCE_FIELDS = {
    "distance_to_upstream_av_meters": ("distanceToAvUS", -1.0),
    "distance_to_downstream_av_meters": ("distanceToAvDS", 1.0),
    "distance_to_upstream_engaged_av_meters": ("distanceToAvUSEng", -1.0),
    "distance_to_downstream_engaged_av_meters": ("distanceToAvDSEng", 1.0),
}

ID_FIELDS = {
    "upstream_av_id": "AvIdUS",
    "downstream_av_id": "AvIdDS",
    "upstream_engaged_av_id": "AvIdUSEng",
    "downstream_engaged_av_id": "AvIdDSEng",
}


def _released_array(value):
    """Released JSON gives [] for empty and null for missing samples."""
    if value is None:
        return None
    values = np.atleast_1d(np.asarray(
        [np.nan if item is None else item for item in np.atleast_1d(value)], dtype=float))
    return None if values.size == 0 else values


def _same(ours, theirs) -> bool:
    if ours is None or theirs is None:
        return ours is None and theirs is None
    if ours.shape != theirs.shape:
        return False
    both_nan = np.isnan(ours) & np.isnan(theirs)
    return bool(np.all(both_nan | (ours == theirs)))


def test_interp1_returns_nan_outside_the_sample_range():
    xp = np.array([10.0, 11.0, 12.0])
    fp = np.array([0.0, 10.0, 20.0])

    values = avdist.interp1_nan_outside(xp, fp, np.array([9.5, 10.5, 12.0, 12.5]))
    assert np.isnan(values[0])
    assert values[1] == 5.0
    assert values[2] == 20.0
    assert np.isnan(values[3])


def test_round_half_away_array_preserves_nan():
    values = avdist.round_half_away_array(np.array([0.5, -0.5, np.nan, 1.4]))
    assert values[0] == 1.0
    assert values[1] == -1.0
    assert np.isnan(values[2])
    assert values[3] == 1.0


def test_av_runs_load(av_runs):
    assert av_runs, "no control-vehicle runs loaded"
    run = av_runs[0]
    assert run.timestamp.size == run.x_position.size
    assert run.first_timestamp <= run.last_timestamp
    assert run.engaged.size == run.timestamp.size


def test_distances_match_released(produced_segments, released_segments, av_runs):
    compared = 0
    for identifier, released in released_segments.items():
        segment = produced_segments[identifier]
        computed = avdist.distance_to_avs(segment, av_runs)
        compared += 1

        for field, (key, sign) in DISTANCE_FIELDS.items():
            ours = computed[key]
            if ours is not None:
                ours = round_decimals(sign * ours, 4)
            assert _same(ours, _released_array(released[field])), \
                f"{identifier}: {field} differs"

    assert compared, "no segments compared"


def test_av_ids_match_released(produced_segments, released_segments, av_runs):
    for identifier, released in released_segments.items():
        computed = avdist.distance_to_avs(produced_segments[identifier], av_runs)
        for field, key in ID_FIELDS.items():
            assert _same(computed[key], _released_array(released[field])), \
                f"{identifier}: {field} differs"


def test_empty_when_no_control_vehicle_qualifies(av_runs):
    """A trajectory in a lane with no control vehicle gets empty fields."""
    segment = {
        "timestamp": np.array([1.0, 2.0, 3.0]),   # far outside the experiment
        "x_position": np.array([310000.0, 310010.0, 310020.0]),
        "lane": 99.0,
        "direction": -1.0,
    }
    computed = avdist.distance_to_avs(segment, av_runs)
    assert all(value is None for value in computed.values())
