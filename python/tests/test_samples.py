"""Tests for per-sample collection near engaged control vehicles.

The integer-cast semantics are pinned with unit tests; full parity against the
released ``.mat`` (88.8M samples, all fields exact) is exercised by
``test_samples_match_released`` when the data and h5py are available.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from mvtpy import samples  # noqa: E402

WORKSPACE = Path(__file__).resolve().parents[3]
SLIM_DIR = WORKSPACE / "results" / "slim" / "2022-11-16"
REFERENCE = (WORKSPACE / "results" / "figures" / "2022-11-16"
             / "samples_for_distance_analysis_16.mat")


def test_matlab_int_rounds_half_away_and_saturates():
    # Round half away from zero, unlike numpy's truncation.
    assert samples.matlab_int([2.5, -2.5, 1.4, -1.6], np.int16).tolist() == [3, -3, 1, -2]
    # Saturate rather than wrap.
    assert samples.matlab_int([40000, -40000], np.int16).tolist() == [32767, -32768]
    assert samples.matlab_int([300, -5], np.uint8).tolist() == [255, 0]


def test_downstream_sample_precedes_upstream_at_a_point():
    """A point near both an up- and downstream AV emits downstream first."""
    record = {
        "direction": -1,
        "timestamp": [1000.0],
        "speed_meters_per_second": [10.0],
        "fuel_rate_grams_per_second": [0.5],
        "x_position_meters": [100.0],
        "coarse_vehicle_class": 0,
        "lane_number": 3,
        "distance_to_downstream_engaged_av_meters": [50.0],
        "distance_to_upstream_engaged_av_meters": [-40.0],
    }
    out = samples.collect_samples([record], 16)
    assert out["samples_dist"].tolist() == [50.0, -40.0]      # downstream, then upstream


def test_samples_outside_max_dist_are_dropped():
    record = {
        "direction": -1,
        "timestamp": [1.0, 2.0],
        "speed_meters_per_second": [10.0, 10.0],
        "fuel_rate_grams_per_second": [0.5, 0.5],
        "x_position_meters": [1.0, 2.0],
        "coarse_vehicle_class": 0,
        "lane_number": 3,
        "distance_to_downstream_engaged_av_meters": [500.0, 1500.0],   # second too far
        "distance_to_upstream_engaged_av_meters": [],
    }
    out = samples.collect_samples([record], 16)
    assert out["samples_dist"].tolist() == [500.0]


def test_fuel_consumption_formula():
    record = {
        "direction": -1,
        "timestamp": [1.0],
        "speed_meters_per_second": [4.0],
        "fuel_rate_grams_per_second": [2.0],
        "x_position_meters": [1.0],
        "coarse_vehicle_class": 0,
        "lane_number": 3,
        "distance_to_downstream_engaged_av_meters": [10.0],
        "distance_to_upstream_engaged_av_meters": [],
    }
    out = samples.collect_samples([record], 16)
    assert out["samples_fcons"][0] == pytest.approx(2.0 / (1e-6 + 4.0))


def test_time_is_seconds_after_six_am_local():
    # 06:00 on 2022-11-16 is epoch 1668600000; a sample 100 s later is t=100.
    record = {
        "direction": -1,
        "timestamp": [1668600100.0],
        "speed_meters_per_second": [10.0],
        "fuel_rate_grams_per_second": [0.5],
        "x_position_meters": [1.0],
        "coarse_vehicle_class": 0,
        "lane_number": 3,
        "distance_to_downstream_engaged_av_meters": [10.0],
        "distance_to_upstream_engaged_av_meters": [],
    }
    out = samples.collect_samples([record], 16)
    assert int(out["samples_t"][0]) == 100


@pytest.mark.skipif(not REFERENCE.is_file() or not SLIM_DIR.is_dir(),
                    reason="released slim data / samples reference not available")
def test_samples_match_released():
    """Every sample array matches the released .mat exactly (needs h5py)."""
    h5py = pytest.importorskip("h5py")
    mine = samples.collect_samples_from_dir(SLIM_DIR, 16)
    with h5py.File(REFERENCE, "r") as ref:
        for field in samples.SAMPLE_DTYPES:
            expected = ref[field][:].ravel()
            assert mine[field].size == expected.size, field
            assert np.array_equal(mine[field], expected), field
