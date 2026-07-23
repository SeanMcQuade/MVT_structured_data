"""Parity tests for lane identification and lane-change clipping.

These decide how a raw trajectory is split into released segments, so the
strongest available check is to reproduce the released segmentation exactly:
same segment identifiers, same sample windows, same corrected lateral position,
same lane number.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from mvtpy import lanes  # noqa: E402
from mvtpy.kinematics import FT_TO_METER  # noqa: E402
from mvtpy.matround import round_decimals  # noqa: E402


# ---------------------------------------------------------------------------
# unit behavior


def test_interp1_extrapolates_beyond_the_sample_range():
    xp = np.array([0.0, 1.0, 2.0])
    fp = np.array([0.0, 2.0, 4.0])

    inside = lanes.interp1_linear_extrap(xp, fp, np.array([0.5, 1.5]))
    assert list(inside) == [1.0, 3.0]

    # NumPy's interp would clamp these to 0 and 4; MATLAB continues the line.
    outside = lanes.interp1_linear_extrap(xp, fp, np.array([-1.0, 3.0]))
    assert list(outside) == [-2.0, 6.0]


def test_interp1_is_exact_at_nodes():
    xp = np.array([0.0, 1.0, 2.0])
    fp = np.array([3.0, -1.0, 7.0])
    assert list(lanes.interp1_linear_extrap(xp, fp, xp)) == list(fp)


def test_median_filter_edge_handling():
    """Ends repeat the first and last computed values, as in the MATLAB loop."""
    values = np.arange(20, dtype=float)
    filtered = lanes._median_filter(values, window=10)

    buffer = 5
    assert np.all(filtered[:buffer - 1] == filtered[buffer - 1])
    assert np.all(filtered[len(values) - buffer:] == filtered[len(values) - buffer - 1])


def test_median_filter_short_trajectory_is_constant():
    values = np.array([4.0, 1.0, 2.0])
    assert list(lanes._median_filter(values, window=10)) == [2.0, 2.0, 2.0]


# ---------------------------------------------------------------------------
# parity against the released segmentation


def test_driving_line_shape(driving_line):
    options = lanes.LaneIdentificationOptions()
    assert driving_line.shift.size == options.n_x_cells
    assert driving_line.x_cells.size == options.n_x_cells
    # The lateral wiggle of the roadway is a fraction of a lane width.
    assert np.abs(driving_line.shift).max() < options.lane_width


def test_every_released_segment_is_reproduced(produced_segments, released_segments):
    missing = sorted(set(released_segments) - set(produced_segments))
    assert not missing, f"{len(missing)} released segments not produced, e.g. {missing[:5]}"


def test_segment_windows_match_released(produced_segments, released_segments):
    for identifier, released in released_segments.items():
        produced = produced_segments[identifier]
        assert len(produced["timestamp"]) == len(released["timestamp"])
        assert np.array_equal(round_decimals(produced["timestamp"], 4),
                              np.asarray(released["timestamp"]))


def test_corrected_lateral_position_matches_released(produced_segments, released_segments):
    """y_position after clipping is the corrected y; released in meters."""
    for identifier, released in released_segments.items():
        produced = produced_segments[identifier]
        ours = round_decimals(np.asarray(produced["y_position"]) * FT_TO_METER, 4)
        assert np.array_equal(ours, np.asarray(released["y_position_corrected_meters"]))


def test_lane_numbers_match_released(produced_segments, released_segments):
    for identifier, released in released_segments.items():
        assert float(produced_segments[identifier]["lane"]) == float(released["lane_number"])


def test_segment_endpoints_match_released(produced_segments, released_segments):
    for identifier, released in released_segments.items():
        produced = produced_segments[identifier]
        assert round_decimals(produced["first_timestamp"], 4) == released["first_timestamp"]
        assert round_decimals(produced["last_timestamp"], 4) == released["last_timestamp"]


def test_segment_identifiers_are_suffixed_in_order(produced_segments):
    """Parts of one trajectory are numbered from 0 upward, as MATLAB does."""
    parts: dict[str, list[int]] = {}
    for identifier in produced_segments:
        oid, _, part = identifier.rpartition("-")
        parts.setdefault(oid, []).append(int(part))

    for oid, numbers in parts.items():
        numbers.sort()
        assert numbers == list(range(len(numbers))), f"gap in parts for {oid}: {numbers}"
