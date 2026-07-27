"""Lane identification and lane-change clipping.

Port of the ``assign_lanes`` and ``clip_lane_changes`` local functions in
``Scripts/generate_data_mvt_slim.m``. Together they decide how each raw I-24
MOTION trajectory becomes one or more released segments (the ``-0``, ``-1``,
... suffixes on ``trajectory_id``), and they produce the corrected lateral
position that the released data reports as ``y_position_corrected_meters``.

The algorithm has two stages:

1. **Driving line** (a whole-file statistic). Westbound trajectories longer
   than five seconds are sampled at 1 Hz, outliers outside the lateral bounds
   are dropped, and the roadway is divided into ``Nr_XCells`` cells along x. In
   each cell with enough samples, the lateral offsets are mapped onto a unit
   circle with period ``LaneWidth`` and averaged circularly, giving the local
   lateral shift of the roadway ("wiggle") in that cell.

2. **Per-trajectory assignment**. Each trajectory's lateral position is
   corrected by the interpolated driving line, converted to a fractional lane
   index, clamped to [0, 5], and median filtered over a 10-sample window.

Clipping then walks each trajectory, emitting a segment for every stretch that
stays within one lane, and discarding the parts spent changing lanes.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, List, Sequence

import numpy as np

from .matround import round_half_away

__all__ = [
    "LaneIdentificationOptions",
    "LaneChangeClippingOptions",
    "DrivingLine",
    "estimate_driving_line",
    "assign_lanes",
    "clip_lane_changes",
    "interp1_linear_extrap",
]


@dataclass(frozen=True)
class LaneIdentificationOptions:
    """Defaults copied from ``processingOpts.laneIdentificationOpts``."""

    lane_width: float = 12.0        # [ft]
    n_x_cells: int = 200
    min_cell_samples: int = 20
    y_up_lim: float = 5.0           # multiples of lane width
    y_low_lim: float = 0.5
    scale_west: float = 0.98        # Sw
    offset_west: float = 1.0        # Cw
    scale_east: float = 0.97        # Se
    offset_east: float = 1.0        # Ce
    median_window: int = 10         # ww


@dataclass(frozen=True)
class LaneChangeClippingOptions:
    """Defaults copied from ``processingOpts.laneChangeClippingOpts``."""

    lane_change_thresh: float = 0.5        # multiples of lane width
    max_lane_change_rate: float = 0.1      # lane widths per second
    min_clip_time: float = 0.5             # [s]
    change_buffer_thresh: float = 0.2      # multiples of lane width


def interp1_linear_extrap(xp: np.ndarray, fp: np.ndarray, x: np.ndarray) -> np.ndarray:
    """``interp1(xp, fp, x, 'linear', 'extrap')``.

    NumPy's ``interp`` clamps outside the sample range; MATLAB continues the
    nearest segment's line, which matters because trajectories reach beyond the
    x-range covered by the driving-line cells.

    Uses the same evaluation MATLAB's interp1 does, checked against it on real
    inputs: the weighted blend ``A*(1-w) + B*w``, except on a flat segment
    (``A == B``) where MATLAB returns the value itself. The algebraically
    equivalent slope form ``A + slope*(x-A)`` that this used before disagrees
    with MATLAB on ~25% of points in the last bit.
    """
    xp = np.asarray(xp, dtype=float)
    fp = np.asarray(fp, dtype=float)
    x = np.asarray(x, dtype=float)

    index = np.searchsorted(xp, x, side="right") - 1
    index = np.clip(index, 0, len(xp) - 2)
    left, right = fp[index], fp[index + 1]
    weight = (x - xp[index]) / (xp[index + 1] - xp[index])
    return np.where(left == right, left, left * (1 - weight) + right * weight)


@dataclass(frozen=True)
class DrivingLine:
    """Per-cell lateral shift of the roadway, with the cell edges it applies to."""

    x_cells: np.ndarray
    shift: np.ndarray

    def __call__(self, x: np.ndarray) -> np.ndarray:
        return interp1_linear_extrap(self.x_cells, self.shift, x)


def estimate_driving_line(
    trajectories: Iterable[dict],
    options: LaneIdentificationOptions = LaneIdentificationOptions(),
    direction: int = -1,
) -> DrivingLine:
    """Stage 1: circular mean of lateral offsets per x cell.

    ``trajectories`` may be any iterable of raw records, including the streaming
    reader, so this never holds the whole file in memory. Only ``x_position``,
    ``y_position``, ``timestamp``, and ``direction`` are used.
    """
    sampled_x: List[np.ndarray] = []
    sampled_y: List[np.ndarray] = []

    for record in trajectories:
        if record["direction"] != direction:
            continue
        time = np.asarray(record["timestamp"], dtype=float)
        if time.size < 2 or (time[-1] - time[0]) <= 5:
            continue
        # ss = ceil(1/ts): sub-sample to about 1 Hz
        step = int(np.ceil(1.0 / (time[1] - time[0])))
        sampled_x.append(np.asarray(record["x_position"], dtype=float)[::step])
        sampled_y.append(np.asarray(record["y_position"], dtype=float)[::step])

    if not sampled_x:
        raise ValueError("no trajectories long enough to estimate the driving line")

    x = np.concatenate(sampled_x)
    y = np.concatenate(sampled_y)

    upper = options.lane_width * options.y_up_lim
    lower = options.lane_width * options.y_low_lim
    keep = (y <= upper) & (y >= lower)
    x, y = x[keep], y[keep]

    edges = np.linspace(x.min(), x.max(), options.n_x_cells + 1)
    shift = np.zeros(options.n_x_cells)

    # Bin membership matches MATLAB's `edges(i) < x & x <= edges(i+1)`: cells are
    # open on the left, so the single sample sitting exactly on the lower edge
    # belongs to no cell (it maps to -1 and is ignored below), as in MATLAB.
    cell_of = np.searchsorted(edges, x, side="left") - 1

    order = np.argsort(cell_of, kind="stable")
    cell_sorted = cell_of[order]
    y_sorted = y[order]
    starts = np.searchsorted(cell_sorted, np.arange(options.n_x_cells), side="left")
    stops = np.searchsorted(cell_sorted, np.arange(options.n_x_cells), side="right")

    for cell in range(options.n_x_cells):
        y_cell = y_sorted[starts[cell]:stops[cell]] - options.lane_width / 2
        if y_cell.size > options.min_cell_samples:
            angle = 2 * np.pi * y_cell / options.lane_width
            shift[cell] = (np.arctan2(np.mean(np.sin(angle)), np.mean(np.cos(angle)))
                           * options.lane_width / 2 / np.pi)

    return DrivingLine(x_cells=edges[:-1], shift=shift)


def assign_lanes(
    record: dict,
    driving_line: DrivingLine,
    options: LaneIdentificationOptions = LaneIdentificationOptions(),
) -> tuple[np.ndarray, np.ndarray]:
    """Stage 2: corrected lateral position and filtered lane index.

    Returns ``(y_corr, lane)`` in feet and fractional lanes, matching the
    ``dataLanes`` struct fields of the MATLAB implementation.
    """
    y = np.asarray(record["y_position"], dtype=float)
    x = np.asarray(record["x_position"], dtype=float)

    if record["direction"] < 0:
        scale, offset = options.scale_west, options.offset_west
    else:
        scale, offset = options.scale_east, options.offset_east

    y_corr = scale * (y - driving_line(x)) + offset

    lane_raw = (np.abs(y_corr) - options.lane_width / 2) / options.lane_width
    lane_raw = np.clip(lane_raw, 0.0, 5.0)

    return y_corr, _median_filter(lane_raw, options.median_window)


def _median_filter(values: np.ndarray, window: int) -> np.ndarray:
    """Centered running median, with the MATLAB edge handling.

    MATLAB fills positions 1..buff-1 with the first computed value and the last
    buff positions with the last computed one, where ``buff = ceil(ww/2)``.
    """
    n = values.size
    buffer = int(np.ceil(window / 2))
    lane = np.zeros(n)

    if n > window:
        for center in range(buffer - 1, n - buffer):
            lane[center] = np.median(values[center - buffer + 1:center + buffer + 1])
        lane[:buffer - 1] = lane[buffer - 1]
        lane[n - buffer:] = lane[n - buffer - 1]
    else:
        lane[:] = np.median(values)
    return lane


def clip_lane_changes(
    record: dict,
    y_corr: np.ndarray,
    lane: np.ndarray,
    options: LaneChangeClippingOptions = LaneChangeClippingOptions(),
) -> List[dict]:
    """Split one trajectory into single-lane segments, dropping lane changes.

    Mirrors ``clip_lane_changes`` line for line, including its index arithmetic.
    Each returned segment carries the raw record's fields with ``timestamp``,
    ``x_position``, ``y_position`` (corrected), ``starting_x``, ``ending_x``,
    ``first_timestamp``, ``last_timestamp``, ``lane``, and the ``-N`` suffixed
    identifier that the released data uses.
    """
    x = np.asarray(record["x_position"], dtype=float)
    t = np.asarray(record["timestamp"], dtype=float)
    y = np.asarray(y_corr, dtype=float)
    lane = np.asarray(lane, dtype=float)

    segments: List[dict] = []
    trajectory_length = t.size
    pointer = 1
    clipped_part = -1

    while pointer < trajectory_length:
        if lane.size < 2:
            break
        rate = np.abs(np.diff(lane) / np.diff(t))
        stable = np.flatnonzero(rate < options.max_lane_change_rate)
        if stable.size == 0:
            break  # lane never stable: discard the remainder

        start = int(stable[0])
        lane, x, y, t = lane[start:], x[start:], y[start:], t[start:]
        temp_lane = round_half_away(lane[0])
        trajectory_length = lane.size

        changed = np.flatnonzero(np.abs(lane - temp_lane) > options.lane_change_thresh)
        # MATLAB index_change is 1-based; keep it that way for the arithmetic.
        index_change = int(changed[0]) + 1 if changed.size else trajectory_length

        mirrored_lane = lane[:index_change][::-1]
        mirrored_time = t[:index_change][::-1]
        if mirrored_lane.size >= 2:
            mirrored_rate = np.abs(np.diff(mirrored_lane) / np.diff(mirrored_time))
            mirrored_rate = np.append(mirrored_rate, mirrored_rate[-1])
        else:
            mirrored_rate = np.zeros(mirrored_lane.size)

        candidates = np.flatnonzero(
            (np.abs(mirrored_lane - temp_lane) < options.change_buffer_thresh)
            & (mirrored_rate < options.max_lane_change_rate)
        )

        if candidates.size:
            second_pointer = int(candidates[0]) + 1          # 1-based
            if (t[index_change - second_pointer] - t[0]) >= options.min_clip_time:
                clipped_part += 1
                stop = index_change - second_pointer          # 1-based count
                segments.append(_make_segment(record, t, x, y, stop, temp_lane, clipped_part))

        # Discard the processed part.
        lane = lane[index_change - 1:]
        x = x[index_change - 1:]
        y = y[index_change - 1:]
        t = t[index_change - 1:]
        temp_lane = round_half_away(lane[0])

        if index_change < trajectory_length:
            within = np.flatnonzero(np.abs(lane - temp_lane) < options.change_buffer_thresh)
            resume = int(within[0]) + 1 if within.size else lane.size

            not_full = np.flatnonzero(
                np.abs(lane - temp_lane) > 1 - options.change_buffer_thresh)
            if not_full.size and (int(not_full[0]) + 1) < resume:
                # False lane change: clip, but do not start a new trajectory.
                resume = int(not_full[0]) + 1
                temp_lane = round_half_away(lane[resume - 1])

            lane, x, y, t = lane[resume:], x[resume:], y[resume:], t[resume:]
            pointer = 1
            trajectory_length = x.size

    return segments


def _make_segment(record: dict, t: np.ndarray, x: np.ndarray, y: np.ndarray,
                  stop: int, lane_number: float, part: int) -> dict:
    segment = dict(record)
    segment["timestamp"] = t[:stop]
    segment["x_position"] = x[:stop]
    segment["y_position"] = y[:stop]
    segment["first_timestamp"] = float(t[0])
    segment["last_timestamp"] = float(t[stop - 1])
    segment["starting_x"] = float(x[0])
    segment["ending_x"] = float(x[stop - 1])
    segment["lane"] = float(lane_number)

    identifier = record["_id"]
    oid = identifier["$oid"] if isinstance(identifier, dict) else identifier
    segment["trajectory_id"] = f"{oid}-{part}"
    return segment
