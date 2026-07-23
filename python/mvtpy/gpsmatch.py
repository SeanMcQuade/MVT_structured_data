"""MOTION-matching bias for control-vehicle x_position.

Final piece of ``assemble_data_GPS.m``: each GPS run's ``x_position`` is
corrected by a per-run median offset (``median_xd``) measured against the I-24
MOTION trajectories the vehicle was observed as. The offset removes a systematic
bias between the on-vehicle GPS and the MOTION observatory (median ~1.5 m,
up to ~5 m across runs).

Computing it means, for each MOTION segment overlapping the day's AV activity:
decode the segment, assign lanes to its trajectories, and for every
(run, candidate trajectory) pair measure the distance between the GPS run and
the trajectory, keeping stretches that stay close in position, lane, and speed
for at least ``min_match_time`` seconds. The median of those distances over all
matched stretches of a run is its bias.

This is the heaviest stage in the port: it re-reads the day's raw MOTION
segments (streamed) and reproduces MATLAB's `smoothdata` gaussian and the
both-direction lane assignment.

Performance note: this implementation is faithful but pure-Python and slow -
tens of minutes for a full day, dominated by the per-trajectory lane median
filter and per-candidate smoothing over the peak-activity segments. It is
correct as written; vectorizing the lane filter and the smoothing would be the
obvious optimization before using it in a production pipeline.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np

from .kinematics import FT_TO_METER, ORIGIN_X_FEET
from .lanes import LaneIdentificationOptions, interp1_linear_extrap
from .gpsassemble import PreprocessedRun
from .rawio import iter_trajectories

__all__ = [
    "MatchOptions",
    "smoothdata_gaussian",
    "assign_lanes_bidirectional",
    "matching_bias",
]


@dataclass(frozen=True)
class MatchOptions:
    """Thresholds from the head of assemble_data_GPS.m."""

    max_match_dist: float = 6.0        # [m]
    max_match_speed_diff: float = 2.0  # [m/s]
    max_match_lane_diff: float = 0.5   # [lane]
    min_match_time: float = 3.0        # [s]
    # Candidate pre-filter (looser than the per-sample thresholds above).
    max_avg_dist: float = 200.0        # [m]
    max_avg_lane_diff: float = 2.0     # [lane]


def smoothdata_gaussian(values: np.ndarray, sample_points: np.ndarray,
                        window: float) -> np.ndarray:
    """Reproduce MATLAB ``smoothdata(values,'gaussian',window,'SamplePoints',t)``.

    MATLAB's gaussian window uses standard deviation ``window/5`` and truncates
    at ``window/2`` (2.5 sigma), then normalizes the weights per point (so edge
    points, whose window is clipped, still sum to one). Verified against a delta
    probe of MATLAB smoothdata.
    """
    values = np.asarray(values, dtype=float)
    t = np.asarray(sample_points, dtype=float)
    sigma = window / 5.0
    half = window / 2.0

    out = np.empty(values.size)
    for i in range(values.size):
        dt = t - t[i]
        within = np.abs(dt) <= half + _EPS
        weights = np.exp(-(dt[within] ** 2) / (2 * sigma ** 2))
        out[i] = np.dot(weights, values[within]) / weights.sum()
    return out


def assign_lanes_bidirectional(trajectories: Sequence[dict],
                               options: LaneIdentificationOptions = LaneIdentificationOptions()
                               ) -> List[np.ndarray]:
    """Lane index per timestep for each MOTION trajectory, both directions.

    The GPS-side ``assign_lanes`` (unlike the westbound-only one in
    ``mvtpy.lanes``) builds a driving line for each direction and corrects
    eastbound trajectories with a sign flip. Returns one lane array per input
    trajectory, in order.
    """
    lane_width = options.lane_width
    upper = lane_width * options.y_up_lim
    lower = lane_width * options.y_low_lim

    west_x, west_y, east_x, east_y = [], [], [], []
    for veh in trajectories:
        t = np.asarray(veh["timestamp"], dtype=float)
        if t.size < 2 or (t[-1] - t[0]) <= 5:
            continue
        step = int(np.ceil(1.0 / (t[1] - t[0])))
        x = np.asarray(veh["x_position"], dtype=float)[::step]
        y = np.asarray(veh["y_position"], dtype=float)[::step]
        if veh["direction"] == -1:
            west_x.append(x); west_y.append(y)
        else:
            east_x.append(x); east_y.append(y)

    west_line = _driving_line(west_x, west_y, lower, upper, options, flip=False)
    east_line = _driving_line(east_x, east_y, lower, upper, options, flip=True)

    lanes: List[np.ndarray] = []
    for veh in trajectories:
        y = np.asarray(veh["y_position"], dtype=float)
        x = np.asarray(veh["x_position"], dtype=float)
        if veh["direction"] > 0:
            y = -y
            corrected = -(options.scale_east
                          * (y - east_line(x)) + options.offset_east)
        else:
            corrected = (options.scale_west
                         * (y - west_line(x)) + options.offset_west)
        lane_raw = np.clip((np.abs(corrected) - lane_width / 2) / lane_width, 0.0, 5.0)
        lanes.append(_median_filter(lane_raw, options.median_window))
    return lanes


def matching_bias(runs: Sequence[PreprocessedRun], raw_motion_dir,
                  day: int, options: MatchOptions = MatchOptions(),
                  segment_files: Optional[Sequence] = None,
                  progress=None) -> Dict[int, float]:
    """Per-run median x offset between each GPS run and its MOTION trajectories.

    Returns ``{run.index: median_xd}`` in meters; runs with no match get 0.0,
    matching MATLAB. ``raw_motion_dir`` is ``data/i24motion/2022-11-DD``.
    """
    from .gpsruns import GpsRunOptions  # noqa: F401  (kept for symmetry)

    raw_motion_dir = Path(raw_motion_dir)
    if segment_files is None:
        segment_files = _relevant_segments(runs, raw_motion_dir, day)

    # Accumulate matched distances per run index.
    matched: Dict[int, List[np.ndarray]] = {run.index: [] for run in runs}

    for seg_index, seg_path in enumerate(segment_files):
        trajectories = list(iter_trajectories(seg_path))
        if not trajectories:
            continue
        first_ts = np.array([tr["first_timestamp"] for tr in trajectories])
        last_ts = np.array([tr["last_timestamp"] for tr in trajectories])
        directions = np.array([tr["direction"] for tr in trajectories])

        concurrent = [run for run in runs
                      if run.timestamp[0] <= first_ts.max()
                      and run.timestamp[-1] >= first_ts.min()]
        if not concurrent:
            if progress:
                progress(seg_index, len(segment_files), 0)
            continue

        lanes = assign_lanes_bidirectional(trajectories)

        for run in concurrent:
            _match_run_in_segment(run, trajectories, lanes, first_ts, last_ts,
                                  directions, options, matched)
        if progress:
            progress(seg_index, len(segment_files), len(concurrent))

    bias: Dict[int, float] = {}
    for run in runs:
        pooled = matched[run.index]
        if pooled:
            values = np.concatenate(pooled)
            values = values[~np.isnan(values)]
            bias[run.index] = float(np.median(values)) if values.size else 0.0
        else:
            bias[run.index] = 0.0
    return bias


# ---------------------------------------------------------------------------
# internals

_EPS = 1e-9


def _match_run_in_segment(run, trajectories, lanes, first_ts, last_ts,
                          directions, options, matched):
    av_t = run.timestamp
    av_x = run.x_position                     # [ft], bumper-shifted, 10 Hz
    av_lane = run.assigned_lane

    candidates = np.flatnonzero((first_ts <= av_t[-1]) & (last_ts >= av_t[0])
                                & (directions == run.direction))
    for idx in candidates:
        traj = trajectories[idx]
        traj_lane = lanes[idx]
        traj_x = np.asarray(traj["x_position"], dtype=float) - ORIGIN_X_FEET
        traj_t = np.asarray(traj["timestamp"], dtype=float)

        # interp1 default: linear, NaN outside the AV's own time range.
        av_at_traj = _interp1_nan(av_t, av_x, traj_t)
        dist_to_av = (av_at_traj - traj_x) * FT_TO_METER

        finite = dist_to_av[~np.isnan(dist_to_av)]
        if finite.size == 0 or np.mean(np.abs(finite)) > options.max_avg_dist:
            continue
        lane_diff = av_lane - traj_lane
        finite_lane = lane_diff[~np.isnan(lane_diff)]
        if finite_lane.size == 0 or np.mean(np.abs(finite_lane)) >= options.max_avg_lane_diff:
            continue

        v_diff = np.diff(dist_to_av) / np.diff(traj_t)
        v_diff = np.concatenate((v_diff[:1], v_diff))
        v_diff_smooth = np.abs(smoothdata_gaussian(v_diff, traj_t, 3.0))

        _collect_matches(run.index, dist_to_av, lane_diff, v_diff_smooth, traj_t,
                         options, matched)


def _collect_matches(run_index, dist_to_av, lane_diff, v_diff_smooth, traj_t,
                     options, matched):
    """Walk the trajectory, emitting stretches matched for >= min_match_time."""
    time_matching = 0.0
    match_start = None
    for step in range(1, traj_t.size):        # MATLAB stepInTraj = 2..end
        ok = (abs(dist_to_av[step]) <= options.max_match_dist
              and abs(lane_diff[step]) <= options.max_match_lane_diff
              and v_diff_smooth[step] <= options.max_match_speed_diff)
        if ok:
            if time_matching == 0:
                match_start = step
            time_matching += traj_t[step] - traj_t[step - 1]
        else:
            if time_matching >= options.min_match_time:
                matched[run_index].append(dist_to_av[match_start:step + 1])
            match_start = None
            time_matching = 0.0
    if time_matching >= options.min_match_time and match_start is not None:
        matched[run_index].append(dist_to_av[match_start:traj_t.size])


def _driving_line(xs, ys, lower, upper, options, flip):
    from .lanes import DrivingLine

    if not xs:
        # No trajectories in this direction: a flat, zero driving line.
        return DrivingLine(x_cells=np.array([0.0, 1.0]), shift=np.array([0.0]))
    x = np.concatenate(xs)
    y = np.concatenate(ys)
    if flip:
        y = -y
    keep = (y <= upper) & (y >= lower)
    x, y = x[keep], y[keep]

    edges = np.linspace(x.min(), x.max(), options.n_x_cells + 1)
    shift = np.zeros(options.n_x_cells)
    cell_of = np.searchsorted(edges, x, side="left") - 1
    order = np.argsort(cell_of, kind="stable")
    cs, yorder = cell_of[order], y[order]
    starts = np.searchsorted(cs, np.arange(options.n_x_cells), side="left")
    stops = np.searchsorted(cs, np.arange(options.n_x_cells), side="right")
    for cell in range(options.n_x_cells):
        y_cell = yorder[starts[cell]:stops[cell]] - options.lane_width / 2
        if y_cell.size > options.min_cell_samples:
            angle = 2 * np.pi * y_cell / options.lane_width
            shift[cell] = (np.arctan2(np.mean(np.sin(angle)), np.mean(np.cos(angle)))
                           * options.lane_width / 2 / np.pi)
    from .lanes import DrivingLine
    return DrivingLine(x_cells=edges[:-1], shift=shift)


def _median_filter(values: np.ndarray, window: int) -> np.ndarray:
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


def _interp1_nan(xp, fp, x):
    """interp1(xp, fp, x) with the weighted-blend form; NaN outside range."""
    xp = np.asarray(xp, dtype=float)
    fp = np.asarray(fp, dtype=float)
    x = np.asarray(x, dtype=float)
    index = np.clip(np.searchsorted(xp, x, side="right") - 1, 0, len(xp) - 2)
    weight = (x - xp[index]) / (xp[index + 1] - xp[index])
    out = fp[index] * (1 - weight) + fp[index + 1] * weight
    out[(x < xp[0]) | (x > xp[-1])] = np.nan
    return out


def _relevant_segments(runs, raw_motion_dir, day):
    """MOTION segments overlapping the AV activity window (minFileNr:maxFileNr)."""
    from datetime import datetime
    from zoneinfo import ZoneInfo

    zone = ZoneInfo("America/Chicago")
    day_start = datetime(2022, 11, day, 6, 0, 0, tzinfo=zone).timestamp()

    starts = np.array([run.timestamp[0] for run in runs])
    ends = np.array([run.timestamp[-1] for run in runs])
    min_av_start = starts[ends > day_start].min()
    max_av_end = ends[starts < day_start + 4 * 3600].max()
    min_file = max(1, int(np.floor((min_av_start - day_start) / 60 / 10)))
    max_file = min(24, int(np.floor((max_av_end - day_start) / 60 / 10)) + 1)

    abbr = {16: "wed", 17: "thu", 18: "fri"}[day]
    files = sorted(raw_motion_dir.glob(f"*_{abbr}_0_*.json"))
    files = [f for f in files if not f.name.startswith(".")]
    # minFileNr:maxFileNr are 1-based positions in the sorted listing.
    return files[min_file - 1:max_file]
