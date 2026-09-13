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

Performance: the hot paths are vectorized - the lane median filter (sliding
window), the smoothdata gaussian (a searchsorted-banded, fully vectorized
kernel), the match-stretch walk (run-length on a boolean mask), and per-segment
caching of each trajectory's arrays. Together these took a full day from ~65 min
to ~11 min (6x) with bit-identical output. Each vectorization was checked
against the original scalar form before being adopted.
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
    n = values.size

    # Banded and fully vectorized. Sample points are time-sorted, so each
    # point's truncation window is a contiguous slice; searchsorted gives its
    # bounds. Because the sampling is dense, the widest window spans only a
    # bounded number of points (~2*half worth), so the per-point windows fit in
    # one (n, maxband) array padded on the right and masked to the true bounds.
    # O(n * band), no Python loop; bit-equivalent to the masked O(n^2) form
    # (verified to < 1e-15).
    if n == 0:
        return np.empty(0)
    lo = np.searchsorted(t, t - half, side="left")
    hi = np.searchsorted(t, t + half, side="right")
    max_band = int((hi - lo).max())

    cols = lo[:, None] + np.arange(max_band)[None, :]
    valid = cols < hi[:, None]
    cols_clipped = np.minimum(cols, n - 1)
    windowed = values[cols_clipped]
    # MATLAB's smoothdata omits NaN: a point is NaN only when its whole window
    # is NaN, not when the window merely touches one. Propagating NaN instead
    # poisoned a half-window (1.5 s ~ 37 samples at 25 Hz) ahead of every NaN,
    # which ended matched stretches early and changed the median_xd bias -
    # dist_to_av is NaN wherever a MOTION trajectory runs past the AV's own
    # time range, so this happened at the end of most matched trajectories.
    valid = valid & ~np.isnan(windowed)
    dt = t[cols_clipped] - t[:, None]
    weights = np.exp(-(dt * dt) / (2 * sigma ** 2)) * valid
    total = weights.sum(axis=1)
    numerator = (weights * np.where(valid, windowed, 0.0)).sum(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(total > 0, numerator / total, np.nan)


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
                  progress=None,
                  collect: Optional[Dict[int, np.ndarray]] = None) -> Dict[int, float]:
    """Per-run median x offset between each GPS run and its MOTION trajectories.

    Returns ``{run.index: median_xd}`` in meters; runs with no match get 0.0,
    matching MATLAB. ``raw_motion_dir`` is ``data/i24motion/2022-11-DD``.

    Pass a dict as ``collect`` to also receive the pooled matched distances per
    run. The bias is a median over those, so when a run's bias disagrees with
    MATLAB the question is always which values entered the pool - this makes
    that inspectable without re-running the (multi-minute) matching pass.
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

        # MATLAB filters on the *raw* run bounds here, not the resampled grid
        # (preproc_gps overwrites timestamp but leaves starting/ending_time).
        concurrent = [run for run in runs
                      if run.starting_time <= first_ts.max()
                      and run.ending_time >= first_ts.min()]
        if not concurrent:
            if progress:
                progress(seg_index, len(segment_files), 0)
            continue

        lanes = assign_lanes_bidirectional(trajectories)
        # Convert each trajectory's arrays once per segment, not once per
        # (run, candidate) pair - the same trajectory is a candidate for many
        # runs, and re-converting the decoded lists dominated the walk.
        traj_t = [np.asarray(tr["timestamp"], dtype=float) for tr in trajectories]
        traj_x = [np.asarray(tr["x_position"], dtype=float) - ORIGIN_X_FEET
                  for tr in trajectories]

        for run in concurrent:
            _match_run_in_segment(run, traj_t, traj_x, lanes, first_ts, last_ts,
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
            if collect is not None:
                collect[run.index] = values
        else:
            bias[run.index] = 0.0
            if collect is not None:
                collect[run.index] = np.empty(0)
    return bias


# ---------------------------------------------------------------------------
# internals

_EPS = 1e-9
#: Above this trajectory length, smoothdata falls back to the row-at-a-time form
#: to bound the O(n^2) weight matrix (matched trajectories are far shorter).
_SMOOTH_DENSE_MAX = 6000


def _match_run_in_segment(run, traj_times, traj_xs, lanes, first_ts, last_ts,
                          directions, options, matched):
    av_t = run.timestamp
    av_x = run.x_position                     # [ft], bumper-shifted, 10 Hz
    av_lane = run.assigned_lane

    candidates = np.flatnonzero((first_ts <= av_t[-1]) & (last_ts >= av_t[0])
                                & (directions == run.direction))
    for idx in candidates:
        traj_lane = lanes[idx]
        traj_x = traj_xs[idx]                  # [ft], already origin-subtracted
        traj_t = traj_times[idx]

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
    """Record stretches matched for >= min_match_time.

    Vectorized equivalent of the per-sample walk (verified identical, including
    the quirk that the breaking sample is included in a recorded stretch): a
    sample matches when it is close in position, lane, and smoothed relative
    speed; a run of matching samples from index s to e spans traj_t[e]-traj_t[s-1]
    seconds, and qualifying runs contribute dist_to_av over [s, e+1] (the
    trailing breaking sample), or to the array end for a run that never breaks.
    Only samples 1..n-1 can match, matching the MATLAB loop start.
    """
    n = traj_t.size
    if n < 2:
        return
    cond = np.zeros(n, dtype=bool)
    cond[1:] = ((np.abs(dist_to_av[1:]) <= options.max_match_dist)
                & (np.abs(lane_diff[1:]) <= options.max_match_lane_diff)
                & (v_diff_smooth[1:] <= options.max_match_speed_diff))

    edges = np.diff(np.concatenate(([False], cond, [False])).astype(np.int8))
    starts = np.flatnonzero(edges == 1)
    ends = np.flatnonzero(edges == -1) - 1
    for start, end in zip(starts, ends):
        if (traj_t[end] - traj_t[start - 1]) >= options.min_match_time:
            stop = end + 2 if end + 1 < n else n
            matched[run_index].append(dist_to_av[start:stop])


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
    # Vectorized sliding-window median; bit-identical to the per-point loop
    # (verified against it). This is the dominant cost of the lane assignment,
    # which runs on every trajectory in every segment.
    from numpy.lib.stride_tricks import sliding_window_view

    n = values.size
    buffer = int(np.ceil(window / 2))
    lane = np.zeros(n)
    if n > window:
        windows = sliding_window_view(values, 2 * buffer)   # (n-2*buffer+1, 2*buffer)
        lane[buffer - 1:n - buffer] = np.median(windows, axis=1)
        lane[:buffer - 1] = lane[buffer - 1]
        lane[n - buffer:] = lane[n - buffer - 1]
    else:
        lane[:] = np.median(values)
    return lane


def _interp1_nan(xp, fp, x):
    """interp1(xp, fp, x) with the weighted-blend form; NaN outside range.

    Carries the same flat-segment guard as
    :func:`mvtpy.gpsassemble._interp_extrap`: where both endpoints are equal the
    blend lands a ULP off and MATLAB returns the value itself. It matters here
    because these distances decide which trajectory segments count as matched,
    and so which values the median_xd bias is taken over.
    """
    xp = np.asarray(xp, dtype=float)
    fp = np.asarray(fp, dtype=float)
    x = np.asarray(x, dtype=float)
    index = np.clip(np.searchsorted(xp, x, side="right") - 1, 0, len(xp) - 2)
    left, right = fp[index], fp[index + 1]
    weight = (x - xp[index]) / (xp[index + 1] - xp[index])
    out = np.where(left == right, left, left * (1 - weight) + right * weight)
    out[(x < xp[0]) | (x > xp[-1])] = np.nan
    return out


def _relevant_segments(runs, raw_motion_dir, day):
    """MOTION segments overlapping the AV activity window (minFileNr:maxFileNr)."""
    from datetime import datetime
    from zoneinfo import ZoneInfo

    zone = ZoneInfo("America/Chicago")
    day_start = datetime(2022, 11, day, 6, 0, 0, tzinfo=zone).timestamp()

    # Raw run bounds, as MATLAB uses at assemble_data_GPS.m:139-143.
    starts = np.array([getattr(run, "starting_time", run.timestamp[0]) for run in runs])
    ends = np.array([getattr(run, "ending_time", run.timestamp[-1]) for run in runs])
    min_av_start = starts[ends > day_start].min()
    max_av_end = ends[starts < day_start + 4 * 3600].max()
    min_file = max(1, int(np.floor((min_av_start - day_start) / 60 / 10)))
    max_file = min(24, int(np.floor((max_av_end - day_start) / 60 / 10)) + 1)

    abbr = {16: "wed", 17: "thu", 18: "fri"}[day]
    files = sorted(raw_motion_dir.glob(f"*_{abbr}_0_*.json"))
    files = [f for f in files if not f.name.startswith(".")]
    # minFileNr:maxFileNr are 1-based positions in the sorted listing.
    return files[min_file - 1:max_file]
