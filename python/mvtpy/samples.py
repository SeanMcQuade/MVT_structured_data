"""Collect per-sample data near engaged control vehicles.

Port of ``generate_data_samples.m``. For one day, it walks the processed
(westbound) MOTION trajectories and, for every sample within 1000 m of an
*engaged* control vehicle, records the distance to that vehicle together with
speed, fuel rate, fuel consumption, vehicle class, x position, lane, and time
after 06:00. The result is the ``.mat`` that ``plot_AV_analysis`` bins to
produce the article's fuel-vs-distance figures (Figure 2, SM2, SM3).

A single sample can contribute twice — once for a downstream engaged AV and once
for an upstream one — and the downstream contribution is appended first, matching
MATLAB's inner loop order.

Output dtypes mirror the MATLAB `.mat` exactly: distance/speed/fuel-rate/fuel-
consumption are double, x position is int16, class and lane are uint8, and time
is uint16. MATLAB's integer casts round to nearest (ties away from zero) and
saturate; ``matlab_int`` reproduces that.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable, List

import numpy as np

from .matround import round_half_away

__all__ = ["MAX_DIST_M", "matlab_int", "collect_samples", "collect_samples_from_dir"]

#: [m] maximum distance from an engaged AV for a sample to be collected.
MAX_DIST_M = 1000.0

#: [s] epoch time of 06:00 on 2022-11-18; other days offset by whole days.
SIX_AM_NOV18 = 1668772800

#: Output field -> numpy dtype, matching the MATLAB .mat.
SAMPLE_DTYPES = {
    "samples_dist": np.float64,
    "samples_speed": np.float64,
    "samples_fr": np.float64,
    "samples_fcons": np.float64,
    "samples_class": np.uint8,
    "samples_xpos": np.int16,
    "samples_lane": np.uint8,
    "samples_t": np.uint16,
}


def matlab_int(values, dtype) -> np.ndarray:
    """Cast to an integer dtype the way MATLAB does: round half away, saturate.

    NumPy truncates toward zero and wraps on overflow; MATLAB rounds to nearest
    (ties away from zero) and clamps to the type's range. The x-position and
    time fields depend on this.
    """
    info = np.iinfo(dtype)
    values = np.asarray(values, dtype=float)
    rounded = np.where(values >= 0, np.floor(values + 0.5), np.ceil(values - 0.5))
    return np.clip(rounded, info.min, info.max).astype(dtype)


def collect_samples(records: Iterable[dict], day: int) -> Dict[str, np.ndarray]:
    """Aggregate samples from an iterable of processed (slim) trajectory records.

    ``records`` may span several files; pass them in file order to reproduce the
    MATLAB sample order. Each record is the released slim schema (a dict with the
    distance-to-AV and per-sample fields).
    """
    six_am_day = SIX_AM_NOV18 - (18 - day) * 24 * 60 * 60

    dist, speed, fuel_rate = [], [], []
    vehicle_class, x_pos, lane, time_after = [], [], [], []

    for veh in records:
        if veh.get("direction", -1) > 0:      # westbound only (slim has no eastbound)
            continue
        downstream = veh["distance_to_downstream_engaged_av_meters"]
        upstream = veh["distance_to_upstream_engaged_av_meters"]
        has_down = _nonempty(downstream)
        has_up = _nonempty(upstream)
        if not (has_down or has_up):
            continue

        timestamps = np.asarray(veh["timestamp"], dtype=float)
        speeds = np.asarray(veh["speed_meters_per_second"], dtype=float)
        rates = np.asarray(veh["fuel_rate_grams_per_second"], dtype=float)
        positions = np.asarray(veh["x_position_meters"], dtype=float)
        klass = int(veh["coarse_vehicle_class"])
        lane_number = veh["lane_number"]

        down = _as_array(downstream) if has_down else None
        up = _as_array(upstream) if has_up else None
        down_hits = (np.flatnonzero(~np.isnan(down) & (down <= MAX_DIST_M))
                     if has_down else np.empty(0, dtype=int))
        up_hits = (np.flatnonzero(~np.isnan(up) & (up >= -MAX_DIST_M))
                   if has_up else np.empty(0, dtype=int))
        if down_hits.size == 0 and up_hits.size == 0:
            continue

        # Preserve MATLAB's per-point order: at each sample the downstream
        # contribution is appended before the upstream one. Keying downstream
        # hits by 2*i and upstream by 2*i+1 and stable-sorting reproduces that
        # interleaving without a per-point Python loop.
        keys = np.concatenate((2 * down_hits, 2 * up_hits + 1))
        order = np.argsort(keys, kind="stable")
        point = np.concatenate((down_hits, up_hits))[order]
        distances = np.concatenate((
            down[down_hits] if has_down else np.empty(0),
            up[up_hits] if has_up else np.empty(0)))[order]

        dist.append(distances)
        speed.append(speeds[point])
        fuel_rate.append(rates[point])
        x_pos.append(positions[point])
        time_after.append(timestamps[point] - six_am_day)
        vehicle_class.append(np.full(point.size, klass))
        lane.append(np.full(point.size, lane_number))

    samples_dist = _cat(dist)
    samples_speed = _cat(speed)
    samples_fr = _cat(fuel_rate)
    # fuel consumption is computed after the integer casts, on full precision.
    samples_fcons = samples_fr / (1e-6 + samples_speed)

    return {
        "samples_dist": samples_dist,
        "samples_speed": samples_speed,
        "samples_fr": samples_fr,
        "samples_fcons": samples_fcons,
        "samples_class": matlab_int(_cat(vehicle_class), np.uint8),
        "samples_xpos": matlab_int(_cat(x_pos), np.int16),
        "samples_lane": matlab_int(_cat(lane), np.uint8),
        "samples_t": matlab_int(_cat(time_after), np.uint16),
    }


def collect_samples_from_dir(slim_dir, day: int) -> Dict[str, np.ndarray]:
    """Collect samples from a day's processed JSON files, in filename order."""
    from .rawio import iter_trajectories

    slim_dir = Path(slim_dir)
    files = sorted(slim_dir.glob(f"I-24MOTION_2022-11-{day}_*.json"))

    def all_records():
        for path in files:
            records = list(iter_trajectories(path))
            # "no AVs on the road" - the whole file is skipped when the field
            # is absent (MATLAB checks isfield on the first record's struct).
            if records and "distance_to_upstream_engaged_av_meters" not in records[0]:
                continue
            yield from records

    return collect_samples(all_records(), day)


# ---------------------------------------------------------------------------


def _cat(chunks: List[np.ndarray]) -> np.ndarray:
    return np.concatenate(chunks) if chunks else np.empty(0)


def _nonempty(value) -> bool:
    if value is None:
        return False
    return np.atleast_1d(np.asarray(value, dtype=object)).size > 0


def _as_array(value) -> np.ndarray:
    return np.asarray([np.nan if v is None else v for v in np.atleast_1d(value)],
                      dtype=float)
