"""Distance from each trajectory to the nearest CIRCLES control vehicle.

Port of ``calculate_distance_to_avs`` in ``Scripts/generate_data_mvt_slim.m``.
For every sample of every trajectory it finds the nearest control vehicle ahead
(downstream) and behind (upstream), both unconditionally and restricted to
vehicles whose controller was engaged, and reports the distance in meters plus
the vehicle's id.

Sign convention, as in the released data: downstream distances are positive,
upstream distances are negative (the caller negates the magnitude computed
here). Where no vehicle qualifies at a sample, the value is NaN, which
``jsonencode`` writes as ``null``; where none qualifies anywhere in the
trajectory, the field is empty.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence

import numpy as np

from .kinematics import FT_TO_METER, ORIGIN_X_FEET
from .matround import round_half_away

__all__ = ["AvRun", "load_av_runs", "interp1_nan_outside", "distance_to_avs"]


@dataclass
class AvRun:
    """One control-vehicle run from the assembled GPS file."""

    av_id: float
    assigned_lane: float
    direction: float
    timestamp: np.ndarray
    x_position: np.ndarray
    engaged: np.ndarray          # control_car when present, else controller_engaged
    first_timestamp: float
    last_timestamp: float

    @classmethod
    def from_record(cls, record: dict) -> "AvRun":
        """Build an AvRun from a released GPS record."""
        timestamp = np.asarray(record["timestamp"], dtype=float)
        # MATLAB prefers control_car and falls back to controller_engaged when
        # that field is empty.
        control = record.get("control_car")
        if control is None or len(np.atleast_1d(control)) == 0:
            control = record["controller_engaged"]
        return cls(
            av_id=float(record["av_id"]),
            assigned_lane=float(record["assigned_lane"]),
            direction=float(record["direction"]),
            timestamp=timestamp,
            x_position=np.asarray(record["x_position"], dtype=float),
            engaged=np.asarray(control, dtype=float),
            first_timestamp=float(record["first_timestamp"]),
            last_timestamp=float(record["last_timestamp"]),
        )


def load_av_runs(path, direction: Optional[float] = None) -> List[AvRun]:
    """Read the assembled GPS file (``results/gps/CIRCLES_GPS_10Hz_*.json``)."""
    from .rawio import iter_trajectories

    runs = []
    for record in iter_trajectories(path):
        if direction is not None and record["direction"] != direction:
            continue
        runs.append(AvRun.from_record(record))
    return runs


def interp1_nan_outside(xp: np.ndarray, fp: np.ndarray, x: np.ndarray) -> np.ndarray:
    """``interp1(xp, fp, x)`` - linear, NaN outside the sample range.

    NumPy's ``interp`` clamps to the end values instead, which would invent a
    control vehicle at times when it was not on the road.
    """
    values = np.interp(x, xp, fp)
    outside = (x < xp[0]) | (x > xp[-1])
    values[outside] = np.nan
    return values


def distance_to_avs(segment: dict, runs: Sequence[AvRun]) -> Dict[str, Optional[np.ndarray]]:
    """Nearest upstream/downstream control vehicle for one clipped trajectory.

    Parameters
    ----------
    segment:
        A clipped trajectory as produced by :func:`mvtpy.lanes.clip_lane_changes`:
        ``timestamp`` (POSIX seconds), ``x_position`` (feet), ``lane``,
        ``direction``.
    runs:
        Control-vehicle runs from :func:`load_av_runs`.

    Returns
    -------
    dict with keys ``distanceToAvUS``, ``AvIdUS``, ``distanceToAvUSEng``,
    ``AvIdUSEng`` and the ``DS`` equivalents. Distances are unsigned magnitudes,
    matching the MATLAB helper; the caller applies the sign. A value of ``None``
    means the field is empty (no qualifying vehicle anywhere).
    """
    time = np.asarray(segment["timestamp"], dtype=float)
    x_meters = FT_TO_METER * (np.asarray(segment["x_position"], dtype=float) - ORIGIN_X_FEET)
    lane = float(segment["lane"])
    direction = float(segment["direction"])

    on_road = [run for run in runs
               if run.assigned_lane == lane
               and run.direction == direction
               and time[0] < run.last_timestamp
               and time[-1] > run.first_timestamp]

    result: Dict[str, Optional[np.ndarray]] = {
        "distanceToAvUS": None, "AvIdUS": None,
        "distanceToAvUSEng": None, "AvIdUSEng": None,
        "distanceToAvDS": None, "AvIdDS": None,
        "distanceToAvDSEng": None, "AvIdDSEng": None,
    }
    if not on_road:
        return result

    # Signed distance to each candidate vehicle at each sample, and the same
    # restricted to samples where its controller was engaged.
    n_runs, n_samples = len(on_road), time.size
    signed = np.empty((n_runs, n_samples))
    engaged_mask = np.zeros((n_runs, n_samples), dtype=bool)

    for index, run in enumerate(on_road):
        projected = interp1_nan_outside(run.timestamp, run.x_position, time)
        signed[index] = (projected - x_meters) * direction

        active = interp1_nan_outside(run.timestamp, run.engaged, time)
        with np.errstate(invalid="ignore"):
            active = np.clip(round_half_away_array(active), 0, 1)
        engaged_mask[index] = active == 1

    for sign, suffix in ((-1.0, "US"), (1.0, "DS")):
        # Wrong side, or no overlap, means "infinitely far" in this direction.
        distance = np.where(np.isnan(signed), sign * np.inf, signed)
        distance = np.where(sign * distance < 0, sign * np.inf, distance)

        engaged_distance = np.where(engaged_mask, distance, sign * np.inf)

        _select_nearest(result, distance, on_road, sign, f"distanceToAv{suffix}",
                        f"AvId{suffix}")
        _select_nearest(result, engaged_distance, on_road, sign,
                        f"distanceToAv{suffix}Eng", f"AvId{suffix}Eng")

    return result


def _select_nearest(result: Dict[str, Optional[np.ndarray]], distance: np.ndarray,
                    runs: Sequence[AvRun], sign: float,
                    distance_key: str, id_key: str) -> None:
    """``min(sn*D, [], 1)`` plus MATLAB's empty/NaN bookkeeping."""
    oriented = sign * distance
    nearest_index = np.argmin(oriented, axis=0)
    nearest = oriented[nearest_index, np.arange(oriented.shape[1])]
    identifiers = np.array([runs[index].av_id for index in nearest_index], dtype=float)

    infinite = np.isinf(nearest)
    if infinite.all():
        result[distance_key] = None
        result[id_key] = None
        return

    nearest = nearest.astype(float)
    nearest[infinite] = np.nan
    identifiers[infinite] = np.nan
    result[distance_key] = nearest
    result[id_key] = identifiers


def round_half_away_array(values: np.ndarray) -> np.ndarray:
    """Vectorized MATLAB ``round`` that leaves NaN untouched."""
    values = np.asarray(values, dtype=float)
    finite = np.isfinite(values)
    out = np.array(values, dtype=float)
    out[finite] = np.where(values[finite] >= 0,
                           np.floor(values[finite] + 0.5),
                           np.ceil(values[finite] - 0.5))
    return out
