"""Parse CIRCLES control-vehicle GPS files into individual runs.

Port of ``parse_gps_data`` in ``Scripts/assemble_data_GPS.m``, the first stage
of GPS assembly. Each vehicle's 10 Hz recording for a day is split into the
separate *runs* it drove through the I-24 MOTION testbed: a run begins when the
vehicle moves inward across the roadway edge within the testbed bounds, and
ends when it leaves those bounds. Runs that are too short in time or distance,
or recorded while the GPS unit was unhealthy, are discarded.

The remaining stages of ``assemble_data_GPS`` (10 Hz resampling, ping merge,
lane assignment, and matching against MOTION trajectories) are not ported yet.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np

__all__ = ["GpsRunOptions", "GpsRun", "parse_gps_data", "day_time_limits"]

#: Vehicle numbers the MATLAB loop scans (``for avID = 1:103``).
MAX_VEHICLE_ID = 103


@dataclass(frozen=True)
class GpsRunOptions:
    """Constants from the parameter block of ``parse_gps_data``."""

    max_rcs_x: float = 2.91e4      # [ft] maximum x position for testbed
    min_rcs_x: float = -3000.0     # [ft] minimum x position for testbed
    max_rcs_y: float = 150.0       # [ft] maximum |y| before the run is over
    max_gps_inactive: float = 0.01  # max fraction of samples with a bad fix
    min_run_length_km: float = 1.2  # [km] minimum extent along x
    min_run_time: float = 60.0      # [s] minimum duration
    km_to_ft: float = 3281.0
    #: Status bytes above this value mark an unhealthy GPS fix.
    status_bad_above: int = 80


@dataclass
class GpsRun:
    """One vehicle run through the testbed, in the MOTION coordinate frame."""

    vin: int                      # vehicle number (car N), not the VIN string
    run_num: int
    timestamp: np.ndarray
    direction: float              # +1 eastbound, -1 westbound
    x_position: np.ndarray        # [ft]
    y_position: np.ndarray        # [ft], sign corrected
    starting_time: float
    ending_time: float
    longitude: np.ndarray
    latitude: np.ndarray
    state_x: np.ndarray
    state_y: np.ndarray
    can_speed: np.ndarray
    control_active: np.ndarray    # bool


def day_time_limits(day: int, hours=(3, 18)) -> tuple:
    """POSIX bounds of a processing day in Nashville local time.

    Mirrors ``datetime(2022,11,day,[3 18],0,0,'TimeZone','America/Chicago')``
    followed by ``convertTo(...,'epochtime')``.
    """
    from datetime import datetime
    from zoneinfo import ZoneInfo

    zone = ZoneInfo("America/Chicago")
    return tuple(
        datetime(2022, 11, day, hour, 0, 0, tzinfo=zone).timestamp() for hour in hours)


def parse_gps_data(gps_folder, day: int,
                   options: GpsRunOptions = GpsRunOptions()) -> List[GpsRun]:
    """Split every vehicle's recording for one day into testbed runs."""
    import pandas as pd

    gps_folder = Path(gps_folder)
    lower, upper = day_time_limits(day)
    runs: List[GpsRun] = []

    for vehicle_id in range(1, MAX_VEHICLE_ID + 1):
        path = gps_folder / f"circles_v2_1_car{vehicle_id}.csv"
        if not path.is_file():
            continue

        # float_precision='round_trip' selects pandas' correctly-rounded (IEEE
        # nearest) float parser, which matches MATLAB's readtable bit-for-bit.
        # The default parser is off by 1 ULP on some values, which flips
        # 6th-decimal roundings in the output. See docs/PYTHON_PORT.md.
        table = pd.read_csv(path, dtype={"control_active": "string"},
                            float_precision="round_trip")
        table = table[(table["Systime"] > lower) & (table["Systime"] < upper)]
        if table.empty:
            continue

        runs.extend(_split_runs(table.reset_index(drop=True), vehicle_id, options))

    return runs


def _split_runs(table, vehicle_id: int, options: GpsRunOptions) -> List[GpsRun]:
    """Walk one vehicle's day, emitting each qualifying run.

    Mirrors the while-loop in parse_gps_data: locate a run start (moving toward
    the roadway centre while inside the testbed), then the first sample outside
    the testbed, evaluate the run, discard it, and repeat.
    """
    runs: List[GpsRun] = []
    run_num = 0

    rows = table
    start = _find_run_start(rows, options)
    if start is None:
        return runs
    rows = rows.iloc[start:].reset_index(drop=True)
    end = _find_run_end(rows, options)

    while len(rows) > 1:
        segment = rows.iloc[:end]
        time = segment["Systime"].to_numpy(dtype=float)
        x = segment["rcs_x"].to_numpy(dtype=float)
        status = segment["Status"].astype(str).str[0].to_numpy()
        unhealthy = np.array([ord(character) > options.status_bad_above
                              for character in status])

        long_enough = (time[-1] - time[0]) > options.min_run_time
        far_enough = abs(x[-1] - x[0]) > options.min_run_length_km * options.km_to_ft
        healthy = unhealthy.sum() < options.max_gps_inactive * x.size

        if long_enough and far_enough and healthy:
            run_num += 1
            control = segment["control_active"].astype("string").fillna("") == "True"
            runs.append(GpsRun(
                vin=vehicle_id,
                run_num=run_num,
                timestamp=time,
                direction=float(np.sign(x[-1] - x[0])),
                x_position=x,
                y_position=-segment["rcs_y"].to_numpy(dtype=float),
                starting_time=float(time[0]),
                ending_time=float(time[-1]),
                longitude=segment["Long"].to_numpy(dtype=float),
                latitude=segment["Lat"].to_numpy(dtype=float),
                state_x=segment["state_x"].to_numpy(dtype=float),
                state_y=segment["state_y"].to_numpy(dtype=float),
                can_speed=segment["can_speed"].to_numpy(dtype=float),
                control_active=control.to_numpy(dtype=bool),
            ))

        if end < len(rows):
            rows = rows.iloc[end:].reset_index(drop=True)
            start = _find_run_start(rows, options)
            if start is None:
                break
            rows = rows.iloc[start:].reset_index(drop=True)
            end = _find_run_end(rows, options)
            if end is None:
                end = len(rows)
        else:
            break

    return runs


def _find_run_start(rows, options: GpsRunOptions) -> Optional[int]:
    """First sample moving toward the roadway centre, inside the testbed."""
    y = rows["rcs_y"].to_numpy(dtype=float)
    x = rows["rcs_x"].to_numpy(dtype=float)
    if y.size == 0:
        return None

    difference = np.diff(y)
    difference = np.concatenate((difference[:1], difference)) if difference.size else np.zeros(1)

    entering = (np.sign(difference) * np.sign(y) < 0) & (np.abs(y) < options.max_rcs_y) \
        & (x > options.min_rcs_x) & (x < options.max_rcs_x)
    found = np.flatnonzero(entering)
    return int(found[0]) if found.size else None


def _find_run_end(rows, options: GpsRunOptions) -> int:
    """First sample outside the testbed bounds, or the end of the data."""
    y = rows["rcs_y"].to_numpy(dtype=float)
    x = rows["rcs_x"].to_numpy(dtype=float)
    outside = (np.abs(y) > options.max_rcs_y) | (x < options.min_rcs_x) | (x > options.max_rcs_x)
    found = np.flatnonzero(outside)
    return int(found[0]) if found.size else len(rows)
