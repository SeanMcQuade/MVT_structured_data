"""Assemble control-vehicle GPS runs into the released per-day GPS records.

Continues the port of ``assemble_data_GPS.m`` past run-splitting (which lives in
``mvtpy.gpsruns``). This module covers the parts that turn a parsed run into a
released record:

* ``preprocess_run`` / ``resample_10hz`` - resample onto a uniform 10 Hz grid
  and shift x from the GPS antenna to the rear bumper (``preproc_gps``,
  ``sample_10hz``);
* ``connection_status`` - server connectivity per vehicle from the ping file
  (``get_connection_status``);
* ``control_car_status`` - the reconstructed control signal and its 30 s
  look-back (``get_control_car_status``);
* ``assemble_run`` - build one output record with the released field order,
  rounding, and testbed clipping.

The one piece not here is the MOTION-matching bias correction (``median_xd``),
which subtracts a per-run constant from ``x_position``. It requires decoding all
24 MOTION segments for the day; ``mvtpy.gpsmatch`` handles it. Every other field
is reproduced and verified without it.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np

from .kinematics import FT_TO_METER
from .gpsruns import GpsRun

__all__ = [
    "IN_VEHICLE_SHIFT_FT",
    "resample_10hz",
    "preprocess_run",
    "PreprocessedRun",
    "connection_status",
    "control_car_status",
    "assemble_run",
    "load_lane_map",
    "assemble_day",
]

#: [ft] x shift from the GPS antenna to the AV rear bumper (``inVehShift``).
IN_VEHICLE_SHIFT_FT = 7.0

#: [ft] testbed x limits for the final clip, in feet before conversion.
TESTBED_X_LIMITS_FT = (-400.0, 2.54e4)


def resample_10hz(original_time: np.ndarray, original_value: np.ndarray,
                  new_time: np.ndarray) -> np.ndarray:
    """Resample onto a 10 Hz grid, dropping samples off the expected cadence.

    Mirrors ``sample_10hz``: keep the first sample and any whose spacing from
    the previous kept sample is within [0.09, 0.11] s, then linearly
    interpolate (with extrapolation) onto ``new_time``. Dropping off-cadence
    samples removes duplicated or stalled records before interpolation.
    """
    original_time = np.asarray(original_time, dtype=float)
    original_value = np.asarray(original_value, dtype=float)

    spacing = np.diff(original_time)
    keep = np.concatenate(([True], (spacing > 0.09) & (spacing < 0.11)))
    kept_time = original_time[keep]
    kept_value = original_value[keep]

    return _interp_extrap(kept_time, kept_value, new_time)


@dataclass
class PreprocessedRun:
    """A run resampled to 10 Hz, with its assigned lane."""

    vin: int
    index: int
    direction: float
    assigned_lane: int
    timestamp: np.ndarray
    x_position: np.ndarray      # [ft], bumper-shifted, 10 Hz
    y_position: np.ndarray      # [ft]
    latitude: np.ndarray
    longitude: np.ndarray
    state_x: np.ndarray
    state_y: np.ndarray
    can_speed: np.ndarray
    control_active: np.ndarray  # bool
    #: Bounds of the *raw* run, before resampling. MATLAB's ``preproc_gps``
    #: overwrites ``timestamp`` with the 10 Hz grid but leaves
    #: ``starting_time``/``ending_time`` at the raw values, and the matching
    #: pass filters on those. The grid is floor/ceil'd outward by up to 0.1 s,
    #: so using it instead admits runs MATLAB excludes.
    starting_time: float = 0.0
    ending_time: float = 0.0


def preprocess_run(run: GpsRun, index: int, lane_map: Dict[int, int]) -> PreprocessedRun:
    """Resample one run to 10 Hz and attach its assigned lane.

    Mirrors the body of ``preproc_gps``. The 10 Hz grid runs from
    ``floor(t0*10)/10`` to ``ceil(tN*10)/10`` in 0.1 s steps. x is shifted from
    the antenna to the rear bumper before resampling: subtract
    ``direction * inVehShift``.
    """
    time = np.asarray(run.timestamp, dtype=float)
    grid = _tenth_second_grid(time[0], time[-1])

    x_shifted = np.asarray(run.x_position, dtype=float) - run.direction * IN_VEHICLE_SHIFT_FT
    control = resample_10hz(time, run.control_active.astype(np.float32), grid)

    return PreprocessedRun(
        vin=run.vin,
        index=index,
        direction=run.direction,
        assigned_lane=int(lane_map[run.vin]),
        timestamp=grid,
        x_position=resample_10hz(time, x_shifted, grid),
        y_position=resample_10hz(time, run.y_position, grid),
        latitude=resample_10hz(time, run.latitude, grid),
        longitude=resample_10hz(time, run.longitude, grid),
        state_x=resample_10hz(time, run.state_x, grid),
        state_y=resample_10hz(time, run.state_y, grid),
        can_speed=resample_10hz(time, run.can_speed, grid),
        # MATLAB: logical(round(sample_10hz(single(control_active)))).
        control_active=np.round(control).astype(bool),
        starting_time=float(run.starting_time),
        ending_time=float(run.ending_time),
    )


def load_lane_map(vins_csv) -> Dict[int, int]:
    """Map car number -> assigned lane from ``cars_vins.csv``.

    MATLAB indexes ``avVINData`` by the car number as a row (1-based), so car N
    is row N; this returns ``{veh_id: lane_num}`` which is equivalent and does
    not depend on row order.
    """
    import pandas as pd

    table = pd.read_csv(vins_csv, float_precision="round_trip")
    return {int(vehicle): int(lane)
            for vehicle, lane in zip(table["veh_id"], table["lane_num"])}


def connection_status(ping_csv, vins_csv) -> Dict[int, dict]:
    """Server-ping timestamps and connectivity per vehicle.

    Mirrors ``get_connection_status``: join the ping records to car numbers by
    VIN, and for each car return the ping times (``gpstime`` in ms -> s) and
    whether each ping carried a non-NaN ``acc_status`` (a live server record).
    """
    import pandas as pd

    # Correctly-rounded float parsing to match MATLAB readtable (see load_lane_map).
    pings = pd.read_csv(ping_csv, float_precision="round_trip")
    vins = pd.read_csv(vins_csv, float_precision="round_trip")
    vin_to_id = {str(vin): int(vehicle) for vin, vehicle in zip(vins["vin"], vins["veh_id"])}

    pings = pings.copy()
    pings["av_id"] = pings["vin"].astype(str).map(vin_to_id).fillna(0).astype(int)

    status: Dict[int, dict] = {}
    for vehicle_id in range(1, 104):
        rows = pings[pings["av_id"] == vehicle_id]
        status[vehicle_id] = {
            "timestamp": rows["gpstime"].to_numpy(dtype=float) / 1e3,
            "status": ~np.isnan(rows["acc_status"].to_numpy(dtype=float)),
        }
    return status


def is_server_connected(run_time: np.ndarray, ping_time: np.ndarray) -> np.ndarray:
    """Per-sample server-connected flag (the loop in the assembly section).

    A sample is connected if the most recent ping at or before it is no more
    than 2 s in the future of... actually no earlier than 2 s before it: the
    MATLAB condition is ``(lastPing - t) >= -2``. Returns 0/1 floats to match
    the released encoding.
    """
    run_time = np.asarray(run_time, dtype=float)
    ping_time = np.asarray(ping_time, dtype=float)
    connected = np.zeros(run_time.size)

    start = _find_last_before(ping_time, run_time[0])
    if start is None:
        return connected
    pings = ping_time[start:]
    if pings.size == 0 or pings[-1] <= run_time[0]:
        return connected

    for index, sample_time in enumerate(run_time):
        last = _find_last_at_or_before(pings, sample_time)
        chosen = pings[last] if last is not None else pings[0]
        connected[index] = float((chosen - sample_time) >= -2)
    return connected


def control_car_status(controller_engaged: np.ndarray, speed: np.ndarray,
                       x_position: np.ndarray, timestamp: np.ndarray):
    """Reconstructed control signal and its 30 s look-back.

    Port of ``get_control_car_status``. Corrects the raw controller-engaged
    signal so that a controller left on while the vehicle is briefly stopped is
    not read as disengaged, and forces disengagement during a manual stop.
    Returns ``(control_car, control_last30)`` as float 0/1 arrays.
    """
    controller_engaged = np.asarray(controller_engaged, dtype=bool)
    speed = np.asarray(speed, dtype=float)

    # Derive speed from position when CAN speed is entirely absent.
    if np.sum(speed) == 0:
        speed = _speed_from_position(x_position, timestamp)

    n = controller_engaged.size
    status = -np.ones(n)
    active = np.flatnonzero(controller_engaged)
    status[active] = 1
    status[(~controller_engaged) & (speed != 0)] = 0

    _resolve_active_while_stopped(status, controller_engaged, speed, active, n)
    _resolve_inactive_while_stopping(status, controller_engaged, speed, n)

    status[status == -1] = 0

    # control_last30: active at any point in the preceding 30 s (300 samples).
    last30 = np.array([float(np.any(status[max(0, i - 300):i + 1] > 0))
                       for i in range(n)])
    return status, last30


def assemble_run(run: PreprocessedRun, status: dict, median_xd: float = 0.0) -> dict:
    """Build one released GPS record from a preprocessed run.

    ``median_xd`` is the MOTION-matching bias (meters) subtracted from
    ``x_position``; pass 0 to get the record up to that per-run offset.
    """
    from .matround import round_decimals

    connected = is_server_connected(run.timestamp, status["timestamp"])

    control_car, control_last30 = control_car_status(
        run.control_active, run.can_speed, FT_TO_METER * run.x_position - median_xd,
        run.timestamp)

    record = {
        "av_id": run.vin,
        "assigned_lane": run.assigned_lane,
        "direction": run.direction,
        "timestamp": round_decimals(run.timestamp, 6),
        "latitude": round_decimals(run.latitude, 6),
        "longitude": round_decimals(run.longitude, 6),
        "x_position": round_decimals(FT_TO_METER * run.x_position - median_xd, 6),
        "y_position": round_decimals(FT_TO_METER * run.y_position, 6),
        # Stays boolean: MATLAB carries this as a logical, so jsonencode writes
        # true/false. Casting to float here encoded 0/1 instead, which is the
        # same information but not the same bytes - and at 3M samples it was
        # 9.8 MB of spurious difference against the released file.
        "controller_engaged": run.control_active.astype(bool),
        "speed": run.can_speed,
        "is_server_connected": connected,
        "first_timestamp": float(round_decimals(run.timestamp[0], 6)),
        "last_timestamp": float(round_decimals(run.timestamp[-1], 6)),
        "control_car": control_car,
        "control_last30": control_last30,
    }
    return _clip_to_testbed(record)


def assemble_day(day: int, data_dir, timings: Optional[dict] = None) -> List[dict]:
    """Full Python equivalent of assemble_data_GPS for one day.

    Ties the whole GPS stage together: parse runs, resample to 10 Hz, load
    server connectivity, compute the MOTION-matching bias, and assemble one
    released record per run (with x_position corrected by the bias). Returns the
    list of records in run order, ready for mvtpy.matjson.dumps.

    Parameters
    ----------
    day:       16, 17, or 18
    data_dir:  the workspace ``data`` folder (contains cars/ and i24motion/)
    timings:   optional dict; populated with per-phase wall-clock seconds

    Notes
    -----
    Byte-parity with MATLAB's output is limited by the CSV float-parse residual
    (see docs/PYTHON_PORT.md); the values are correct to ~1e-6. This is the
    heaviest stage - the matching pass streams the day's raw MOTION segments.
    """
    import time

    from . import gpsmatch, gpsruns

    data_dir = Path(data_dir)
    cars = data_dir / "cars"
    motion = data_dir / "i24motion" / f"2022-11-{day}"
    clock = {} if timings is None else timings

    def phase(name, fn):
        start = time.time()
        result = fn()
        clock[name] = time.time() - start
        return result

    runs = phase("parse_runs", lambda: gpsruns.parse_gps_data(cars / "cars_gps", day))
    lane_map = load_lane_map(cars / "cars_vins.csv")
    pre = phase("preprocess", lambda: [preprocess_run(run, index + 1, lane_map)
                                       for index, run in enumerate(runs)])
    status = phase("connection_status",
                   lambda: connection_status(cars / f"veh_ping_202211{day}.csv",
                                             cars / "cars_vins.csv"))
    bias = phase("matching_bias", lambda: gpsmatch.matching_bias(pre, motion, day))
    records = phase("assemble", lambda: [assemble_run(run, status[run.vin],
                                                      bias.get(run.index, 0.0))
                                         for run in pre])
    return records


# ---------------------------------------------------------------------------
# internals


def _tenth_second_grid(start: float, stop: float) -> np.ndarray:
    lo = np.floor(start * 10) / 10
    hi = np.ceil(stop * 10) / 10
    return _colon(lo, 0.1, hi)


def _colon(a: float, step: float, b: float) -> np.ndarray:
    """Reproduce MATLAB's ``a:step:b`` bit-for-bit.

    MATLAB's colon operator does NOT compute ``a + k*step``; it builds the
    vector from both ends to stay accurate at each end - the first half from
    ``a + k*step``, the second from ``b - (n-k)*step``. For a large base like a
    POSIX timestamp the two forms disagree in the last bit on ~20% of points,
    which shifts the interp1 query points and flips 6th-decimal roundings in the
    resampled GPS fields.

    When the number of steps is **even** there is an exact middle element, and
    MATLAB sets it to ``(a+b)/2`` rather than to either one-sided form. Missing
    that was worth 104 wrong grid points per day: neither a forward nor a
    backward split reproduces it, because the correct value is sometimes one and
    sometimes the other, and always the average.

    Verified bit-for-bit against `lo:0.1:hi` evaluated in MATLAB over the 772
    real 2022-11-18 run grids: 0 of 3,739,194 points differ.
    """
    n = int(round((b - a) / step))
    k = np.arange(n + 1)
    half = n // 2
    grid = np.empty(n + 1)
    if n % 2 == 0:
        grid[:half] = a + k[:half] * step
        grid[half] = (a + b) / 2
        grid[half + 1:] = b - (n - k[half + 1:]) * step
    else:
        grid[:half + 1] = a + k[:half + 1] * step
        grid[half + 1:] = b - (n - k[half + 1:]) * step
    return grid


def _interp_extrap(xp: np.ndarray, fp: np.ndarray, x: np.ndarray) -> np.ndarray:
    """Linear interpolation with linear extrapolation, bit-identical to interp1.

    MATLAB's ``interp1`` evaluates a segment as the weighted blend
    ``A(i)*(1-w) + A(i+1)*w`` with ``w = (x-B(i))/(B(i+1)-B(i))``, NOT the
    algebraically equivalent ``A(i) + slope*(x-B(i))``. The two disagree in the
    last bit (~1e-14), which is invisible at the 4-decimal rounding of the MOTION
    data but flips 6th-decimal roundings in the GPS data.

    On a **flat segment** (``A(i) == A(i+1)``) the blend is not exact: the two
    products round independently and their sum lands a ULP off the common value,
    whereas MATLAB returns the value itself. Speed data is full of flat runs, so
    this was the single remaining source of 1-ULP differences in the resampled
    fields. Verified against ``interp1(B,A,newT,'linear','extrap')`` on a real
    9,449-point run: the plain blend differs on 88 values, all of them flat
    segments; with this guard, 0 differ.
    """
    xp = np.asarray(xp, dtype=float)
    fp = np.asarray(fp, dtype=float)
    x = np.asarray(x, dtype=float)
    index = np.clip(np.searchsorted(xp, x, side="right") - 1, 0, len(xp) - 2)
    left, right = fp[index], fp[index + 1]
    weight = (x - xp[index]) / (xp[index + 1] - xp[index])
    return np.where(left == right, left, left * (1 - weight) + right * weight)


def _speed_from_position(x_position: np.ndarray, timestamp: np.ndarray) -> np.ndarray:
    x = np.asarray(x_position, dtype=float)
    t = np.asarray(timestamp, dtype=float)
    raw = (x[1:] - x[:-1]) / (t[1:] - t[:-1])
    raw = np.round(np.abs(np.concatenate((raw, raw[-1:]))) * 2) / 2

    window = 10
    padded = np.concatenate((np.zeros(window // 2), raw, np.zeros(window // 2)))
    filtered = _movmedian(padded, window)
    return filtered[window // 2:len(filtered) - window // 2]


def _movmedian(values: np.ndarray, window: int) -> np.ndarray:
    """Centered moving median matching MATLAB's movmedian window convention.

    For even window k, MATLAB uses the k+1 samples centered on each point:
    floor(k/2) before and floor(k/2) after (edges shrink the window).
    """
    n = values.size
    half = window // 2
    out = np.empty(n)
    for i in range(n):
        lo = max(0, i - half)
        hi = min(n, i + half + 1)
        out[i] = np.median(values[lo:hi])
    return out


def _resolve_active_while_stopped(status, controller_engaged, speed, active, n):
    """Reinstate control that was dropped only because the vehicle stopped."""
    stopped_active = np.flatnonzero(controller_engaged & (speed == 0))
    if stopped_active.size == 0:
        return
    # Edges: a stopped-active index whose successor is not the next index.
    successors = np.append(stopped_active[1:], n)
    edges = stopped_active[np.diff(np.append(stopped_active, n)) > 1] \
        if stopped_active.size else stopped_active
    # Recompute edges exactly as MATLAB: diff([indActiveStopped; n]) > 1.
    extended = np.append(stopped_active, n)
    edges = stopped_active[np.diff(extended) > 1]

    for index in edges:
        after = active[active > index]
        end = (after[0] - 1) if after.size else (n - 1)
        if end == index:
            continue
        max_steps = 15
        moving = np.sum(speed[index:end + 1] > 0)
        if moving < min(max_steps, end - index + 2):
            status[index:end + 1] = 1
        else:
            status[index:end + 1] = 0


def _resolve_inactive_while_stopping(status, controller_engaged, speed, n):
    """Force disengagement while the vehicle comes to a manual stop."""
    coming_to_stop = np.flatnonzero(
        (~controller_engaged[:-1]) & (speed[:-1] != 0) & (speed[1:] == 0))
    if coming_to_stop.size == 0:
        return
    resume_points = np.flatnonzero(controller_engaged | (speed != 0))
    for index in coming_to_stop:
        after = resume_points[resume_points > index]
        end = (after[0] - 1) if after.size else (n - 1)
        status[index:end + 1] = 0


def _clip_to_testbed(record: dict) -> dict:
    """Keep only samples inside the testbed x limits (in meters)."""
    lo, hi = (limit * FT_TO_METER for limit in TESTBED_X_LIMITS_FT)
    x = record["x_position"]
    inside = (x > lo) & (x < hi)
    if not np.any(inside):
        return record

    clipped = dict(record)
    for key, value in record.items():
        if isinstance(value, np.ndarray) and value.size > 1:
            clipped[key] = value[inside]
    clipped["first_timestamp"] = float(clipped["timestamp"][0])
    clipped["last_timestamp"] = float(clipped["timestamp"][-1])
    return clipped


def _find_last_before(times: np.ndarray, value: float) -> Optional[int]:
    found = np.flatnonzero(times < value)
    return int(found[-1]) if found.size else None


def _find_last_at_or_before(times: np.ndarray, value: float) -> Optional[int]:
    found = np.flatnonzero(times <= value)
    return int(found[-1]) if found.size else None
