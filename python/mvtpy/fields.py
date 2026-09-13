"""Macroscopic traffic-state fields from processed MOTION trajectories.

Port of ``generate_macroscopic_fields.m``. Converts one day of westbound
trajectories into fields on a regular time-space grid (5 s x 50 m, 06:00-10:00
local, x in [0, 6500] m): density ``Rho``, flow ``Q``, fuel-rate density ``F``,
and the derived ``U = Q/Rho`` (speed), ``Phi = F/Rho`` (bulk fuel rate), and
``Psi = F/Q`` (bulk fuel consumption). These underlie Figure 3 and SM5.

Each trajectory point contributes ``dt * (1, v, f)`` to every grid cell within a
box window (+/- ht in time, +/- hx in space; box kernel, so weight 1 inside).
The accumulation is done with a 2D difference array (add to a rectangle in O(1),
integrate at the end) rather than MATLAB's per-point rectangle add. That is
mathematically identical but sums in a different order, so the fields agree with
MATLAB's to floating-point rounding (~1e-9 relative), not bit-for-bit - which is
the right target for a computed physical field compared with tolerance.
"""

from __future__ import annotations

from dataclasses import dataclass, field as _dc_field
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable
from zoneinfo import ZoneInfo

import numpy as np

__all__ = ["FieldOptions", "macroscopic_fields", "macroscopic_fields_from_dir"]


@dataclass(frozen=True)
class FieldOptions:
    """Constants from the head of generate_macroscopic_fields.m."""

    hx: float = 100.0            # [m] spatial window half-width
    v_char: float = 20.0         # [m/s] characteristic speed -> temporal window
    direction: int = -1          # westbound
    lane: int = 0                # all lanes
    t_start: str = "06:00:00"
    t_end: str = "10:00:00"
    x_min: float = 0.0
    x_max: float = 6500.0
    t_res: float = 5.0           # [s]
    x_res: float = 50.0          # [m]
    skip_t_compute: int = 5      # sub-sampling of trajectory points

    @property
    def ht(self) -> float:
        """[s] temporal window half-width, hx / v_char."""
        return self.hx / self.v_char

    @property
    def factor(self) -> float:
        """Box-kernel normalization, 1 / (4 * ht * hx)."""
        return 1.0 / (4 * self.ht * self.hx)


def macroscopic_fields(records: Iterable[dict], day: int,
                       options: FieldOptions = FieldOptions()) -> dict:
    """Accumulate the macroscopic fields from processed trajectory records.

    Returns a dict with ``t``, ``x`` (grid axes), ``direction``, ``lane``, and
    ``field`` -> {Rho, Q, F, U, Phi, Psi} arrays of shape (len(t), len(x)).
    """
    zone = ZoneInfo("America/Chicago")
    t0 = datetime.fromisoformat(f"2022-11-{day:02d} {options.t_start}").replace(tzinfo=zone).timestamp()
    t1 = datetime.fromisoformat(f"2022-11-{day:02d} {options.t_end}").replace(tzinfo=zone).timestamp()

    t = _colon(t0, options.t_res, t1)
    x = _colon(options.x_min, options.x_res, options.x_max)
    ht, hx = options.ht, options.hx

    # Difference arrays, one per accumulated field (Rho, Q, F). One extra row and
    # column so a rectangle's far corner can be marked without bounds trouble.
    shape = (t.size + 1, x.size + 1)
    diff_rho = np.zeros(shape)
    diff_q = np.zeros(shape)
    diff_f = np.zeros(shape)

    for veh in records:
        if veh.get("direction", options.direction) * options.direction <= 0:
            continue
        if not (options.lane == 0 or veh["lane_number"] == options.lane):
            continue

        traj_t_full = np.asarray(veh["timestamp"], dtype=float)
        traj_x_full = np.asarray(veh["x_position_meters"], dtype=float)
        if traj_t_full.size < 2:
            continue

        # Central-difference speed on the full trajectory, then sub-sample.
        forward_x = np.concatenate((traj_x_full[1:], traj_x_full[-1:]))
        backward_x = np.concatenate((traj_x_full[:1], traj_x_full[:-1]))
        forward_t = np.concatenate((traj_t_full[1:], traj_t_full[-1:]))
        backward_t = np.concatenate((traj_t_full[:1], traj_t_full[:-1]))
        speed = (forward_x - backward_x) / (forward_t - backward_t) * options.direction

        step = options.skip_t_compute
        traj_t = traj_t_full[::step]
        traj_x = traj_x_full[::step]
        traj_v = speed[::step]
        traj_f = np.asarray(veh["fuel_rate_grams_per_second"], dtype=float)[::step]

        # Trapezoidal dt: average of the gaps on each side of a point.
        gaps = np.diff(traj_t)
        traj_dt = (np.concatenate((gaps[:1], gaps))
                   + np.concatenate((gaps, gaps[-1:]))) / 2

        # Box bounds per point (inclusive |grid - point| <= h).
        t_lo = np.searchsorted(t, traj_t - ht, side="left")
        t_hi = np.searchsorted(t, traj_t + ht, side="right")   # exclusive end
        x_lo = np.searchsorted(x, traj_x - hx, side="left")
        x_hi = np.searchsorted(x, traj_x + hx, side="right")

        _scatter_rectangles(diff_rho, t_lo, t_hi, x_lo, x_hi, traj_dt)
        _scatter_rectangles(diff_q, t_lo, t_hi, x_lo, x_hi, traj_dt * traj_v)
        _scatter_rectangles(diff_f, t_lo, t_hi, x_lo, x_hi, traj_dt * traj_f)

    rho = _integrate(diff_rho, t.size, x.size) * options.factor
    q = _integrate(diff_q, t.size, x.size) * options.factor
    f = _integrate(diff_f, t.size, x.size) * options.factor

    with np.errstate(invalid="ignore", divide="ignore"):
        u = q / rho
        phi = f / rho
        psi = f / q
    u[rho < 1e-3] = np.nan
    phi[rho < 1e-3] = np.nan
    psi[q < 1e-2] = np.nan

    return {
        "t": t, "x": x, "direction": options.direction, "lane": options.lane,
        "field": {"Rho": rho, "Q": q, "F": f, "U": u, "Phi": phi, "Psi": psi},
    }


def macroscopic_fields_from_dir(slim_dir, day: int,
                                options: FieldOptions = FieldOptions()) -> dict:
    """Accumulate fields from a day's processed JSON files."""
    from .rawio import iter_trajectories

    slim_dir = Path(slim_dir)
    files = sorted(slim_dir.glob("I-24MOTION_*.json"))

    def all_records():
        for path in files:
            yield from iter_trajectories(path)

    return macroscopic_fields(all_records(), day, options)


# ---------------------------------------------------------------------------


def _scatter_rectangles(diff, t_lo, t_hi, x_lo, x_hi, values):
    """Add each value to a rectangle of the difference array (four corners)."""
    np.add.at(diff, (t_lo, x_lo), values)
    np.add.at(diff, (t_hi, x_lo), -values)
    np.add.at(diff, (t_lo, x_hi), -values)
    np.add.at(diff, (t_hi, x_hi), values)


def _integrate(diff, n_t, n_x):
    """Two cumulative sums turn the difference array into the accumulated field."""
    return np.cumsum(np.cumsum(diff, axis=0), axis=1)[:n_t, :n_x]


def _colon(a: float, step: float, b: float) -> np.ndarray:
    """MATLAB a:step:b (both-ends construction; see gpsassemble._colon)."""
    n = int(round((b - a) / step))
    k = np.arange(n + 1)
    half = n // 2
    grid = np.empty(n + 1)
    grid[:half + 1] = a + k[:half + 1] * step
    grid[half + 1:] = b - (n - k[half + 1:]) * step
    return grid
