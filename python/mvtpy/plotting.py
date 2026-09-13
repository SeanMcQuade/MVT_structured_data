"""Matplotlib renderings of the MVT figures.

These are *not* pixel-for-pixel reproductions of the MATLAB figures. They use the
same colors, colormap (MATLAB's ``parula``), color-scale limits, and overall
layout so the output reads as the same figure, while leaving exact geometry to
matplotlib. They consume the data structures the ported stages produce
(``mvtpy.fields`` and ``mvtpy.avanalysis``).

Two figures are provided:

* :func:`plot_field` - a macroscopic field heatmap with the control-vehicle GPS
  overlay (Figure 3 / SM5 family), matching ``plot_macroscopic_fields.m``.
* :func:`plot_av_fuel` - the effective fuel consumption vs. distance to the
  nearest engaged AV, ahead and behind, per day (the fuel core of Figure 2),
  matching ``plot_AV_analysis.m``.
"""

from __future__ import annotations

from datetime import datetime, timezone
from typing import Dict, Optional, Sequence
from zoneinfo import ZoneInfo

import numpy as np

__all__ = ["PARULA", "parula_cmap", "FIELD_COLOR_LIMITS", "DAY_COLORS",
           "plot_field", "plot_av_fuel"]

#: MATLAB parula colormap anchor points (0..1), interpolated to a smooth map.
PARULA = np.array([
    [0.2422, 0.1504, 0.6603],
    [0.2803, 0.3220, 0.9570],
    [0.1786, 0.5289, 0.9682],
    [0.0689, 0.6948, 0.8394],
    [0.2161, 0.7843, 0.5923],
    [0.4668, 0.7772, 0.3798],
    [0.6720, 0.7793, 0.2227],
    [0.9970, 0.7659, 0.2199],
    [0.9769, 0.9839, 0.0805],
])

#: Per-field color-scale upper limits (caxis [0, limit]) from
#: plot_macroscopic_fields.m, and the display unit.
FIELD_COLOR_LIMITS = {
    "Rho": (220.0, "veh/km"),
    "Q": (8000.0, "veh/h"),
    "F": (0.25, "g/s"),
    "U": (35.0, "m/s"),
    "Phi": (4.5, "g/s"),
    "Psi": (0.22, "g/m"),
}

#: Field display factors (field.factor in the MATLAB struct).
FIELD_FACTORS = {"Rho": 1000.0, "Q": 3600.0, "F": 1.0, "U": 1.0, "Phi": 1.0, "Psi": 1.0}

#: Per-day line colors (plotCol in plot_AV_analysis.m): blue, red, green.
DAY_COLORS = {16: (0.0, 0.0, 1.0), 17: (1.0, 0.0, 0.0), 18: (0.0, 0.6, 0.0)}

_CENTRAL = ZoneInfo("America/Chicago")


def parula_cmap():
    """Return MATLAB's parula as a matplotlib colormap."""
    from matplotlib.colors import LinearSegmentedColormap

    return LinearSegmentedColormap.from_list("parula", PARULA)


def plot_field(fields: dict, name: str, gps_records: Optional[Sequence[dict]] = None,
               ax=None, zoom_military=(610, 950)):
    """Heatmap of a macroscopic field, matching plot_macroscopic_fields.m.

    Parameters
    ----------
    fields:       output of ``mvtpy.fields.macroscopic_fields``
    name:         one of Rho, Q, F, U, Phi, Psi
    gps_records:  control-vehicle GPS records to overlay (engaged red, disengaged
                  white), or None to skip the overlay
    ax:           axis to draw on (a new figure is made if None)
    zoom_military: time window (military) to zoom the x axis to

    Returns the matplotlib Axes.
    """
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates

    t = np.asarray(fields["t"], dtype=float)
    x = np.asarray(fields["x"], dtype=float)
    value = np.asarray(fields["field"][name], dtype=float) * FIELD_FACTORS[name]
    limit, unit = FIELD_COLOR_LIMITS[name]

    if ax is None:
        _, ax = plt.subplots(figsize=(15, 5))
    times = [datetime.fromtimestamp(ts, _CENTRAL) for ts in t]

    # value is (t, x); show x (km) on the vertical axis, time on the horizontal.
    mesh = ax.pcolormesh(mdates.date2num(times), x / 1000.0, value.T,
                         cmap=parula_cmap(), vmin=0.0, vmax=limit, shading="nearest")
    colorbar = ax.figure.colorbar(mesh, ax=ax)
    colorbar.set_label(unit)

    if gps_records:
        _overlay_av_trajectories(ax, gps_records)
        engaged = ax.plot([], [], "r.", markersize=12, label="AV engaged")[0]
        disengaged = ax.plot([], [], "o", markerfacecolor="none",
                             markeredgecolor="0.4", markersize=6, label="AV disengaged")[0]
        ax.legend(handles=[engaged, disengaged], fontsize=11, loc="upper right")

    ax.invert_yaxis()   # position 0 at top, as in MATLAB's imagesc
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M:%S", tz=_CENTRAL))
    ax.set_xlabel("time (America/Chicago)")
    ax.set_ylabel("position / km")
    ax.set_title(_field_title(times[0], name, fields))

    _apply_time_zoom(ax, times[0], zoom_military)
    return ax


def plot_av_fuel(bin_stats_by_day: Dict[int, dict], statistic: str = "effective",
                 masking_dist: float = 30.0, axes=None):
    """Fuel consumption vs. distance to the nearest engaged AV, per day.

    Reproduces the fuel core of Figure 2: two panels (ahead of the AV, behind the
    AV), one colored line per day (blue/red/green), with the near-AV band shown
    faded, as in plot_AV_analysis.m.

    Parameters
    ----------
    bin_stats_by_day: {day: output of mvtpy.avanalysis.bin_samples}
    statistic:        'effective', 'mean', or 'median'
    masking_dist:     [m] band around the AV drawn faded
    axes:             a pair of axes (ahead, behind), created if None
    """
    import matplotlib.pyplot as plt

    if axes is None:
        _, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    ahead_ax, behind_ax = axes

    for day, stats in sorted(bin_stats_by_day.items()):
        centers = np.asarray(stats["centers"], dtype=float)
        values = np.asarray(stats[statistic], dtype=float)
        color = DAY_COLORS.get(day, (0.3, 0.3, 0.3))

        ahead = centers < 0
        behind = centers > 0
        _plot_side(ahead_ax, -centers[ahead], values[ahead], color, masking_dist, day)
        _plot_side(behind_ax, centers[behind], values[behind], color, masking_dist, None)

    ahead_ax.set_title("ahead of AV")
    behind_ax.set_title("behind AV")
    for ax in axes:
        ax.set_xlabel("distance to nearest engaged AV / m")
    ahead_ax.set_ylabel(f"{statistic} fuel consumption / (g/m)")
    ahead_ax.legend(title="day", fontsize=10)
    return axes


# ---------------------------------------------------------------------------


def _plot_side(ax, distance, values, color, masking_dist, day):
    """Draw one day's curve; the near-AV band (< masking_dist) is faded."""
    order = np.argsort(distance)
    distance, values = distance[order], values[order]
    near = distance < masking_dist
    label = f"Nov {day}" if day is not None else None
    ax.plot(distance[~near], values[~near], color=color, linewidth=1.5, label=label)
    ax.plot(distance[near], values[near], color=color, linewidth=1.5, alpha=0.3)


def _overlay_av_trajectories(ax, gps_records):
    import matplotlib.dates as mdates

    for rec in gps_records:
        if rec.get("direction", -1) >= 0:
            continue
        control = np.asarray(rec.get("control_car", rec.get("controller_engaged")))
        t = np.asarray(rec["timestamp"], dtype=float)
        x = np.asarray(rec["x_position"], dtype=float)
        times = mdates.date2num([datetime.fromtimestamp(ts, _CENTRAL) for ts in t])
        disengaged = control == 0
        ax.plot(times[disengaged], x[disengaged] / 1000.0, "w.", markersize=1)
        ax.plot(times[~disengaged], x[~disengaged] / 1000.0, "r.", markersize=1)


def _field_title(first_time, name, fields):
    names = {"Rho": "vehicle density", "Q": "flow rate", "F": "fuel rate density",
             "U": "bulk velocity", "Phi": "bulk fuel rate", "Psi": "bulk fuel consumption"}
    weekday = first_time.strftime("%A")
    date = first_time.strftime("%d-%b-%Y")
    return f"Westbound (all lanes) on {weekday} {date}: {names.get(name, name)}"


def _apply_time_zoom(ax, first_time, zoom_military):
    """Zoom the x axis to the [lo, hi] military-time window within the morning."""
    import matplotlib.dates as mdates

    lo, hi = zoom_military
    base = first_time.replace(hour=0, minute=0, second=0, microsecond=0)

    def at(military):
        return base.replace(hour=military // 100, minute=military % 100)

    ax.set_xlim(mdates.date2num(at(lo)), mdates.date2num(at(hi)))
