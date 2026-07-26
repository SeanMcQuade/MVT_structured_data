"""Time-space plots of individual trajectories (stage 5b).

Ports ``Scripts/plot_microscopic_trajectories.m``: every trajectory is drawn as
a ribbon whose vertical extent is the vehicle's occupied space, colored by
speed, over a time-space plot of the westbound testbed. Three figures come out
of one render pass - the full field, the same field with the zoom window marked,
and the zoom itself.

Memory
------
The MATLAB original called ``patch()`` once per trajectory, which builds ~10^5
persistent graphics objects for a full day - each with kilobyte-scale fixed
overhead independent of vertex count - and exhausted RAM. It was fixed by
batching each file into a single ``patch`` with a NaN-padded ``Faces`` matrix.

This port applies the same fix in matplotlib terms, and the choice of primitive
is the whole point:

* One artist **per file**, not per trajectory (24 artists for a day, not
  ~370,000).
* That artist is a ``TriMesh`` (via :func:`~matplotlib.axes.Axes.tripcolor` with
  ``shading="gouraud"``), which stores the geometry as three flat arrays and
  creates no per-triangle Python object. A ``PolyCollection`` would be the more
  obvious translation of ``patch``, but ``set_verts`` materializes one ``Path``
  per polygon, reintroducing the object-count problem this fix exists to avoid.
  Gouraud shading also reproduces MATLAB's ``'FaceColor','interp'`` directly.
* Trajectories are converted to geometry as they stream off disk and the decoded
  record is dropped immediately, so a file's decoded JSON is never held whole -
  peak memory is the compact float arrays (tens of MB), not the multi-GB decode
  the MATLAB version needs its ``*_reduced.mat`` caches to avoid.
* The zoom is produced by changing the axes limits on the *same* artists. The
  MATLAB original walked every patch object and deleted the ones outside the
  window, a second O(#trajectories) pass that spiked RAM; axes clipping does the
  same job for free.

Figures are close in color and layout to the MATLAB output, not pixel-exact,
matching the convention of :mod:`mvtpy.plotting`.
"""

from __future__ import annotations

import zlib
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import List, Optional, Sequence, Tuple
from zoneinfo import ZoneInfo

import numpy as np

from .rawio import iter_trajectories

__all__ = ["MicroOptions", "speed_cmap", "trajectory_mesh", "file_mesh",
           "add_mesh", "setup_axes", "apply_limits", "render_day",
           "save_figures", "plot_microscopic_trajectories"]

_CENTRAL = ZoneInfo("America/Chicago")

#: [m/ft] conversion factor from feet to meter (vehicle length arrives in feet).
FT2M = 0.3048


@dataclass(frozen=True)
class MicroOptions:
    """Tunables, mirroring the parameter block of the MATLAB script.

    The MATLAB defaults are reproduced exactly, except ``skip_t_plot``, which is
    documented there as "use 5 to 50" and left at the value the released figures
    were made with.
    """

    direction: int = -1                                   # -1 westbound, 1 eastbound
    lane: int = 0                                         # 0 = all lanes
    time_window: Tuple[str, str] = ("06:00:00", "10:00:00")
    x_window_m: Tuple[float, float] = (0.0, 6500.0)
    skip_t_plot: int = 50                                 # sub-sample for plotting
    x_field: str = "x_position_meters"
    speed_limit: float = 35.0                             # [m/s] color-scale top
    fig_size_in: Tuple[float, float] = (25.0, 8.0)        # MATLAB fig_res/100
    dpi: int = 200                                        # -> ~5000x1600 px
    font_size: int = 20
    n_xticks: int = 16
    trim_minutes: float = 10.0                            # crop each end of the window
    zoom_time_military: Tuple[int, int] = (612, 619)
    zoom_x_m: Tuple[float, float] = (845.0, 1495.0)
    #: Axes box as a fraction of the figure, from the MATLAB `set(gca,'Position')`.
    #: Nearly full-bleed; the default matplotlib box is far narrower and crowds
    #: the 17 time labels into each other.
    axes_position: Tuple[float, float, float, float] = (0.023, 0.093, 0.90, 0.866)
    colorbar_position: Tuple[float, float, float, float] = (0.935, 0.093, 0.012, 0.866)
    #: Width of the zoom figure relative to the main one. MATLAB computes a
    #: factor from the window's aspect ratio, but under ``flag_use_int32_vars``
    #: that expression mixes scaled and unscaled units, comes out ~350, and the
    #: window is silently clamped to the screen - so the released zoom is simply
    #: main-figure width. None reproduces that; a float applies the factor.
    zoom_width_factor: Optional[float] = None

    def window_bounds(self, day: int) -> Tuple[float, float]:
        """[s] posix bounds of the plotted time window for a November 2022 day."""
        base = datetime(2022, 11, day, tzinfo=_CENTRAL)

        def at(hms: str) -> float:
            hour, minute, second = (int(part) for part in hms.split(":"))
            return base.replace(hour=hour, minute=minute, second=second).timestamp()

        return at(self.time_window[0]), at(self.time_window[1])


def speed_cmap():
    """Red (stopped) to green (free-flow), as in the MATLAB ``cMap``."""
    from matplotlib.colors import LinearSegmentedColormap

    return LinearSegmentedColormap.from_list("mvt_speed", [(1.0, 0.0, 0.0),
                                                           (0.0, 1.0, 0.0)])


# --- geometry -------------------------------------------------------------


def trajectory_mesh(record: dict, opts: MicroOptions, t0: float):
    """Triangulated ribbon for one trajectory.

    Returns ``(points, triangles, speeds)`` with ``points`` in (seconds since
    ``t0``, km), ``triangles`` indexing into ``points``, and ``speeds`` per
    point in m/s - or None if the trajectory is too short to draw.

    The ribbon spans the vehicle's occupied space: the sampled positions on one
    edge and those positions offset by the vehicle length on the other, matching
    the MATLAB polygon.
    """
    t = np.asarray(record["timestamp"], dtype=float)
    x = np.asarray(record[opts.x_field], dtype=float)
    if t.size < 2 or t.size != x.size:
        return None

    # Centered difference with one-sided ends, expression-for-expression from the
    # MATLAB source (which recomputes speed here rather than reading the released
    # speed_meters_per_second field).
    t_next, t_prev = np.concatenate((t[1:], t[-1:])), np.concatenate((t[:1], t[:-1]))
    x_next, x_prev = np.concatenate((x[1:], x[-1:])), np.concatenate((x[:1], x[:-1]))
    with np.errstate(invalid="ignore", divide="ignore"):
        v = (x_next - x_prev) / (t_next - t_prev) * opts.direction
    v[~np.isfinite(v)] = 0.0

    index = _subsample_indices(t.size, opts.skip_t_plot)
    ts, xs, vs = t[index], x[index], v[index]
    n = index.size
    length_m = float(record.get("length", 0.0)) * FT2M

    points = np.empty((2 * n, 2), dtype=np.float64)
    points[:, 0] = np.concatenate((ts, ts)) - t0
    points[:, 1] = np.concatenate((xs, xs + length_m * opts.direction)) / 1000.0
    speeds = np.concatenate((vs, vs))

    # Two triangles per segment, stitching the near edge (0..n-1) to the far
    # edge (n..2n-1).
    near = np.arange(n - 1)
    far = near + n
    triangles = np.concatenate((np.stack((near, near + 1, far), axis=1),
                                np.stack((near + 1, far + 1, far), axis=1)))
    return points, triangles, speeds


def _subsample_indices(n: int, skip: int) -> np.ndarray:
    """MATLAB's ``round(linspace(1,n,ns))``, 0-based.

    ``np.floor(v + 0.5)`` is MATLAB's round-half-away-from-zero for the
    non-negative values here, where numpy's banker's rounding would differ.
    """
    count = max(int(np.ceil(n / skip)), 2)
    return np.floor(np.linspace(0, n - 1, count) + 0.5).astype(np.intp)


def _keep(record: dict, opts: MicroOptions) -> bool:
    """Direction and lane selection, as in the MATLAB index set.

    Slim segments carry no ``direction`` field (they are westbound-only), in
    which case direction filtering is skipped and only the lane applies.
    """
    if opts.lane and record.get("lane_number") != opts.lane:
        return False
    direction = record.get("direction")
    return direction is None or direction * opts.direction > 0


def file_mesh(path: "str | Path", opts: MicroOptions, t0: float):
    """Triangulate every selected trajectory in one segment file.

    Streams the file, converts each record to geometry, and drops the decoded
    record immediately, so peak memory is the compact arrays rather than the
    whole decoded segment.

    Returns ``(points, triangles, speeds)`` for the file, or None if nothing was
    selected.
    """
    pieces: List[Tuple[float, tuple]] = []
    for record in iter_trajectories(path):
        if not _keep(record, opts):
            continue
        mesh = trajectory_mesh(record, opts, t0)
        if mesh is not None:
            pieces.append((float(record.get("length", 0.0)), mesh))
    if not pieces:
        return None

    # Longest vehicles first, so shorter ones end up drawn on top; triangle
    # order is draw order within the mesh. (MATLAB sorts the index set the same
    # way before building patches.)
    pieces.sort(key=lambda item: -item[0])
    return _concatenate([mesh for _, mesh in pieces])


def _concatenate(meshes: Sequence[tuple]):
    points = np.concatenate([mesh[0] for mesh in meshes])
    speeds = np.concatenate([mesh[2] for mesh in meshes])
    triangles = np.empty((sum(mesh[1].shape[0] for mesh in meshes), 3), dtype=np.int32)
    row = 0
    offset = 0
    for mesh_points, mesh_triangles, _ in meshes:
        end = row + mesh_triangles.shape[0]
        triangles[row:end] = mesh_triangles + offset
        row, offset = end, offset + mesh_points.shape[0]
    return points, triangles, speeds


# --- rendering ------------------------------------------------------------


def add_mesh(ax, mesh, opts: MicroOptions):
    """Add one file's geometry to the axes as a single ``TriMesh`` artist."""
    from matplotlib.tri import Triangulation

    points, triangles, speeds = mesh
    triangulation = Triangulation(points[:, 0], points[:, 1], triangles)
    return ax.tripcolor(triangulation, speeds, shading="gouraud",
                        cmap=speed_cmap(), vmin=0.0, vmax=opts.speed_limit)


def setup_axes(ax, day: int, opts: MicroOptions) -> None:
    """Title, ticks, limits and colors of the time-space axes."""
    t_start, t_end = opts.window_bounds(day)
    span = t_end - t_start

    ax.set_facecolor("black")
    ax.tick_params(labelsize=opts.font_size)

    ticks = np.linspace(0.0, span, opts.n_xticks + 1)
    labels = [datetime.fromtimestamp(t_start + offset, _CENTRAL).strftime("%H:%M:%S")
              for offset in ticks]
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels)
    ax.set_yticks(np.arange(0, 8))

    ax.set_xlabel("time", fontsize=opts.font_size)
    ax.set_ylabel("position / km", fontsize=opts.font_size)
    ax.set_title(_title(t_start, opts), fontsize=opts.font_size)
    apply_limits(ax, day, opts)


def apply_limits(ax, day: int, opts: MicroOptions) -> None:
    """Set the plotted window. Idempotent, so it can re-assert the limits after
    ``tripcolor`` autoscales them."""
    t_start, t_end = opts.window_bounds(day)
    trim = opts.trim_minutes * 60.0
    ax.set_xlim(trim, (t_end - t_start) - trim)
    # High position first: MATLAB's `axis ij` puts position 0 at the top. Setting
    # the limits reversed is idempotent where invert_yaxis() would toggle.
    low, high = np.asarray(opts.x_window_m) / 1000.0
    ax.set_ylim(high, low)


def _title(t_start: float, opts: MicroOptions) -> str:
    start = datetime.fromtimestamp(t_start, _CENTRAL)
    direction_name = "Westbound" if opts.direction < 0 else "Eastbound"
    lane_name = f"lane {opts.lane}" if opts.lane else "all lanes"
    date = f"{start.day}-{start.strftime('%b-%Y')}"
    return (f"{direction_name} ({lane_name}) on {start.strftime('%A')} {date} "
            f"in UTC{start.strftime('%z')}: trajectories")


def render_day(slim_dir: "str | Path", day: int, opts: Optional[MicroOptions] = None,
               cache_dir: "str | Path | None" = None, log=None):
    """Build the full time-space figure for one day.

    One ``TriMesh`` is added per segment file, keeping the artist count at ~24
    for a day. Returns ``(figure, axes)``; the caller renders and saves.

    ``cache_dir`` stores each file's triangulation as an ``.npz`` so later runs
    skip the JSON decode, which dominates the runtime. A cache is reused only
    when it is newer than its source segment *and* was built with the same
    geometry options (see :func:`_cache_key`).
    """
    import matplotlib.pyplot as plt

    opts = opts or MicroOptions()
    log = log or (lambda *_: None)
    t_start, _ = opts.window_bounds(day)
    files = sorted(Path(slim_dir).glob("I-24MOTION_*.json"))
    if not files:
        raise FileNotFoundError(f"no I-24MOTION_*.json segments in {slim_dir}")

    figure = plt.figure(figsize=opts.fig_size_in)
    figure.patch.set_facecolor("white")
    ax = figure.add_axes(opts.axes_position)
    setup_axes(ax, day, opts)

    mappable = None
    for path in files:
        mesh = _cached_file_mesh(path, opts, t_start, cache_dir, log)
        if mesh is None:
            continue
        log(f"    adding {path.name}: {mesh[1].shape[0]} triangles")
        mappable = add_mesh(ax, mesh, opts)

    if mappable is not None:
        # Its own axes, so the main axes keeps the full-bleed geometry above
        # (colorbar(ax=...) would shrink it instead).
        bar = figure.colorbar(mappable, cax=figure.add_axes(opts.colorbar_position))
        bar.set_label("vehicle speed / m/s", fontsize=opts.font_size)
        bar.ax.tick_params(labelsize=opts.font_size)

    apply_limits(ax, day, opts)     # tripcolor autoscales; restore the window
    return figure, ax


def _cache_key(opts: MicroOptions, t0: float) -> np.ndarray:
    """The options that change the cached geometry.

    Stored alongside the arrays and checked on load, so that changing e.g.
    ``skip_t_plot`` rebuilds instead of silently reusing a mesh built at a
    different sub-sampling.
    """
    return np.array([opts.skip_t_plot, opts.direction, opts.lane, t0,
                     zlib.crc32(opts.x_field.encode())], dtype=float)


def _cached_file_mesh(path: Path, opts: MicroOptions, t0: float,
                      cache_dir: "str | Path | None", log):
    if cache_dir is None:
        log(f"    reading {path.name}")
        return file_mesh(path, opts, t0)

    cache_dir = Path(cache_dir)
    cache = cache_dir / f"{path.stem}_micro.npz"
    key = _cache_key(opts, t0)
    if cache.is_file() and cache.stat().st_mtime >= path.stat().st_mtime:
        with np.load(cache) as data:
            if "key" in data and np.array_equal(data["key"], key):
                log(f"    reusing {cache.name}")
                return data["points"], data["triangles"], data["speeds"]
        log(f"    rebuilding {cache.name}: built with different options")

    log(f"    reading {path.name}")
    mesh = file_mesh(path, opts, t0)
    if mesh is not None:
        cache_dir.mkdir(parents=True, exist_ok=True)
        np.savez(cache, points=mesh[0], triangles=mesh[1], speeds=mesh[2], key=key)
    return mesh


def save_figures(figure, ax, out_dir: "str | Path", day: int,
                 opts: Optional[MicroOptions] = None, log=None) -> List[Path]:
    """Write the full, zoom-window and zoom PNGs from one rendered figure.

    All three come from the same artists - the zoom only changes the axes
    limits - so no trajectory geometry is rebuilt or traversed a second time.

    The figure is mutated in place and is left in the zoomed state (limits
    tightened, chrome hidden); it is meant to be closed afterwards, as
    :func:`plot_microscopic_trajectories` does.
    """
    opts = opts or MicroOptions()
    log = log or (lambda *_: None)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    direction_name = "west" if opts.direction < 0 else "east"
    lane_name = f"lane{opts.lane}" if opts.lane else "laneall"
    # `_py` keeps these from overwriting the MATLAB figures of the same stage,
    # which live in this same folder.
    date = datetime(2022, 11, day).strftime("%Y%m%d")
    stem = f"fig_motion_trajectories_{date}_{direction_name}_{lane_name}_py"
    written: List[Path] = []

    def write(suffix: str) -> None:
        path = out_dir / f"{stem}_{suffix}.png"
        log(f"    writing {path.name}")
        figure.savefig(path, dpi=opts.dpi, facecolor=figure.get_facecolor())
        written.append(path)

    write("lowres")

    zoom_t, zoom_x = _zoom_window(opts)
    rectangle, = ax.plot([zoom_t[0], zoom_t[1], zoom_t[1], zoom_t[0], zoom_t[0]],
                         [zoom_x[0], zoom_x[0], zoom_x[1], zoom_x[1], zoom_x[0]],
                         "-", linewidth=2.5, color="yellow")
    write("zoomwin")

    # The zoom itself: same artists, tighter limits, chrome removed.
    rectangle.set_linewidth(10.0)
    ax.set_xlim(*zoom_t)
    ax.set_ylim(zoom_x[1], zoom_x[0])      # stay inverted, as `axis ij`
    ax.set_axis_off()
    ax.set_position([0.0, 0.0, 1.0, 1.0])
    for colorbar_ax in [a for a in figure.axes if a is not ax]:
        colorbar_ax.set_visible(False)
    figure.patch.set_facecolor("black")
    if opts.zoom_width_factor is not None:
        width, height = opts.fig_size_in
        figure.set_size_inches(width * opts.zoom_width_factor, height)
    write("zoom")

    return written


def _zoom_window(opts: MicroOptions) -> Tuple[Tuple[float, float], Tuple[float, float]]:
    """Zoom limits in plot coordinates: (seconds from window start, km)."""
    start_hhmm, end_hhmm = opts.zoom_time_military
    window_minutes = _minutes(opts.time_window[1]) - _minutes(opts.time_window[0])

    def offset(hhmm: int) -> float:
        minutes = (hhmm // 100) * 60 + hhmm % 100 - _minutes(opts.time_window[0])
        return min(window_minutes, max(0.0, minutes)) * 60.0

    return ((offset(start_hhmm), offset(end_hhmm)),
            (opts.zoom_x_m[0] / 1000.0, opts.zoom_x_m[1] / 1000.0))


def _minutes(hms: str) -> float:
    hour, minute, _ = (int(part) for part in hms.split(":"))
    return hour * 60 + minute


def plot_microscopic_trajectories(slim_dir: "str | Path", day: int,
                                  out_dir: "str | Path",
                                  opts: Optional[MicroOptions] = None,
                                  cache_dir: "str | Path | None" = None,
                                  log=None) -> List[Path]:
    """Render and save the stage-5b figures for one day.

    Convenience wrapper over :func:`render_day` and :func:`save_figures`;
    returns the paths written.
    """
    import matplotlib.pyplot as plt

    opts = opts or MicroOptions()
    figure, ax = render_day(slim_dir, day, opts, cache_dir=cache_dir, log=log)
    try:
        return save_figures(figure, ax, out_dir, day, opts, log=log)
    finally:
        plt.close(figure)
