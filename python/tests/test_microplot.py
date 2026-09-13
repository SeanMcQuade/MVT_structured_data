"""Tests for the microscopic trajectory figures (stage 5b).

The load-bearing property here is the memory fix: one artist per file rather
than one per trajectory. That is asserted directly, because it is the thing that
regressed in MATLAB and rebooted a machine. The rest are geometry and
labelling checks on small synthetic inputs - these figures are visually close to
the MATLAB output, not pixel-exact.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from mvtpy import microplot  # noqa: E402

DAY = 18
OPTS = microplot.MicroOptions()
T0, _ = OPTS.window_bounds(DAY)


def make_record(start=600.0, duration=60.0, n=200, x0=5000.0, speed=25.0,
                length_ft=18.0, lane=2.0):
    """A westbound trajectory: position decreases at `speed` m/s."""
    t = T0 + start + np.linspace(0.0, duration, n)
    x = x0 - speed * (t - t[0])
    return {"timestamp": t.tolist(), "x_position_meters": x.tolist(),
            "length": length_ft, "lane_number": lane}


# --- the memory fix -------------------------------------------------------


def test_one_artist_per_batch_not_per_trajectory():
    """The whole point of the fix: artist count is O(files), not O(trajectories)."""
    records = [make_record(start=600.0 + i) for i in range(500)]
    mesh = microplot._concatenate(
        [microplot.trajectory_mesh(r, OPTS, T0) for r in records])

    fig, ax = plt.subplots()
    microplot.add_mesh(ax, mesh, OPTS)
    assert len(ax.collections) == 1
    plt.close(fig)


def test_mesh_holds_no_per_triangle_python_objects():
    """Geometry lives in flat arrays; a PolyCollection would make one Path each."""
    records = [make_record(start=600.0 + i) for i in range(50)]
    points, triangles, speeds = microplot._concatenate(
        [microplot.trajectory_mesh(r, OPTS, T0) for r in records])

    assert isinstance(points, np.ndarray) and points.ndim == 2
    assert isinstance(triangles, np.ndarray) and triangles.dtype == np.int32
    assert speeds.shape[0] == points.shape[0]
    # Indices stay within the concatenated point array.
    assert triangles.min() >= 0 and triangles.max() < points.shape[0]


# --- geometry -------------------------------------------------------------


def test_subsampling_follows_skip_t_plot():
    n = 200
    for skip in (5, 50, 500):
        index = microplot._subsample_indices(n, skip)
        assert index.size == max(int(np.ceil(n / skip)), 2)
        assert index[0] == 0 and index[-1] == n - 1
        assert np.all(np.diff(index) >= 0)


def test_ribbon_spans_vehicle_length_and_is_speed_colored():
    length_ft, speed = 18.0, 25.0
    points, triangles, speeds = microplot.trajectory_mesh(
        make_record(length_ft=length_ft, speed=speed), OPTS, T0)

    half = points.shape[0] // 2
    near, far = points[:half, 1], points[half:, 1]
    # Westbound (direction -1): the far edge trails by the vehicle length, in km.
    assert np.allclose(near - far, length_ft * microplot.FT2M / 1000.0)
    # Speed is positive travelling westbound, and both edges carry it.
    assert np.allclose(speeds, speed, rtol=1e-6)
    assert triangles.shape[0] == 2 * (half - 1)


def test_time_is_seconds_from_window_start():
    start = 900.0
    points, _, _ = microplot.trajectory_mesh(make_record(start=start), OPTS, T0)
    assert points[:, 0].min() == pytest.approx(start)


def test_short_and_malformed_trajectories_are_dropped():
    assert microplot.trajectory_mesh(
        {"timestamp": [1.0], "x_position_meters": [2.0], "length": 18.0}, OPTS, T0) is None
    assert microplot.trajectory_mesh(
        {"timestamp": [1.0, 2.0], "x_position_meters": [2.0], "length": 18.0},
        OPTS, T0) is None


def test_zero_duration_samples_do_not_produce_nan_speeds():
    """Repeated timestamps divide by zero; those speeds must not reach the colormap."""
    t = [T0 + 600.0] * 4
    record = {"timestamp": t, "x_position_meters": [10.0, 9.0, 8.0, 7.0], "length": 18.0}
    _, _, speeds = microplot.trajectory_mesh(record, OPTS, T0)
    assert np.all(np.isfinite(speeds))


def test_lane_and_direction_selection():
    from dataclasses import replace

    lane_two = make_record(lane=2.0)
    assert microplot._keep(lane_two, OPTS)                       # lane 0 = all
    assert microplot._keep(lane_two, replace(OPTS, lane=2))
    assert not microplot._keep(lane_two, replace(OPTS, lane=3))
    # Slim segments carry no direction field, so direction filtering is skipped.
    assert microplot._keep(dict(lane_two, direction=-1.0), OPTS)
    assert not microplot._keep(dict(lane_two, direction=1.0), OPTS)


# --- figure ---------------------------------------------------------------


def test_axes_orientation_ticks_and_window():
    fig, ax = plt.subplots()
    microplot.setup_axes(ax, DAY, OPTS)

    low, high = ax.get_ylim()
    assert low > high                       # MATLAB `axis ij`: position 0 on top
    assert ax.get_xlim() == (600.0, 13800.0)   # 10 min trimmed from each end
    assert len(ax.get_xticks()) == OPTS.n_xticks + 1
    assert "Westbound (all lanes)" in ax.get_title()
    assert "18-Nov-2022" in ax.get_title()
    plt.close(fig)


def test_apply_limits_is_idempotent():
    """render_day re-asserts limits after tripcolor autoscales; twice must equal once."""
    fig, ax = plt.subplots()
    microplot.setup_axes(ax, DAY, OPTS)
    once = (ax.get_xlim(), ax.get_ylim())
    microplot.apply_limits(ax, DAY, OPTS)
    assert (ax.get_xlim(), ax.get_ylim()) == once
    plt.close(fig)


def test_speed_colormap_runs_red_to_green():
    cmap = microplot.speed_cmap()
    assert cmap(0.0)[:3] == pytest.approx((1.0, 0.0, 0.0))
    assert cmap(1.0)[:3] == pytest.approx((0.0, 1.0, 0.0))


def test_colorbar_limits_match_its_label():
    """MATLAB's int32 path labels a 0-3.5e4 bar 'm/s'; this port must not."""
    fig, ax = plt.subplots()
    mesh = microplot._concatenate([microplot.trajectory_mesh(make_record(), OPTS, T0)])
    mappable = microplot.add_mesh(ax, mesh, OPTS)
    assert mappable.get_clim() == (0.0, OPTS.speed_limit)
    plt.close(fig)


def test_zoom_window_maps_military_time_to_plot_coordinates():
    (t_lo, t_hi), (x_lo, x_hi) = microplot._zoom_window(OPTS)
    assert (t_lo, t_hi) == (720.0, 1140.0)      # 06:12 and 06:19 past 06:00
    assert (x_lo, x_hi) == (0.845, 1.495)


def test_save_figures_writes_three_distinct_pngs(tmp_path):
    fig, ax = plt.subplots(figsize=(4, 2))
    microplot.setup_axes(ax, DAY, OPTS)
    mesh = microplot._concatenate(
        [microplot.trajectory_mesh(make_record(start=700.0 + 3 * i, x0=1400.0, speed=5.0),
                                   OPTS, T0) for i in range(40)])
    microplot.add_mesh(ax, mesh, OPTS)

    written = microplot.save_figures(fig, ax, tmp_path, DAY, OPTS)
    plt.close(fig)

    assert [p.name for p in written] == [
        f"fig_motion_trajectories_20221118_west_laneall_py_{suffix}.png"
        for suffix in ("lowres", "zoomwin", "zoom")]
    assert all(p.stat().st_size > 0 for p in written)
    # The three renders must differ: the zoom window is added, then zoomed into.
    assert len({p.read_bytes() for p in written}) == 3


def test_rendered_zoom_is_not_blank(tmp_path):
    """The MATLAB failure wrote a background-only zoom; guard against that here."""
    fig, ax = plt.subplots(figsize=(4, 2))
    microplot.setup_axes(ax, DAY, OPTS)
    mesh = microplot._concatenate(
        [microplot.trajectory_mesh(make_record(start=700.0 + 3 * i, x0=1400.0, speed=5.0),
                                   OPTS, T0) for i in range(40)])
    microplot.add_mesh(ax, mesh, OPTS)
    written = microplot.save_figures(fig, ax, tmp_path, DAY, OPTS)
    plt.close(fig)

    image = plt.imread(written[-1])[:, :, :3]
    colors = np.unique(image.reshape(-1, 3), axis=0)
    assert len(colors) > 8, "zoom render carries only background colors"


# --- end to end on a synthetic segment tree -------------------------------


def test_render_day_one_artist_per_file_and_uses_cache(tmp_path):
    from mvtpy import matjson

    slim = tmp_path / "slim"
    slim.mkdir()
    for index in range(3):
        records = [make_record(start=600.0 + 30 * index + i, x0=6000.0 - 100 * i)
                   for i in range(5)]
        (slim / f"I-24MOTION_2022-11-18_0{index}-00-00.json").write_text(
            matjson.dumps(records), encoding="utf-8")

    cache = tmp_path / "cache"
    fig, ax = microplot.render_day(slim, DAY, OPTS, cache_dir=cache)
    assert len(ax.collections) == 3          # one artist per file, not per trajectory
    plt.close(fig)

    caches = sorted(cache.glob("*_micro.npz"))
    assert len(caches) == 3

    # Second pass reads the caches and produces the same geometry.
    fig, ax = microplot.render_day(slim, DAY, OPTS, cache_dir=cache)
    assert len(ax.collections) == 3
    plt.close(fig)


def test_cache_is_rebuilt_when_geometry_options_change(tmp_path):
    """A skip_t_plot change must not silently reuse a mesh built at another one."""
    from dataclasses import replace

    from mvtpy import matjson

    slim = tmp_path / "slim"
    slim.mkdir()
    (slim / "I-24MOTION_2022-11-18_00-00-00.json").write_text(
        matjson.dumps([make_record(n=500)]), encoding="utf-8")
    cache = tmp_path / "cache"

    coarse = microplot._cached_file_mesh(
        next(slim.glob("*.json")), OPTS, T0, cache, log=lambda *_: None)
    dense = microplot._cached_file_mesh(
        next(slim.glob("*.json")), replace(OPTS, skip_t_plot=5), T0, cache,
        log=lambda *_: None)

    assert dense[1].shape[0] > coarse[1].shape[0]


def test_file_mesh_orders_longest_vehicles_first(tmp_path):
    from mvtpy import matjson

    path = tmp_path / "I-24MOTION_2022-11-18_00-00-00.json"
    path.write_text(matjson.dumps([make_record(length_ft=10.0, x0=5000.0),
                                   make_record(length_ft=40.0, x0=4000.0),
                                   make_record(length_ft=25.0, x0=3000.0)]),
                    encoding="utf-8")

    points, _, _ = microplot.file_mesh(path, OPTS, T0)
    half = points.shape[0] // 6          # points per trajectory / 2 edges
    # Ribbon thickness in km, in the order the trajectories were meshed.
    thickness = [points[start:start + half, 1][0] - points[start + half:start + 2 * half, 1][0]
                 for start in range(0, points.shape[0], 2 * half)]
    assert thickness == sorted(thickness, reverse=True)


def test_render_day_without_segments_is_an_error(tmp_path):
    with pytest.raises(FileNotFoundError):
        microplot.render_day(tmp_path, DAY, OPTS)
