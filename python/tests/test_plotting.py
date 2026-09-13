"""Smoke tests for the matplotlib figure routines.

These check the figures build without error and carry the intended colors and
structure on small synthetic inputs. They are not pixel comparisons - the
figures are deliberately visually close rather than geometrically exact.
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

from mvtpy import plotting  # noqa: E402


def test_parula_colormap_endpoints():
    cmap = plotting.parula_cmap()
    # Parula runs dark blue -> yellow.
    assert cmap(0.0)[:3] == pytest.approx(tuple(plotting.PARULA[0]), abs=1e-6)
    assert cmap(1.0)[:3] == pytest.approx(tuple(plotting.PARULA[-1]), abs=1e-6)


def test_day_colors_are_blue_red_green():
    assert plotting.DAY_COLORS[16] == (0.0, 0.0, 1.0)
    assert plotting.DAY_COLORS[17] == (1.0, 0.0, 0.0)
    assert plotting.DAY_COLORS[18] == (0.0, 0.6, 0.0)


def test_plot_field_builds_with_overlay():
    t = 1668600000 + np.arange(0, 200, 5.0)
    x = np.arange(0, 300, 50.0)
    fields = {
        "t": t, "x": x, "direction": -1, "lane": 0,
        "field": {name: np.random.default_rng(0).random((t.size, x.size))
                  for name in ("Rho", "Q", "F", "U", "Phi", "Psi")},
    }
    gps = [{
        "direction": -1,
        "control_car": [0, 1, 1, 0],
        "timestamp": (1668600050 + np.arange(4) * 5.0).tolist(),
        "x_position": [100.0, 150.0, 200.0, 250.0],
    }]
    ax = plotting.plot_field(fields, "Psi", gps_records=gps)
    # colorbar limit and unit come from FIELD_COLOR_LIMITS
    assert ax.collections  # a pcolormesh was drawn
    assert ax.yaxis_inverted()   # position 0 at top, MATLAB orientation
    plt.close(ax.figure)


def test_plot_av_fuel_builds_and_colors_days():
    def stats():
        centers = np.linspace(-350, 350, 71)
        return {"centers": centers,
                "effective": 0.1 + 0.01 * np.sin(centers / 50),
                "mean": np.full(71, 0.1), "median": np.full(71, 0.1),
                "count": np.full(71, 100)}

    axes = plotting.plot_av_fuel({16: stats(), 17: stats(), 18: stats()})
    ahead, behind = axes
    assert ahead.get_title() == "ahead of AV"
    assert behind.get_title() == "behind AV"
    # one labeled line per day in the ahead panel
    labels = [line.get_label() for line in ahead.get_lines() if line.get_label().startswith("Nov")]
    assert set(labels) == {"Nov 16", "Nov 17", "Nov 18"}
    plt.close(ahead.figure)
