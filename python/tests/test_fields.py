"""Tests for the macroscopic-field accumulation.

The difference-array accumulation is checked against a direct box-add on a
synthetic case (bit-close), the grid and derived-field rules with unit tests,
and the whole field against the released ``.mat`` in
``test_fields_match_released`` (to floating-point roundoff, with identical NaN
layout) when the data is available.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from mvtpy import fields  # noqa: E402

WORKSPACE = Path(__file__).resolve().parents[3]
SLIM_DIR = WORKSPACE / "results" / "slim" / "2022-11-16"
REFERENCE = (WORKSPACE / "results" / "figures" / "2022-11-16"
             / "fields_motion_2022-11-16.mat")


def test_grid_dimensions():
    out_t = fields._colon(0.0, 5.0, 14400.0)
    assert out_t.size == 2881
    out_x = fields._colon(0.0, 50.0, 6500.0)
    assert out_x.size == 131


def test_difference_array_equals_direct_box_add():
    rng = np.random.default_rng(0)
    t = np.arange(0, 100, 5.0)
    x = np.arange(0, 300, 50.0)
    ht, hx = 5.0, 100.0
    tp = rng.uniform(5, 95, 40)
    xp = rng.uniform(0, 300, 40)
    val = rng.uniform(-2, 3, 40)

    direct = np.zeros((t.size, x.size))
    for i in range(tp.size):
        direct[np.ix_(np.abs(t - tp[i]) <= ht, np.abs(x - xp[i]) <= hx)] += val[i]

    diff = np.zeros((t.size + 1, x.size + 1))
    fields._scatter_rectangles(
        diff,
        np.searchsorted(t, tp - ht, "left"), np.searchsorted(t, tp + ht, "right"),
        np.searchsorted(x, xp - hx, "left"), np.searchsorted(x, xp + hx, "right"),
        val)
    got = fields._integrate(diff, t.size, x.size)
    assert np.allclose(got, direct, atol=1e-12, rtol=0)


def test_box_bounds_are_inclusive():
    t = np.arange(0, 100, 5.0)
    ht = 5.0
    # A point on a grid node includes the node and one node either side.
    lo = int(np.searchsorted(t, 10.0 - ht, "left"))
    hi = int(np.searchsorted(t, 10.0 + ht, "right"))
    assert list(range(lo, hi)) == [1, 2, 3]      # t = 5, 10, 15


def test_derived_fields_and_nan_masks():
    records = [{
        "direction": -1,
        "lane_number": 3,
        "timestamp": [0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
        "x_position_meters": [10.0, 12.0, 14.0, 16.0, 18.0, 20.0],
        "fuel_rate_grams_per_second": [0.5] * 6,
    }]
    out = fields.macroscopic_fields(records, 16)
    rho, q, f = out["field"]["Rho"], out["field"]["Q"], out["field"]["F"]
    u, phi, psi = out["field"]["U"], out["field"]["Phi"], out["field"]["Psi"]

    # U, Phi are NaN exactly where Rho < 1e-3; Psi where Q < 1e-2.
    assert np.array_equal(np.isnan(u), rho < 1e-3)
    assert np.array_equal(np.isnan(phi), rho < 1e-3)
    assert np.array_equal(np.isnan(psi), q < 1e-2)
    # Where defined, the derived fields are the documented ratios.
    good = rho >= 1e-3
    assert np.allclose(u[good], q[good] / rho[good], equal_nan=True)


@pytest.mark.skipif(not REFERENCE.is_file() or not SLIM_DIR.is_dir(),
                    reason="released slim data / fields reference not available")
def test_fields_match_released():
    """Fields match the released .mat to float roundoff, NaN layout identical."""
    import scipy.io as sio

    out = fields.macroscopic_fields_from_dir(SLIM_DIR, 16)
    reference = sio.loadmat(REFERENCE, struct_as_record=False, squeeze_me=True)["field"].value

    assert np.array_equal(out["t"], sio.loadmat(REFERENCE, squeeze_me=True)["t"])
    for name in ("Rho", "Q", "F", "U", "Phi", "Psi"):
        mine = out["field"][name]
        ref = np.asarray(getattr(reference, name))
        assert np.array_equal(np.isnan(mine), np.isnan(ref)), f"{name}: NaN layout"
        diff = np.abs(np.nan_to_num(mine) - np.nan_to_num(ref))
        assert diff.max() < 1e-9, f"{name}: max abs diff {diff.max()}"
