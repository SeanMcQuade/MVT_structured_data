"""Tests for the AV-effect binning (numerical core of plot_AV_analysis).

The per-bin statistics are checked against a committed MATLAB reference
(``tests/fixtures/av_bin_stats_2022-11-16.npz``, produced by the same binning in
MATLAB). Median and count are bit-exact; the summed statistics agree to
floating-point roundoff.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from mvtpy import avanalysis  # noqa: E402

WORKSPACE = Path(__file__).resolve().parents[3]
from conftest import RESULTS_DIR  # noqa: E402
SAMPLES_MAT = (RESULTS_DIR / "figures" / "2022-11-16"
               / "samples_for_distance_analysis_16.mat")
BIN_FIXTURE = (Path(__file__).resolve().parent / "fixtures"
               / "av_bin_stats_2022-11-16.npz")


def test_bin_edges_structure():
    edges = avanalysis.bin_edges()
    # 71 bins -> 72 edges, symmetric, with +/- buffer straddling zero.
    assert edges.size == 72
    assert edges[0] == -350.0 and edges[-1] == 350.0
    assert -1.0 in edges and 1.0 in edges
    centers = (edges[:-1] + edges[1:]) / 2
    assert centers.size == 71


def test_time_window_conversion():
    # 06:45-09:15 military -> seconds after 06:00.
    lo, hi = avanalysis._time_window_seconds(avanalysis.AnalysisOptions())
    assert lo == 2700 and hi == 11700


def test_binning_on_synthetic_samples():
    # Two samples in known bins with a controlled x/time window.
    samples = {
        "samples_dist": np.array([-5.0, 200.0]),
        "samples_speed": np.array([10.0, 20.0]),
        "samples_fr": np.array([1.0, 4.0]),
        "samples_fcons": np.array([0.1, 0.2]),
        "samples_xpos": np.array([300, 300], dtype=np.int16),
        "samples_t": np.array([3000, 3000], dtype=np.uint16),
    }
    out = avanalysis.bin_samples(samples)
    assert out["count"].sum() == 2
    # The -5 m sample lands in the [-10, -1) bin; effective = fr/speed there.
    b = np.searchsorted(out["edges"], -5.0, side="right") - 1
    assert out["effective"][b] == pytest.approx(1.0 / 10.0)


@pytest.mark.skipif(not BIN_FIXTURE.is_file() or not SAMPLES_MAT.is_file(),
                    reason="samples .mat or bin reference not available")
def test_bin_stats_match_matlab():
    """Per-bin stats match the MATLAB reference (needs h5py + the samples .mat)."""
    pytest.importorskip("h5py")
    out = avanalysis.bin_samples_from_mat(SAMPLES_MAT)
    ref = np.load(BIN_FIXTURE)

    assert np.array_equal(out["centers"], ref["centers"])
    assert np.array_equal(out["count"], ref["count"])            # exact
    assert np.array_equal(out["median"], ref["median"])          # exact
    for name in ("effective", "mean"):
        mine, expected = out[name], ref[name]
        assert np.array_equal(np.isnan(mine), np.isnan(expected)), name
        diff = np.abs(np.nan_to_num(mine) - np.nan_to_num(expected))
        assert diff.max() < 1e-12, f"{name}: max diff {diff.max()}"
