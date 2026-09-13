"""AV-effect fuel analysis: bin samples by distance to the nearest engaged AV.

Port of the numerical core of ``plot_AV_analysis.m`` — the binning that produces
the article's fuel-consumption-versus-distance results (Figure 2, SM2, SM3). The
figure layout and regression overlays are presentation and are not ported here;
this reproduces the per-bin statistics the figures are drawn from.

For one day it filters the pooled samples (from ``generate_data_samples``) by x
position and time window, bins them by signed distance to the nearest engaged
AV, and for each bin computes:

* ``effective`` fuel consumption = sum(fuel rate) / max(sum(speed), 1e-6)
* ``mean`` fuel consumption over samples with speed >= min_speed_mean
* ``median`` fuel consumption
* the sample count
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict

import numpy as np

__all__ = ["AnalysisOptions", "bin_edges", "bin_samples", "bin_samples_from_mat"]


@dataclass(frozen=True)
class AnalysisOptions:
    """Constants from the head of plot_AV_analysis.m."""

    max_dist: float = 350.0           # [m] analysis half-window
    av_loc_buffer: float = 1.0        # [m] exclusion around the AV itself
    x_window: tuple = (350.0 - 122.0, 7000.0)   # [m] XWINDOW = [MAXDIST-122, 7000]
    time_window: tuple = (645, 915)   # military time
    min_speed_mean: float = 1.0       # [m/s] for the mean statistic


def bin_edges(options: AnalysisOptions = AnalysisOptions()) -> np.ndarray:
    """Signed-distance bin edges, with the +/- buffer inserted around zero.

    Mirrors the ``distToA`` construction: a 10 m grid over [-max_dist, max_dist]
    with -av_loc_buffer and +av_loc_buffer spliced in at the centre.
    """
    md = options.max_dist
    grid = np.round(np.linspace(-md, md, int(np.ceil(md * 2 / 10) + 1)))
    half = grid.size // 2
    return np.concatenate((grid[:half], [-options.av_loc_buffer, options.av_loc_buffer],
                           grid[half + 1:]))


def _time_window_seconds(options: AnalysisOptions) -> tuple:
    lo, hi = options.time_window
    return tuple(((value // 100 - 6) * 60 + value % 100) * 60 for value in (lo, hi))


def bin_samples(samples: Dict[str, np.ndarray],
                options: AnalysisOptions = AnalysisOptions()) -> dict:
    """Filter and bin one day's samples; return per-bin statistics.

    ``samples`` is the dict of arrays from ``mvtpy.samples`` (or a released
    ``.mat``). Returns ``edges``, ``centers``, and length-(n_bins) arrays
    ``effective``, ``mean``, ``median``, ``count``.
    """
    edges = bin_edges(options)
    centers = (edges[:-1] + edges[1:]) / 2

    x = np.asarray(samples["samples_xpos"], dtype=float)
    t = np.asarray(samples["samples_t"], dtype=float)
    keep = (x >= options.x_window[0]) & (x <= options.x_window[1])
    t_lo, t_hi = _time_window_seconds(options)
    keep &= (t >= t_lo) & (t <= t_hi)

    dist = np.asarray(samples["samples_dist"], dtype=float)[keep]
    speed = np.asarray(samples["samples_speed"], dtype=float)[keep]
    rate = np.asarray(samples["samples_fr"], dtype=float)[keep]
    fcons = np.asarray(samples["samples_fcons"], dtype=float)[keep]

    n_bins = centers.size
    effective = np.zeros(n_bins)
    mean = np.zeros(n_bins)
    median = np.zeros(n_bins)
    count = np.zeros(n_bins, dtype=np.int64)

    # bin index per sample: edges[b] <= d < edges[b+1]
    which = np.searchsorted(edges, dist, side="right") - 1
    order = np.argsort(which, kind="stable")
    which_sorted = which[order]
    starts = np.searchsorted(which_sorted, np.arange(n_bins), side="left")
    stops = np.searchsorted(which_sorted, np.arange(n_bins), side="right")

    for b in range(n_bins):
        idx = order[starts[b]:stops[b]]
        count[b] = idx.size
        bin_speed = speed[idx]
        bin_rate = rate[idx]
        bin_fcons = fcons[idx]
        effective[b] = bin_rate.sum() / max(bin_speed.sum(), 1e-6)
        fast = bin_fcons[bin_speed >= options.min_speed_mean]
        mean[b] = fast.mean() if fast.size else np.nan
        median[b] = np.median(bin_fcons) if bin_fcons.size else np.nan

    return {"edges": edges, "centers": centers, "effective": effective,
            "mean": mean, "median": median, "count": count}


def bin_samples_from_mat(mat_path, options: AnalysisOptions = AnalysisOptions()) -> dict:
    """Bin a released ``samples_for_distance_analysis_*.mat`` (needs h5py)."""
    import h5py

    with h5py.File(Path(mat_path), "r") as handle:
        samples = {name: handle[name][:].ravel()
                   for name in ("samples_dist", "samples_speed", "samples_fr",
                                "samples_fcons", "samples_xpos", "samples_t")}
    return bin_samples(samples, options)
