"""Python implementation of the CIRCLES MegaVanderTest pipeline.

The MATLAB pipeline under ``Scripts/`` remains the reference implementation.
This package reproduces it stage by stage, with the requirement that the JSON it
writes is byte-identical to MATLAB's output for the same inputs.

Modules
-------
matjson
    MATLAB-compatible JSON encoding (``jsonencode`` byte parity).
matround
    MATLAB rounding semantics (round-half-away-from-zero), used by the 4-decimal
    rounding the pipeline applies before encoding.
kinematics
    Speed, acceleration, road grade, and quadrature, ported expression by
    expression from generate_data_mvt_slim.m.
fuel
    The CIRCLES simplified fuel-consumption models from Models/.
lanes
    Lane identification and lane-change clipping, which decide how a raw
    trajectory becomes released segments.
rawio
    Streaming reader for the multi-gigabyte raw MOTION segment files.
avdist
    Distance from each trajectory to the nearest control vehicle.
slim
    Whole-segment assembly and writing; byte-identical to MATLAB's output.
gpsruns
    Control-vehicle GPS files split into testbed runs (first stage of
    assemble_data_GPS).
gpsassemble
    GPS resampling, control-signal reconstruction, and record assembly (the
    rest of assemble_data_GPS, bar the MOTION-matching bias).
gpsmatch
    The MOTION-matching x_position bias (median_xd), the last piece of
    assemble_data_GPS.
samples
    Per-sample data near engaged control vehicles (generate_data_samples).
fields
    Macroscopic traffic-state fields on a time-space grid
    (generate_macroscopic_fields).
avanalysis
    Distance-binning of fuel samples (numerical core of plot_AV_analysis).
plotting
    Matplotlib figures (field heatmaps, AV fuel curves) - same colors and
    layout as the MATLAB figures, not pixel-exact. Requires matplotlib.
microplot
    Microscopic trajectory time-space figures (plot_microscopic_trajectories),
    with the per-file batching that keeps the artist count bounded.
segments
    Raw segment -> processed filename mapping, resolved before any decode
    (the MATLAB build manifest).
build
    A make for the pipeline: dependency graph, staleness, and a -j scheduler
    over per-day and per-segment work units.
progress
    Live progress display for a build, degrading to plain lines when piped.
verify
    Checks generated outputs against expected md5s, or against another
    results tree.
"""

from . import (avanalysis, avdist, build, fields, fuel, gpsassemble, gpsmatch,
               gpsruns, kinematics, lanes, matjson, matround, microplot,
               plotting, progress, rawio, samples, segments, slim, verify)

__all__ = ["avanalysis", "avdist", "build", "fields", "fuel", "gpsassemble",
           "gpsmatch", "gpsruns", "kinematics", "lanes", "matjson", "matround",
           "microplot", "plotting", "progress", "rawio", "samples", "segments",
           "slim", "verify"]
__version__ = "0.1.0"
