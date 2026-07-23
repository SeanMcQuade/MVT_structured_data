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
"""

from . import fuel, kinematics, lanes, matjson, matround, rawio

__all__ = ["fuel", "kinematics", "lanes", "matjson", "matround", "rawio"]
__version__ = "0.1.0"
