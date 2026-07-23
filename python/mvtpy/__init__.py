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
"""

from . import fuel, kinematics, matjson, matround

__all__ = ["fuel", "kinematics", "matjson", "matround"]
__version__ = "0.1.0"
