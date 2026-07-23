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
"""

from . import matjson, matround

__all__ = ["matjson", "matround"]
__version__ = "0.1.0"
