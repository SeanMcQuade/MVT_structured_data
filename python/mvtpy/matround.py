"""MATLAB rounding semantics.

The pipeline rounds every numeric field to four decimals immediately before
encoding (``round(x, 4, 'decimals')`` in ``generate_data_mvt_slim.m``). That one
step is what makes byte-identical output between MATLAB and Python realistic:
it erases floating-point differences below 5e-5 that would otherwise appear in
the last digits. Reproducing it exactly requires two departures from the
obvious implementation.

**Ties go away from zero.** MATLAB's ``round`` rounds halves away from zero
(``round(0.5) == 1``), while Python's ``round`` and NumPy's ``round`` use
banker's rounding (``round(0.5) == 0``).

**Near-ties are snapped to the tie.** MATLAB compensates for binary
representation error: a value whose scaled form lands within one ULP below the
midpoint is treated as a tie and rounded away from zero. This is the same
behavior that makes ``round(2.675, 2)`` return ``2.68`` in MATLAB even though
2.675 is stored as 2.67499999999999982236431605997495353221893310546875.

That second rule is not a rounding error to be corrected - it decides real
bytes in the released data. For the POSIX timestamp 1668600000.1863499 the
exact scaled value is 16686000001863.4986877, one ULP (0.001953125) below the
midpoint, and MATLAB writes .1864 where exact decimal rounding gives .1863.
The tolerance was fitted against 80 values evaluated by MATLAB R2025b itself
and is bracketed on both sides: 0.5 ULP mismatches 40 of them, 2 ULP mismatches
14, 1 ULP matches all 80.
"""

from __future__ import annotations

import math
from typing import Any

try:  # NumPy is optional for this module
    import numpy as _np
except ImportError:  # pragma: no cover - exercised only without NumPy
    _np = None

__all__ = ["round_half_away", "round_decimals", "TIE_TOLERANCE_ULPS"]

#: Distance below the midpoint, in ULPs of the scaled value, that MATLAB still
#: treats as a tie. Fitted against MATLAB R2025b; see the module docstring.
TIE_TOLERANCE_ULPS = 1.0


def round_half_away(value: float) -> float:
    """Round to the nearest integer, halves away from zero (MATLAB ``round``)."""
    return math.floor(value + 0.5) if value >= 0 else math.ceil(value - 0.5)


def round_decimals(value: Any, decimals: int = 4) -> Any:
    """MATLAB's ``round(x, decimals)`` for scalars, sequences, and NumPy arrays.

    Scales by ``10**decimals``, rounds halves - and near-ties within one ULP -
    away from zero, then scales back.
    """
    scale = 10.0**decimals

    if _np is not None and isinstance(value, _np.ndarray):
        values = _np.asarray(value, dtype=float)
        scaled = values * scale
        magnitude = _np.abs(scaled)

        with _np.errstate(invalid="ignore"):
            floor = _np.floor(magnitude)
            fraction = magnitude - floor
            ulp = _np.nextafter(magnitude, _np.inf) - magnitude
            away = (fraction > 0.5) | (_np.abs(fraction - 0.5) <= TIE_TOLERANCE_ULPS * ulp)
            rounded = _np.copysign(floor + away, scaled)

        result = rounded / scale
        return _np.where(_np.isfinite(scaled), result, values)

    if isinstance(value, (list, tuple)):
        return type(value)(round_decimals(item, decimals) for item in value)

    number = float(value)
    if not math.isfinite(number):
        return number

    scaled = number * scale
    magnitude = abs(scaled)
    floor = math.floor(magnitude)
    fraction = magnitude - floor
    ulp = math.nextafter(magnitude, math.inf) - magnitude

    if fraction > 0.5 or abs(fraction - 0.5) <= TIE_TOLERANCE_ULPS * ulp:
        floor += 1
    return math.copysign(floor, scaled) / scale
