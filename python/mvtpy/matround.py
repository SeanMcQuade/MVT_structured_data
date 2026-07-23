"""MATLAB rounding semantics.

The pipeline rounds every numeric field to four decimals immediately before
encoding (``round(x, 4, 'decimals')`` in ``generate_data_mvt_slim.m``). That one
step is what makes byte-identical output between MATLAB and Python realistic:
it erases floating-point differences below 5e-5 that would otherwise appear in
the last digits.

Reproducing it requires care, because the two languages break ties differently:

* MATLAB's ``round`` rounds halves **away from zero** (round(0.5) == 1,
  round(-0.5) == -1).
* Python's ``round`` and NumPy's ``round`` use banker's rounding, so
  ``round(0.5) == 0`` and ``np.round(2.5) == 2``.

Ties are rare but real in this data, and a single mismatched tie changes a byte,
so the MATLAB rule is implemented explicitly here.
"""

from __future__ import annotations

import math
from typing import Any

try:  # NumPy is optional for this module
    import numpy as _np
except ImportError:  # pragma: no cover - exercised only without NumPy
    _np = None

__all__ = ["round_half_away", "round_decimals"]


def round_half_away(value: float) -> float:
    """Round to the nearest integer, halves away from zero (MATLAB ``round``)."""
    return math.floor(value + 0.5) if value >= 0 else math.ceil(value - 0.5)


def round_decimals(value: Any, decimals: int = 4) -> Any:
    """MATLAB's ``round(x, decimals)`` for scalars and NumPy arrays.

    Scales by 10**decimals, rounds halves away from zero, and scales back - the
    same operation order MATLAB uses, so the results agree bit for bit.
    """
    scale = 10.0**decimals

    if _np is not None and isinstance(value, _np.ndarray):
        scaled = value * scale
        rounded = _np.where(scaled >= 0,
                            _np.floor(scaled + 0.5),
                            _np.ceil(scaled - 0.5))
        return rounded / scale

    if isinstance(value, (list, tuple)):
        return type(value)(round_decimals(item, decimals) for item in value)

    number = float(value)
    if math.isnan(number) or math.isinf(number):
        return number
    return round_half_away(number * scale) / scale
