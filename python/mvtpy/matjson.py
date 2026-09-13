"""MATLAB-compatible JSON encoding.

The MVT pipeline writes its released data with MATLAB's ``jsonencode``. To let a
Python port produce *byte-identical* files, this module reproduces that
encoder's output exactly, rather than relying on :mod:`json`, whose number
formatting and spacing differ.

Rules implemented here, each calibrated against the released data set (see
``docs/MATLAB_JSON_FORMAT.md`` for the evidence, and
``python/matlab_probes/probe_jsonencode.m`` for the script that re-verifies
them inside MATLAB):

* no whitespace anywhere: ``{"a":1,"b":[1,2]}``
* object keys keep insertion order (MATLAB struct field order)
* ``NaN`` and ``Inf`` encode as ``null``
* ``True``/``False`` encode as ``true``/``false``
* integral doubles print without a decimal point (``840``, ``-0``)
* digits come from ``%.15g``; if that does not round-trip, from ``%.17g``.
  MATLAB never uses 16 significant digits, which is why the shortest
  round-trip representation (Python's ``repr``) is *not* a substitute: for
  1668691841.800001 MATLAB writes ``1.6686918418000009E+9``, not the shorter
  ``1.668691841800001E+9``. Verified against ~10 million released values.
* values print in fixed notation, switching to scientific
  ``1.6685999999836E+9`` outside the exponent range below

The scientific-notation switch is expressed as a decimal-exponent band. The
lower bound is pinned by the data (``0.0001`` prints fixed, so the switch is at
exponent < -4, matching C's ``%g``). The upper bound is pinned from above at 9
(``1.6686E+9``, an integral value, prints scientific) and from below at 5
(``268983.3831`` prints fixed); values between 1e6 and 1e9 do not occur in the
released data, so 9 is the calibrated value and the probe script confirms it.
"""

from __future__ import annotations

import json
import math
from typing import Any, Mapping, Sequence

__all__ = ["dumps", "loads", "encode_number", "SCI_EXP_MIN", "SCI_EXP_MAX"]

#: Scientific notation is used when the decimal exponent is below this value.
SCI_EXP_MIN = -4
#: Scientific notation is used when the decimal exponent is at or above this.
#: Confirmed by probe: 123456.789 prints fixed, 1e6 prints as "1.0E+6".
SCI_EXP_MAX = 6

_ESCAPES = {
    '"': '\\"',
    "\\": "\\\\",
    "\b": "\\b",
    "\f": "\\f",
    "\n": "\\n",
    "\r": "\\r",
    "\t": "\\t",
}


def dumps(value: Any) -> str:
    """Encode ``value`` exactly as MATLAB's ``jsonencode`` would.

    Accepts the usual JSON-compatible Python types: ``dict`` (object, key order
    preserved), ``list``/``tuple`` (array), ``str``, ``bool``, ``int``,
    ``float``, and ``None``. NumPy scalars and arrays are accepted when NumPy is
    installed, and are converted through their Python equivalents.
    """
    out: list[str] = []
    _write(value, out)
    return "".join(out)


def loads(text: str) -> Any:
    """Decode JSON the way the pipeline needs it for byte-exact round-trips.

    Differs from :func:`json.loads` in one respect that matters here: every
    number is decoded as a ``float``, so ``-0`` survives as ``-0.0``. MATLAB
    writes negative zero (road grade fields are full of it) and the standard
    decoder would turn it into the integer ``0``, changing a byte on re-encode.
    """
    return json.loads(text, parse_int=float, parse_float=float)


def encode_number(value: float) -> str:
    """Format a single numeric value the way MATLAB's ``jsonencode`` does."""
    if isinstance(value, bool):
        return "true" if value else "false"

    number = float(value)

    if math.isnan(number) or math.isinf(number):
        # MATLAB has no JSON spelling for these and emits null.
        return "null"

    if number == 0.0:
        # Negative zero keeps its sign: MATLAB writes -0.
        return "-0" if math.copysign(1.0, number) < 0 else "0"

    digits, exponent = _shortest_digits(number)

    if exponent < SCI_EXP_MIN or exponent >= SCI_EXP_MAX:
        return _scientific(number < 0, digits, exponent)
    return _fixed(number < 0, digits, exponent)


# ---------------------------------------------------------------------------
# structure


def _write(value: Any, out: list[str]) -> None:
    if value is None:
        out.append("null")
        return
    if isinstance(value, bool):
        out.append("true" if value else "false")
        return
    if isinstance(value, str):
        out.append(_string(value))
        return
    if isinstance(value, Mapping):
        out.append("{")
        first = True
        for key, item in value.items():
            if not first:
                out.append(",")
            first = False
            out.append(_string(str(key)))
            out.append(":")
            _write(item, out)
        out.append("}")
        return
    if isinstance(value, (int, float)):
        out.append(encode_number(value))
        return

    converted = _from_numpy(value)
    if converted is not None:
        _write(converted, out)
        return

    if isinstance(value, Sequence):
        out.append("[")
        for index, item in enumerate(value):
            if index:
                out.append(",")
            _write(item, out)
        out.append("]")
        return

    raise TypeError(f"cannot encode object of type {type(value).__name__}")


def _from_numpy(value: Any) -> Any:
    """Convert NumPy scalars/arrays to plain Python, or return None."""
    if type(value).__module__.split(".")[0] != "numpy":
        return None
    if hasattr(value, "ndim") and getattr(value, "ndim") == 0:
        return value.item()
    if hasattr(value, "tolist"):
        return value.tolist()
    return None


def _string(text: str) -> str:
    out = ['"']
    for char in text:
        escape = _ESCAPES.get(char)
        if escape is not None:
            out.append(escape)
        elif ord(char) < 0x20:
            out.append(f"\\u{ord(char):04x}")
        else:
            out.append(char)
    out.append('"')
    return "".join(out)


# ---------------------------------------------------------------------------
# numbers


#: Significant-digit counts MATLAB tries, in order. It uses 15 when that
#: round-trips and 17 otherwise; 16 is never used.
PRECISION_LADDER = (15, 17)


def _shortest_digits(number: float) -> tuple[str, int]:
    """Return (significant digits, decimal exponent) as MATLAB would print them.

    Formats at each precision in :data:`PRECISION_LADDER` and takes the first
    that reproduces the value exactly, then strips trailing zeros. For
    1668599999.9836 this returns ("16685999999836", 9).
    """
    magnitude = abs(number)

    text = ""
    for precision in PRECISION_LADDER:
        text = f"{magnitude:.{precision - 1}e}"
        if float(text) == magnitude:
            break

    mantissa, _, exponent_text = text.partition("e")
    exponent = int(exponent_text)
    digits = mantissa.replace(".", "").rstrip("0") or "0"
    return digits, exponent


def _fixed(negative: bool, digits: str, exponent: int) -> str:
    sign = "-" if negative else ""
    if exponent >= 0:
        if len(digits) <= exponent + 1:
            # Integral value: MATLAB writes no decimal point.
            return sign + digits + "0" * (exponent + 1 - len(digits))
        return sign + digits[: exponent + 1] + "." + digits[exponent + 1 :]
    return sign + "0." + "0" * (-exponent - 1) + digits


def _scientific(negative: bool, digits: str, exponent: int) -> str:
    sign = "-" if negative else ""
    mantissa = digits[0]
    if len(digits) > 1:
        mantissa += "." + digits[1:]
    elif exponent >= 0:
        # Quirk confirmed by probe: a one-digit mantissa is padded for positive
        # exponents ("1.0E+6") but not for negative ones ("1E-5").
        mantissa += ".0"
    exponent_sign = "+" if exponent >= 0 else "-"
    return f"{sign}{mantissa}E{exponent_sign}{abs(exponent)}"
