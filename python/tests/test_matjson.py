"""Tests for MATLAB-compatible JSON encoding.

Two layers:

* unit tests for the formatting rules, using values taken from the released
  data set (these always run);
* parity tests against the real MATLAB output under ``results/`` (skipped
  automatically when that data is not present, e.g. on a machine that only has
  the repository).
"""

from __future__ import annotations

import math
import re
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from mvtpy import matjson, matround  # noqa: E402
from tools.json_parity import extract_prefix_records  # noqa: E402

RESULTS = Path(__file__).resolve().parents[3] / "results"
STRING_LITERAL = re.compile(r'"(?:[^"\\]|\\.)*"')
NUMBER_LITERAL = re.compile(r"-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?")


# ---------------------------------------------------------------------------
# formatting rules


@pytest.mark.parametrize(
    "value,expected",
    [
        (0.0, "0"),
        (-0.0, "-0"),          # road_grade_radians is full of these
        (840.0, "840"),
        (-1.0, "-1"),
        (5177.3133, "5177.3133"),
        (268983.3831, "268983.3831"),      # largest fixed-notation value seen
        (0.0001, "0.0001"),                # smallest; exponent -4 stays fixed
        (1668599999.9836, "1.6685999999836E+9"),
        (1668600000.0, "1.6686E+9"),       # integral, but still scientific
        (1668691841.800001, "1.6686918418000009E+9"),  # 15 fails -> 17, not 16
    ],
)
def test_number_formats(value, expected):
    assert matjson.encode_number(value) == expected


def test_non_finite_becomes_null():
    assert matjson.encode_number(math.nan) == "null"
    assert matjson.encode_number(math.inf) == "null"
    assert matjson.encode_number(-math.inf) == "null"


def test_structure_has_no_whitespace():
    document = [{"trajectory_id": {"x_oid": "abc-0"}, "timestamp": [1.5, 2.5], "flag": True}]
    assert matjson.dumps(document) == (
        '[{"trajectory_id":{"x_oid":"abc-0"},"timestamp":[1.5,2.5],"flag":true}]'
    )


def test_empty_array_and_null():
    assert matjson.dumps({"downstream_av_id": [], "missing": None}) == (
        '{"downstream_av_id":[],"missing":null}'
    )


def test_loads_preserves_negative_zero():
    decoded = matjson.loads('{"road_grade_radians":[-0,0]}')
    assert matjson.dumps(decoded) == '{"road_grade_radians":[-0,0]}'


def test_round_half_away_from_zero():
    # MATLAB rounds halves away from zero; Python's round() would give 0 and 2.
    assert matround.round_half_away(0.5) == 1
    assert matround.round_half_away(-0.5) == -1
    assert matround.round_decimals(1.00005, 4) == pytest.approx(1.0001, abs=0)


# ---------------------------------------------------------------------------
# parity against real MATLAB output


def _released_files():
    candidates = [
        RESULTS / "slim" / "2022-11-16" / "I-24MOTION_2022-11-16_05-59-59.json",
        RESULTS / "full" / "2022-11-17" / "I-24MOTION_2022-11-17_06-49-59.json",
        RESULTS / "gps" / "CIRCLES_GPS_10Hz_2022-11-16.json",
    ]
    return [path for path in candidates if path.is_file()]


@pytest.mark.parametrize("path", _released_files() or [pytest.param(None, marks=pytest.mark.skip(
    reason="released results/ data not available"))])
def test_number_tokens_match_matlab(path):
    """Every numeric literal MATLAB wrote must re-encode to the same text."""
    with path.open("rb") as handle:
        text = handle.read(8 * 1024 * 1024).decode("utf-8", errors="ignore")
    text = STRING_LITERAL.sub('""', text).rsplit(",", 1)[0]

    mismatches = [
        (token, matjson.encode_number(float(token)))
        for token in NUMBER_LITERAL.findall(text)
        if matjson.encode_number(float(token)) != token
    ]
    assert not mismatches[:5], f"{len(mismatches)} tokens differ, e.g. {mismatches[:5]}"


@pytest.mark.parametrize("path", _released_files() or [pytest.param(None, marks=pytest.mark.skip(
    reason="released results/ data not available"))])
def test_document_roundtrip_is_byte_identical(path):
    """Decoding and re-encoding real MATLAB output must reproduce it exactly."""
    with path.open("rb") as handle:
        raw = handle.read(90 * 1024 * 1024).decode("utf-8", errors="ignore")

    prefix = extract_prefix_records(raw, 3)
    assert matjson.dumps(matjson.loads(prefix)) == prefix
