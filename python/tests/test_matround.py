"""Tests for MATLAB rounding semantics.

The near-tie rule is the subtle one: MATLAB treats a scaled value within one
ULP below the midpoint as a tie and rounds it away from zero, compensating for
binary representation error. It decides real bytes in the released data, so it
is pinned here against values evaluated by MATLAB itself
(``matlab_probes/round_probe_2025b.json``, produced by ``probe_round.m``).
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from mvtpy.matround import TIE_TOLERANCE_ULPS, round_decimals, round_half_away  # noqa: E402

PROBE_DIR = Path(__file__).resolve().parents[1] / "matlab_probes"


def _probe_files():
    return sorted(PROBE_DIR.glob("round_probe_*.json"))


def test_halves_go_away_from_zero():
    assert round_half_away(0.5) == 1
    assert round_half_away(-0.5) == -1
    assert round_half_away(2.5) == 3      # NumPy would give 2
    assert round_decimals(1.00005, 4) == pytest.approx(1.0001, abs=0)
    assert round_decimals(-1.00005, 4) == pytest.approx(-1.0001, abs=0)


def test_representation_error_is_compensated():
    """The documented MATLAB behavior: round(2.675, 2) is 2.68, not 2.67."""
    assert round_decimals(2.675, 2) == 2.68
    assert round_decimals(-2.675, 2) == -2.68


def test_non_finite_values_pass_through():
    values = round_decimals(np.array([np.nan, np.inf, -np.inf, 0.0]), 4)
    assert math.isnan(values[0])
    assert values[1] == np.inf
    assert values[2] == -np.inf
    assert values[3] == 0.0


def test_scalar_and_array_paths_agree():
    values = np.array([1668600000.1863499, 5177.31334, -0.00005, 268983.38314])
    from_array = round_decimals(values, 4)
    from_scalars = [round_decimals(float(value), 4) for value in values]
    assert list(from_array) == from_scalars


@pytest.mark.parametrize("probe_path", _probe_files() or [pytest.param(None, marks=pytest.mark.skip(
    reason="no MATLAB round probe recorded"))])
def test_matches_matlab_round(probe_path):
    probe = json.loads(probe_path.read_text())
    inputs = np.asarray(probe["input"], dtype=float)
    expected = np.asarray(probe["matlab_round"], dtype=float)

    ours = round_decimals(inputs, 4)
    mismatches = np.flatnonzero(ours != expected)
    assert mismatches.size == 0, (
        f"{mismatches.size}/{inputs.size} differ, e.g. input {inputs[mismatches[0]]!r}: "
        f"MATLAB {expected[mismatches[0]]!r}, mvtpy {ours[mismatches[0]]!r}")


@pytest.mark.parametrize("probe_path", _probe_files() or [pytest.param(None, marks=pytest.mark.skip(
    reason="no MATLAB round probe recorded"))])
def test_tie_tolerance_is_bracketed(probe_path):
    """A tolerance of half or twice one ULP does not reproduce MATLAB.

    This is what makes the fitted value meaningful rather than arbitrary: the
    probe set discriminates between them.
    """
    probe = json.loads(probe_path.read_text())
    inputs = np.asarray(probe["input"], dtype=float)
    expected = np.asarray(probe["matlab_round"], dtype=float)

    def round_with(tolerance_ulps):
        scaled = inputs * 1e4
        magnitude = np.abs(scaled)
        floor = np.floor(magnitude)
        fraction = magnitude - floor
        ulp = np.nextafter(magnitude, np.inf) - magnitude
        away = (fraction > 0.5) | (np.abs(fraction - 0.5) <= tolerance_ulps * ulp)
        return np.copysign(floor + away, scaled) / 1e4

    assert np.array_equal(round_with(TIE_TOLERANCE_ULPS), expected)
    assert not np.array_equal(round_with(0.5), expected)
    assert not np.array_equal(round_with(2.0), expected)
