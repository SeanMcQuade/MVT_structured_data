"""Parity tests for the trajectory math ported from generate_data_mvt_slim.m.

Each test runs the Python implementation on the raw I-24 MOTION samples and
requires the result, rounded to four decimals the way the pipeline does, to
equal exactly what MATLAB wrote into the released data.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from mvtpy import kinematics as kin  # noqa: E402
from mvtpy.matround import round_decimals  # noqa: E402


def _relative_distance(x_feet: np.ndarray) -> np.ndarray:
    """``x = ft2meterFactor*abs(x-x(1))`` - distance travelled, always positive."""
    return kin.FT_TO_METER * np.abs(x_feet - x_feet[0])


def test_position_conversion_matches_released(segment_pairs):
    for pair in segment_pairs:
        ours = round_decimals(kin.FT_TO_METER * (pair["x_feet"] - kin.ORIGIN_X_FEET), 4)
        assert np.array_equal(ours, np.asarray(pair["segment"]["x_position_meters"]))


def test_speed_matches_released(segment_pairs):
    for pair in segment_pairs:
        time = pair["time"] - pair["time"][0]
        ours = round_decimals(kin.speed(_relative_distance(pair["x_feet"]), time), 4)
        assert np.array_equal(ours, np.asarray(pair["segment"]["speed_meters_per_second"]))


def test_acceleration_matches_released(segment_pairs):
    for pair in segment_pairs:
        time = pair["time"] - pair["time"][0]
        ours = round_decimals(kin.acceleration(_relative_distance(pair["x_feet"]), time), 4)
        assert np.array_equal(
            ours,
            np.asarray(pair["segment"]["acceleration_meters_per_second_per_second"]),
        )


def test_road_grade_matches_released(segment_pairs, grade_map):
    for pair in segment_pairs:
        x_meters = kin.FT_TO_METER * (pair["x_feet"] - kin.ORIGIN_X_FEET)
        ours = round_decimals(grade_map(x_meters, pair["direction"]), 4)
        assert np.array_equal(ours, np.asarray(pair["segment"]["road_grade_radians"]))


def test_grade_map_shape():
    """The fit has one slope/intercept per interval, and one more edge than cells."""
    grade = kin.GradeMap.from_csv(
        Path(__file__).resolve().parents[2] / "Models" / "Eastbound_grade_fit.csv"
    )
    assert len(grade.points) == len(grade.slope) + 1
    assert len(grade.slope) == len(grade.intercept)
    assert grade.points[0] == pytest.approx(0.0)


def test_trapezoid_integral_matches_matlab_dot_form():
    time = np.array([0.0, 0.5, 1.5, 3.0])
    values = np.array([1.0, 2.0, 0.0, 4.0])
    expected = float(np.dot(np.diff(time), (values[:-1] + values[1:]) / 2))
    assert kin.trapezoid_integral(time, values) == expected


def test_neumaier_dot_matches_matlab_bit_for_bit():
    """Compensated summation is a fixed sequence of IEEE ops, so MATLAB's
    mvt.neumaierDot and this loop agree exactly.

    Verified on 60 real trajectories (60/60 bit-identical). The case below is a
    regression guard: the naive product-then-sum disagrees with it, which is the
    platform-dependence the flag exists to remove.
    """
    from mvtpy.kinematics import neumaier_dot

    # Terms spanning many magnitudes, where accumulation order matters.
    a = np.array([1e16, 1.0, -1e16, 1.0])
    b = np.ones(4)
    assert neumaier_dot(a, b) == 2.0          # exact
    assert float(np.sum(a * b)) != 2.0        # naive summation loses both ones


def test_quadrature_flag_selects_the_implementation():
    from mvtpy import kinematics

    t = np.arange(500) * 0.04
    v = np.abs(np.sin(t)) * 2.0
    deterministic = kinematics.trapezoid_integral(t, v, deterministic=True)
    blas = kinematics.trapezoid_integral(t, v, deterministic=False)

    assert deterministic == kinematics.neumaier_dot(t[1:] - t[:-1],
                                                    (v[:-1] + v[1:]) / 2)
    # Same value to within rounding; the point is that only one of them is
    # reproducible across platforms.
    assert deterministic == pytest.approx(blas, rel=1e-12)


def test_default_flag_is_deterministic():
    """Must stay in step with flag_deterministic_quadrature in the .m files."""
    from mvtpy import kinematics

    assert kinematics.DETERMINISTIC_QUADRATURE is True
