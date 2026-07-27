"""Trajectory kinematics and road grade, ported from ``generate_data_mvt_slim.m``.

Everything here reproduces a specific expression in the MATLAB source, in the
same operation order, so that results agree to the last bit wherever the inputs
do. Each function names the lines it mirrors.
"""

from __future__ import annotations

from typing import Optional

import numpy as np

__all__ = [
    "FT_TO_METER",
    "METER_TO_MILE",
    "GRAM_TO_GALLON",
    "ORIGIN_X_FEET",
    "MILL_CREEK_OFFSET_MILES",
    "speed",
    "acceleration",
    "road_grade",
    "trapezoid_integral",
    "GradeMap",
]

#: [m/ft] conversion factor from feet to meter
FT_TO_METER = 0.3048
#: [mile/m] conversion factor from meter to mile
METER_TO_MILE = 6.213712e-04
#: [gallon/g] conversion factor from fuel gram to gallon
GRAM_TO_GALLON = 3.522294e-04
#: [ft] location of the origin for the x-coordinates (Mill Creek Bridge)
ORIGIN_X_FEET = 309804.0625
#: [mile] distance between the Mill Creek origin (MM58.675) and MM58.9, the
#: estimated origin of the road grade map
MILL_CREEK_OFFSET_MILES = 0.225


def speed(x: np.ndarray, t: np.ndarray) -> np.ndarray:
    """Central-difference speed with duplicated end points.

    Mirrors::

        v = (x([2:end,end])-x([1,1:end-1]))./(t([2:end,end])-t([1,1:end-1]))
    """
    x = np.asarray(x, dtype=float)
    t = np.asarray(t, dtype=float)
    forward = np.concatenate((x[1:], x[-1:]))
    backward = np.concatenate((x[:1], x[:-1]))
    dt = np.concatenate((t[1:], t[-1:])) - np.concatenate((t[:1], t[:-1]))
    return (forward - backward) / dt


def acceleration(x: np.ndarray, t: np.ndarray) -> np.ndarray:
    """Second central difference, with the first/last values repeated.

    Mirrors::

        a = (x(1:end-2)-2*x(2:end-1)+x(3:end))./((t(3:end)-t(1:end-2))/2).^2;
        a = a([1,1:end,end]);
    """
    x = np.asarray(x, dtype=float)
    t = np.asarray(t, dtype=float)
    interior = (x[:-2] - 2 * x[1:-1] + x[2:]) / ((t[2:] - t[:-2]) / 2) ** 2
    return np.concatenate((interior[:1], interior, interior[-1:]))


class GradeMap:
    """Piecewise-linear road grade fit from ``Models/Eastbound_grade_fit.csv``.

    Columns 2 and 3 are the start and end of each fitted cell (in miles from the
    grade-map origin); columns 4 and 5 are the slope and intercept, in percent.
    """

    def __init__(self, grade_data: np.ndarray):
        grade_data = np.asarray(grade_data, dtype=float)
        starts = grade_data[:, 1]
        ends = grade_data[:, 2]
        self.points = np.concatenate((starts, ends[-1:]))
        self.slope = grade_data[:, 3]
        self.intercept = grade_data[:, 4]

    @classmethod
    def from_csv(cls, path) -> "GradeMap":
        """Read the grade fit CSV, skipping its header as ``readmatrix`` does.

        Columns: interval_number, interval_start, interval_end, slope, intercept.
        """
        return cls(np.loadtxt(path, delimiter=",", skiprows=1))

    def __call__(self, x_meters: np.ndarray, direction: float) -> np.ndarray:
        return road_grade(x_meters, direction, self.points, self.slope, self.intercept)


def road_grade(x_meters: np.ndarray, direction: float, points: np.ndarray,
               slope: np.ndarray, intercept: np.ndarray) -> np.ndarray:
    """Road grade in radians at each position.

    Mirrors the block that starts ``xNew = ...*meter2mileFactor - mcDist`` in
    generate_data_mvt_slim.m: locate each position in the fitted cells, clamp to
    the mapped range, evaluate the percent-grade line, and take the arcsine.
    Westbound trajectories get the negated grade.
    """
    x_meters = np.asarray(x_meters, dtype=float)
    x_new = x_meters * METER_TO_MILE - MILL_CREEK_OFFSET_MILES

    # cellInd: index of the last cell start at or below x_new (MATLAB loops over
    # cells assigning j wherever x_new - points(j) >= 0).
    cell_index = np.searchsorted(points, x_new, side="right")
    cell_index = np.clip(cell_index, 1, len(points) - 1)

    x_local = np.clip(x_new, points[0], points[-1])
    grade_percent = (slope[cell_index - 1] * x_local / 100
                     + intercept[cell_index - 1] / 100)
    theta = np.arcsin(grade_percent)
    return theta if direction > 0 else -theta


#: Match ``flag_deterministic_quadrature`` in generate_data_mvt_{slim,full}.m.
#: True uses compensated summation, which is bit-identical on every platform and
#: in both languages. False reproduces MATLAB's ``dot`` (a BLAS call) as closely
#: as numpy can, for comparison against pre-2026-07 outputs - but note the two
#: BLAS libraries do not agree with each other either, so "false" is not a
#: well-defined target. See docs/REPRODUCIBLE_QUADRATURE.md.
DETERMINISTIC_QUADRATURE = True


def neumaier_dot(a: np.ndarray, b: np.ndarray) -> float:
    """Sum of ``a*b`` by compensated (Kahan-Babuska-Neumaier) summation.

    A fixed sequence of IEEE-754 double operations, so it returns identical bits
    everywhere. Deliberately a scalar Python loop: the recurrence is sequential,
    and any vectorized reassociation would reintroduce exactly the
    order-dependence this exists to remove.

    Verified against the same loop in MATLAB (``mvt.neumaierDot``) on 60 real
    trajectories: 60 of 60 bit-identical, and equal to the exactly-rounded sum
    on all 60.
    """
    total = 0.0
    compensation = 0.0
    for x, y in zip(np.asarray(a, dtype=float).tolist(),
                    np.asarray(b, dtype=float).tolist()):
        product = x * y
        running = total + product
        if abs(total) >= abs(product):
            compensation += (total - running) + product
        else:
            compensation += (product - running) + total
        total = running
    return total + compensation


def trapezoid_integral(t: np.ndarray, values: np.ndarray,
                       deterministic: Optional[bool] = None) -> float:
    """Trapezoidal quadrature.

    Mirrors ``integrate = @(t,v) <dot>(t(2:end)-t(1:end-1), (v(1:end-1)+v(2:end))/2)``
    where ``<dot>`` is ``mvt.neumaierDot`` or MATLAB's ``dot``, selected by
    ``flag_deterministic_quadrature`` there and :data:`DETERMINISTIC_QUADRATURE`
    here. The two flags must agree or the fuel totals will differ in the 4th
    decimal on a small fraction of trajectories.
    """
    t = np.asarray(t, dtype=float)
    values = np.asarray(values, dtype=float)
    weights = t[1:] - t[:-1]
    heights = (values[:-1] + values[1:]) / 2
    if DETERMINISTIC_QUADRATURE if deterministic is None else deterministic:
        return neumaier_dot(weights, heights)
    return float(np.dot(weights, heights))
