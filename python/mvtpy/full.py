"""Build the "full" data set for one MOTION segment.

Python equivalent of ``Scripts/generate_data_mvt_full.m``. It is the ``slim``
stage plus three things:

* **both directions.** ``slim`` keeps westbound trajectories only; ``full``
  keeps every trajectory and records ``direction``. The lane machinery already
  handles eastbound - it estimates a driving line per direction and applies a
  different scale/offset - so this stage estimates both lines up front.
* **a reference trajectory.** A two-phase constant-acceleration drive covering
  the same distance in the same time, used as the counterfactual "what would an
  undisturbed vehicle have burned?". See :func:`reference_trajectory`.
* **flat-road and reference fuel.** The same fuel model is evaluated four times
  per trajectory: as-driven, as-driven on zero grade, reference, and reference
  on zero grade. All four totals use the deterministic quadrature.

The three helpers the MATLAB stages share - ``assign_lanes``,
``clip_lane_changes`` and ``calculate_distance_to_avs`` - are byte-identical
between ``generate_data_mvt_slim.m`` and ``generate_data_mvt_full.m``, so this
module reuses :mod:`mvtpy.lanes` and :mod:`mvtpy.avdist` unchanged.

Field order follows the 41-field ``init_data_struct`` of the MATLAB source
followed by the four upstream fields in assignment order, verified against a
released ``full`` file. :data:`FIELD_ORDER` records it.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

from . import avdist, fuel, kinematics as kin, lanes
from .rawio import iter_trajectories
from .slim import write_segment as _write_segment

__all__ = ["FIELD_ORDER", "reference_trajectory", "build_record", "build_segment",
           "write_segment"]

#: Output key order: init_data_struct (41 fields) then the four upstream fields
#: in the order generate_data_mvt_full.m assigns them. Verified by scanning the
#: first record of a released full segment.
FIELD_ORDER = (
    "trajectory_id",
    "timestamp",
    "x_position_meters",
    "y_position_corrected_meters",
    "coarse_vehicle_class",
    "direction",
    "first_timestamp",
    "last_timestamp",
    "starting_x",
    "ending_x",
    "length",
    "width",
    "height",
    "total_distance_traversed_meters",
    "speed_meters_per_second",
    "acceleration_meters_per_second_per_second",
    "road_grade_radians",
    "lane_number",
    "reference_a1_meters_per_second_per_second",
    "reference_a2_meters_per_second_per_second",
    "energy_model",
    "fuel_rate_grams_per_second",
    "percent_infeasibility",
    "total_fuel_consumed_grams",
    "total_fuel_consumed_gallons",
    "total_fuel_economy_mpg",
    "fuel_rate_flat_road_grams_per_second",
    "percent_infeasibility_flat_road",
    "total_fuel_consumed_flat_road_grams",
    "total_fuel_consumed_flat_road_gallons",
    "total_fuel_economy_flat_road_mpg",
    "reference_fuel_rate_grams_per_second",
    "percent_reference_infeasibility",
    "total_reference_fuel_consumed_grams",
    "reference_fuel_rate_flat_road_grams_per_second",
    "percent_reference_infeasibility_flat_road",
    "total_reference_fuel_consumed_flat_road_grams",
    "downstream_av_id",
    "distance_to_downstream_av_meters",
    "downstream_engaged_av_id",
    "distance_to_downstream_engaged_av_meters",
    "upstream_engaged_av_id",
    "distance_to_upstream_engaged_av_meters",
    "upstream_av_id",
    "distance_to_upstream_av_meters",
)


def reference_trajectory(distance: np.ndarray, speed: np.ndarray,
                         time: np.ndarray) -> Tuple[float, float, np.ndarray,
                                                    np.ndarray, np.ndarray]:
    """The two-phase reference drive: ``(a1, a2, v_ref, a_ref, x_ref)``.

    Mirrors ``generate_data_mvt_full.m`` lines 232-251. The primary form
    accelerates at ``a1`` until the midpoint ``t1 = t(end)/2`` then at ``a2``,
    chosen to cover the same distance in the same time as the real drive.

    That form can imply a negative speed on a trajectory that nearly stops. When
    it does, MATLAB switches to a three-phase form that decelerates to rest,
    waits, then accelerates away. **The switch is a branch on an exact
    floating-point comparison** (``any(vRef < 0)``), so a trajectory sitting an
    ULP either side of zero takes a different path and every reference field
    changes together rather than in the last decimal. Nothing here can soften
    that; it is inherent to the MATLAB formulation and is noted in
    docs/PYTHON_PORT.md.

    ``distance``, ``speed`` and ``time`` are all relative to the trajectory
    start, exactly as the MATLAB code has them at this point.
    """
    x, v, t = distance, speed, time
    t_end = t[-1]
    # a2 as a function of the switch time, then a1 given a2. Written to match
    # the MATLAB expression term for term rather than simplified: the algebra is
    # equivalent but the floating-point result of a rearranged form is not.
    def a2_of(t1: float) -> float:
        return ((x[-1] - x[0] - 0.5 * t1 * v[0] - v[-1] * t_end + 0.5 * t1 * v[-1])
                / (-0.5 * t_end ** 2 + 0.5 * t1 * t_end))

    def a1_of(t1: float) -> float:
        return (v[-1] + a2_of(t1) * (t1 - t_end) - v[0]) / t1

    t1 = t_end / 2
    a1, a2 = a1_of(t1), a2_of(t1)

    early = t <= t1
    late = t > t1
    v_ref = (v[0] + a1 * t) * early + (v[-1] + a2 * (t - t_end)) * late
    a_ref = a1 * early + a2 * late
    x_ref = ((x[0] + v[0] * t + a1 * t ** 2 / 2) * early
             + ((x[-1] - v[-1] * t_end + a2 * t_end ** 2 / 2)
                + v[-1] * t + a2 * (t / 2 - t_end) * t) * late)

    if np.any(v_ref < 0):
        span = x[-1] - x[0]
        mean_speed = (v[0] + v[-1]) / 2
        tt1 = span / mean_speed
        tt2 = t_end - tt1
        a1 = -mean_speed / span * v[0]
        a2 = mean_speed / span * v[-1]

        first = t <= tt1
        middle = (t > tt1) & (t < tt2)
        last = t > tt2
        v_ref = (v[0] + a1 * t) * first + 0.0 * middle + (v[-1] + a2 * (t - t_end)) * last
        a_ref = a1 * first + 0.0 * middle + a2 * last
        x_ref = ((x[0] + v[0] * t + a1 * t ** 2 / 2) * first
                 + (x[0] + v[0] * tt1 + a1 * tt1 ** 2 / 2) * middle
                 + ((x[-1] - v[-1] * t_end + a2 * t_end ** 2 / 2)
                    + v[-1] * t + a2 * (t / 2 - t_end) * t) * last)

    return a1, a2, v_ref, a_ref, x_ref


def build_record(segment: dict, grade_map: kin.GradeMap,
                 av_runs: Sequence[avdist.AvRun]) -> Dict[str, object]:
    """Assemble one full-data-set trajectory record from a clipped segment."""
    time = np.asarray(segment["timestamp"], dtype=float)
    x_feet = np.asarray(segment["x_position"], dtype=float)
    y_feet = np.asarray(segment["y_position"], dtype=float)
    direction = segment["direction"]

    x_meters = kin.FT_TO_METER * (x_feet - kin.ORIGIN_X_FEET)
    distance = kin.FT_TO_METER * np.abs(x_feet - x_feet[0])
    relative_time = time - time[0]

    speed = kin.speed(distance, relative_time)
    acceleration = kin.acceleration(distance, relative_time)
    grade = grade_map(x_meters, direction)

    a1, a2, speed_ref, accel_ref, x_ref = reference_trajectory(
        distance, speed, relative_time)
    # The reference runs forward from the same starting point; eastbound x
    # increases, westbound decreases.
    x_ref_absolute = (x_meters[0] + x_ref) if direction > 0 else (x_meters[0] - x_ref)
    grade_ref = grade_map(x_ref_absolute, direction)

    model = fuel.model_for_coarse_class(segment["coarse_vehicle_class"])
    flat = np.zeros_like(grade)
    flat_ref = np.zeros_like(grade_ref)

    rate, _, infeasible = model(speed, acceleration, grade, project=True)
    rate_flat, _, infeasible_flat = model(speed, acceleration, flat, project=True)
    rate_ref, _, infeasible_ref = model(speed_ref, accel_ref, grade_ref, project=True)
    rate_ref_flat, _, infeasible_ref_flat = model(
        speed_ref, accel_ref, flat_ref, project=True)

    total_grams = kin.trapezoid_integral(relative_time, rate)
    total_flat_grams = kin.trapezoid_integral(relative_time, rate_flat)
    total_ref_grams = kin.trapezoid_integral(relative_time, rate_ref)
    total_ref_flat_grams = kin.trapezoid_integral(relative_time, rate_ref_flat)

    total_gallons = kin.GRAM_TO_GALLON * total_grams
    total_flat_gallons = kin.GRAM_TO_GALLON * total_flat_grams
    travelled = abs(distance[-1] - distance[0])
    with np.errstate(divide="ignore", invalid="ignore"):
        economy = (travelled * kin.METER_TO_MILE) / total_gallons
        economy_flat = (travelled * kin.METER_TO_MILE) / total_flat_gallons

    distances = avdist.distance_to_avs(segment, av_runs)

    return {
        "trajectory_id": {"x_oid": segment["trajectory_id"]},
        "timestamp": time,
        "x_position_meters": x_meters,
        "y_position_corrected_meters": y_feet * kin.FT_TO_METER,
        "coarse_vehicle_class": segment["coarse_vehicle_class"],
        "direction": direction,
        "first_timestamp": segment["first_timestamp"],
        "last_timestamp": segment["last_timestamp"],
        "starting_x": x_meters[0],
        "ending_x": x_meters[-1],
        "length": segment["length"],
        "width": segment["width"],
        "height": segment["height"],
        "total_distance_traversed_meters": travelled,
        "speed_meters_per_second": speed,
        "acceleration_meters_per_second_per_second": acceleration,
        "road_grade_radians": grade,
        "lane_number": segment["lane"],
        "reference_a1_meters_per_second_per_second": a1,
        "reference_a2_meters_per_second_per_second": a2,
        "energy_model": model.name,
        "fuel_rate_grams_per_second": rate,
        "percent_infeasibility": _percent(infeasible),
        "total_fuel_consumed_grams": total_grams,
        "total_fuel_consumed_gallons": total_gallons,
        "total_fuel_economy_mpg": economy,
        "fuel_rate_flat_road_grams_per_second": rate_flat,
        "percent_infeasibility_flat_road": _percent(infeasible_flat),
        "total_fuel_consumed_flat_road_grams": total_flat_grams,
        "total_fuel_consumed_flat_road_gallons": total_flat_gallons,
        "total_fuel_economy_flat_road_mpg": economy_flat,
        "reference_fuel_rate_grams_per_second": rate_ref,
        "percent_reference_infeasibility": _percent(infeasible_ref),
        "total_reference_fuel_consumed_grams": total_ref_grams,
        "reference_fuel_rate_flat_road_grams_per_second": rate_ref_flat,
        "percent_reference_infeasibility_flat_road": _percent(infeasible_ref_flat),
        "total_reference_fuel_consumed_flat_road_grams": total_ref_flat_grams,
        # Downstream distances are positive, upstream negative.
        "downstream_av_id": distances["AvIdDS"],
        "distance_to_downstream_av_meters": distances["distanceToAvDS"],
        "downstream_engaged_av_id": distances["AvIdDSEng"],
        "distance_to_downstream_engaged_av_meters": distances["distanceToAvDSEng"],
        "upstream_engaged_av_id": distances["AvIdUSEng"],
        "distance_to_upstream_engaged_av_meters": _negate(distances["distanceToAvUSEng"]),
        "upstream_av_id": distances["AvIdUS"],
        "distance_to_upstream_av_meters": _negate(distances["distanceToAvUS"]),
    }


def build_segment(raw_path, gps_path, grade_csv, limit: Optional[int] = None) -> List[dict]:
    """Process one raw segment into full-data-set records, in MATLAB's order."""
    grade_map = kin.GradeMap.from_csv(grade_csv)
    av_runs = avdist.load_av_runs(gps_path)

    # A driving line per direction, each estimated from that direction's
    # trajectories only - as MATLAB's assign_lanes does with its separate
    # xWest/yWest and xEast/yEast samples. Estimating both needs two passes over
    # the file, which the streaming reader supports.
    driving_lines = {}
    for heading in (-1, 1):
        try:
            driving_lines[heading] = lanes.estimate_driving_line(
                iter_trajectories(raw_path), direction=heading)
        except ValueError:
            # A segment with no long trajectory in one direction: leave it out
            # and skip those records below rather than failing the whole file.
            driving_lines[heading] = None

    records: List[dict] = []
    for record in iter_trajectories(raw_path):
        if limit is not None and len(records) >= limit:
            break
        driving_line = driving_lines[-1 if record["direction"] < 0 else 1]
        if driving_line is None:
            continue
        y_corr, lane = lanes.assign_lanes(record, driving_line)
        for segment in lanes.clip_lane_changes(record, y_corr, lane):
            records.append(build_record(segment, grade_map, av_runs))
    return records


def write_segment(records: Iterable[dict], output_path) -> Path:
    """Round, encode, and write records exactly as the MATLAB stage does."""
    return _write_segment(records, output_path, field_order=FIELD_ORDER)


# ---------------------------------------------------------------------------


def _percent(infeasible: np.ndarray) -> float:
    """Share of samples flagged infeasible, as ``length(x(x>0))/length(x)*100``."""
    return float(np.count_nonzero(infeasible > 0)) / infeasible.size * 100


def _negate(values):
    return None if values is None else -values
