"""Build the released "slim" data set for one MOTION segment.

Python equivalent of the body of ``Scripts/generate_data_mvt_slim.m``: it takes
a raw 10-minute I-24 MOTION segment, keeps westbound trajectories, assigns
lanes and clips lane changes, computes kinematics, road grade, fuel, and
distances to the CIRCLES control vehicles, and writes the result as JSON that
is byte-identical to MATLAB's.

Field order matters and is not alphabetical: it follows the MATLAB struct
declared in ``init_data_struct`` (27 fields), after which the four upstream
fields appear in the order the code first assigns them. :data:`FIELD_ORDER`
records that, and the writer emits keys in exactly that sequence.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence

import numpy as np

from . import avdist, fuel, kinematics as kin, lanes, matjson
from .matround import round_decimals
from .rawio import iter_trajectories

__all__ = ["FIELD_ORDER", "build_record", "build_segment", "write_segment"]

#: Output key order, from init_data_struct plus the dynamically added upstream
#: fields. Verified against the released files.
FIELD_ORDER = (
    "trajectory_id",
    "timestamp",
    "x_position_meters",
    "y_position_corrected_meters",
    "coarse_vehicle_class",
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
    "energy_model",
    "fuel_rate_grams_per_second",
    "percent_infeasibility",
    "total_fuel_consumed_grams",
    "total_fuel_consumed_gallons",
    "total_fuel_economy_mpg",
    "downstream_av_id",
    "distance_to_downstream_av_meters",
    "downstream_engaged_av_id",
    "distance_to_downstream_engaged_av_meters",
    "upstream_engaged_av_id",
    "distance_to_upstream_engaged_av_meters",
    "upstream_av_id",
    "distance_to_upstream_av_meters",
)

#: Everything except the identifier is rounded to four decimals before
#: encoding, matching the loop over fields 2..end in the MATLAB source.
DECIMALS = 4


def build_record(segment: dict, grade_map: kin.GradeMap,
                 av_runs: Sequence[avdist.AvRun]) -> Dict[str, object]:
    """Assemble one released trajectory record from a clipped segment."""
    time = np.asarray(segment["timestamp"], dtype=float)
    x_feet = np.asarray(segment["x_position"], dtype=float)
    y_feet = np.asarray(segment["y_position"], dtype=float)

    x_meters = kin.FT_TO_METER * (x_feet - kin.ORIGIN_X_FEET)
    # Distance travelled, measured from the segment start (always increasing).
    distance = kin.FT_TO_METER * np.abs(x_feet - x_feet[0])
    relative_time = time - time[0]

    speed = kin.speed(distance, relative_time)
    acceleration = kin.acceleration(distance, relative_time)
    grade = grade_map(x_meters, segment["direction"])

    model = fuel.model_for_coarse_class(segment["coarse_vehicle_class"])
    rate, _, infeasible = model(speed, acceleration, grade, project=True)

    total_grams = kin.trapezoid_integral(relative_time, rate)
    total_gallons = kin.GRAM_TO_GALLON * total_grams
    travelled = abs(distance[-1] - distance[0])
    with np.errstate(divide="ignore", invalid="ignore"):
        economy = (travelled * kin.METER_TO_MILE) / total_gallons

    distances = avdist.distance_to_avs(segment, av_runs)

    return {
        "trajectory_id": {"x_oid": segment["trajectory_id"]},
        "timestamp": time,
        "x_position_meters": x_meters,
        "y_position_corrected_meters": y_feet * kin.FT_TO_METER,
        "coarse_vehicle_class": segment["coarse_vehicle_class"],
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
        "energy_model": model.name,
        "fuel_rate_grams_per_second": rate,
        "percent_infeasibility": float(np.count_nonzero(infeasible > 0)) / infeasible.size * 100,
        "total_fuel_consumed_grams": total_grams,
        "total_fuel_consumed_gallons": total_gallons,
        "total_fuel_economy_mpg": economy,
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
    """Process one raw segment into released records, in MATLAB's order."""
    grade_map = kin.GradeMap.from_csv(grade_csv)
    av_runs = avdist.load_av_runs(gps_path)

    # The driving line is a whole-file statistic over westbound trajectories.
    driving_line = lanes.estimate_driving_line(
        record for record in iter_trajectories(raw_path) if record["direction"] < 0)

    records: List[dict] = []
    for count, record in enumerate(iter_trajectories(raw_path)):
        if record["direction"] >= 0:
            continue
        if limit is not None and len(records) >= limit:
            break
        y_corr, lane = lanes.assign_lanes(record, driving_line)
        for segment in lanes.clip_lane_changes(record, y_corr, lane):
            records.append(build_record(segment, grade_map, av_runs))
    return records


def write_segment(records: Iterable[dict], output_path) -> Path:
    """Round, encode, and write records exactly as the MATLAB stage does."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    payload = matjson.dumps([_prepare(record) for record in records])
    # Write to a temporary name and rename, mirroring mvt.atomicWrite.
    temporary = output_path.with_name(f".{output_path.name}.tmp")
    temporary.write_text(payload, encoding="utf-8")
    temporary.replace(output_path)
    return output_path


# ---------------------------------------------------------------------------


def _negate(values):
    return None if values is None else -values


def _prepare(record: Dict[str, object]) -> Dict[str, object]:
    """Round numeric fields and convert to JSON-ready values, in field order."""
    prepared: Dict[str, object] = {}
    for name in FIELD_ORDER:
        value = record[name]
        if name == "trajectory_id" or isinstance(value, (str, dict)):
            prepared[name] = value
        else:
            prepared[name] = _to_json(round_decimals(value, DECIMALS)
                                      if value is not None else None)
    return prepared


def _to_json(value):
    """Mirror MATLAB's array-to-JSON mapping.

    ``jsonencode`` writes a 1x1 array as a bare number, an empty array as ``[]``,
    and anything longer as an array.
    """
    if value is None:
        return []
    if np.isscalar(value):
        return float(value)
    array = np.asarray(value, dtype=float)
    if array.size == 0:
        return []
    if array.size == 1:
        return float(array.reshape(-1)[0])
    return array.tolist()
