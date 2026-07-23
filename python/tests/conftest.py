"""Shared fixtures: real raw/released data pairs for parity testing.

The MATLAB pipeline is the reference implementation, and the data it already
produced is the oracle. These fixtures pair each released "slim" trajectory
segment with the window of raw I-24 MOTION samples it was computed from, so a
test can run the Python implementation on the raw inputs and demand the exact
values MATLAB wrote.

Everything skips cleanly when the data trees are not present (the repository is
a sibling of ``data/`` and ``results/``, which are distributed separately).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from mvtpy import matjson  # noqa: E402
from mvtpy.matround import round_decimals  # noqa: E402
from tools.json_parity import extract_prefix_records  # noqa: E402

WORKSPACE = Path(__file__).resolve().parents[3]
RAW_SEGMENT = (WORKSPACE / "data" / "i24motion" / "2022-11-16"
               / "64888dc1d2834b01ce52a162__wed_0_00.json")
SLIM_SEGMENT = (WORKSPACE / "results" / "slim" / "2022-11-16"
                / "I-24MOTION_2022-11-16_05-59-59.json")
GRADE_FIT = Path(__file__).resolve().parents[2] / "Models" / "Eastbound_grade_fit.csv"

#: How much of each (multi-GB) file to read, and how many records to decode.
READ_BYTES = 60 * 1024 * 1024
RAW_RECORDS = 400
SLIM_RECORDS = 120


def _read_records(path: Path, count: int):
    with path.open("rb") as handle:
        text = handle.read(READ_BYTES).decode("utf-8", errors="ignore")
    return matjson.loads(extract_prefix_records(text, count))


#: Trajectories to run through lane assignment in the parity fixtures. The
#: driving line always uses the whole file (it is a whole-file statistic); only
#: the per-trajectory stage is limited, to keep the suite quick.
LANE_TRAJECTORIES = 200


@pytest.fixture(scope="session")
def raw_segment_path():
    if not RAW_SEGMENT.is_file():
        pytest.skip("raw data/ not available")
    return RAW_SEGMENT


@pytest.fixture(scope="session")
def driving_line(raw_segment_path):
    """Stage-1 driving line, estimated by streaming the entire raw segment."""
    from mvtpy import lanes
    from mvtpy.rawio import iter_trajectories

    westbound = (record for record in iter_trajectories(raw_segment_path)
                 if record["direction"] < 0)
    return lanes.estimate_driving_line(westbound)


CARS_DIR = WORKSPACE / "data" / "cars"


@pytest.fixture(scope="session")
def gps_assembly_case():
    """Preprocessed runs, connection status, and released records for day 16.

    Pairs each resampled control-vehicle run with the released GPS record it
    produced (aligned on the 10 Hz grid), for field-by-field parity checks.
    Skips cleanly without the raw vehicle CSVs or the released GPS file.
    """
    pytest.importorskip("pandas")
    cars_gps = CARS_DIR / "cars_gps"
    vins = CARS_DIR / "cars_vins.csv"
    pings = CARS_DIR / "veh_ping_20221116.csv"
    released_file = WORKSPACE / "results" / "gps" / "CIRCLES_GPS_10Hz_2022-11-16.json"
    if not (cars_gps.is_dir() and vins.is_file() and pings.is_file()
            and released_file.is_file()):
        pytest.skip("raw vehicle GPS data or released GPS file not available")

    from mvtpy import gpsassemble as ga, gpsruns
    from mvtpy.matround import round_decimals
    from mvtpy.rawio import iter_trajectories

    runs = gpsruns.parse_gps_data(cars_gps, 16)
    lane_map = ga.load_lane_map(vins)
    preprocessed = [ga.preprocess_run(run, index + 1, lane_map)
                    for index, run in enumerate(runs)]
    status = ga.connection_status(pings, vins)
    released = list(iter_trajectories(released_file))

    def released_for(pre):
        best, best_overlap = None, 0.0
        for record in released:
            if int(record["av_id"]) != pre.vin:
                continue
            overlap = (min(pre.timestamp[-1], record["timestamp"][-1])
                       - max(pre.timestamp[0], record["timestamp"][0]))
            if overlap > best_overlap:
                best, best_overlap = record, overlap
        return best

    pairs = []
    for pre in preprocessed:
        record = released_for(pre)
        if record is None:
            continue
        grid = round_decimals(pre.timestamp, 6)
        released_time = np.asarray(record["timestamp"])
        start = int(np.searchsorted(grid, released_time[0] - 1e-7))
        if start + len(released_time) > len(grid):
            continue
        window = slice(start, start + len(released_time))
        if not np.allclose(grid[window], released_time, atol=1e-6):
            continue
        pairs.append({"pre": pre, "record": record, "window": window,
                      "status": status[pre.vin]})

    if not pairs:
        pytest.skip("no GPS runs could be aligned to released records")
    return pairs


MOTION_2022_11_16 = WORKSPACE / "data" / "i24motion" / "2022-11-16"
LANE_FIXTURE = (Path(__file__).resolve().parent / "fixtures"
                / "motion_lanes_2022-11-16_seg08.npz")


@pytest.fixture(scope="session")
def motion_segment_lanes():
    """(python lanes, MATLAB lanes) for one MOTION segment's trajectories.

    MATLAB reference is committed (tests/fixtures); the Python side is computed
    from the raw segment, so it skips when the raw MOTION data is not present.
    """
    if not LANE_FIXTURE.is_file():
        pytest.skip("committed lane reference missing")
    reference = np.load(LANE_FIXTURE, allow_pickle=True)
    segment = MOTION_2022_11_16 / str(reference["segment"])
    if not segment.is_file():
        pytest.skip("raw MOTION segment not available")

    from mvtpy.gpsmatch import assign_lanes_bidirectional
    from mvtpy.rawio import iter_trajectories

    trajectories = list(iter_trajectories(segment))
    ours = assign_lanes_bidirectional(trajectories)
    n = int(reference["n"])
    matlab = [reference[f"lane_{k}"] for k in range(n)]
    return ours[:n], matlab


@pytest.fixture(scope="session")
def matching_bias_case():
    """(computed bias, recovered target) per run for 2022-11-16.

    Skips unless both the raw data and a current-code GPS reference are present.
    The recovered target is median(ft2m*x_python - x_matlab) per run, which is
    exactly what the matching pass computes. Running the matching is slow
    (streams the day's MOTION segments), so this is opt-in via the data.
    """
    cars = WORKSPACE / "data" / "cars"
    reference = WORKSPACE / "results" / "gps" / "CIRCLES_GPS_10Hz_2022-11-16.json"
    if not ((cars / "cars_gps").is_dir() and MOTION_2022_11_16.is_dir()
            and reference.is_file()):
        pytest.skip("raw data or GPS reference not available")

    from mvtpy import gpsassemble as ga, gpsmatch, gpsruns
    from mvtpy.kinematics import FT_TO_METER
    from mvtpy.matround import round_decimals
    from mvtpy.rawio import iter_trajectories

    runs = gpsruns.parse_gps_data(cars / "cars_gps", 16)
    lane_map = ga.load_lane_map(cars / "cars_vins.csv")
    preprocessed = [ga.preprocess_run(run, index + 1, lane_map)
                    for index, run in enumerate(runs)]

    released = list(iter_trajectories(reference))
    target = {}
    for pre in preprocessed:
        best, best_overlap = None, 0.0
        for record in released:
            if int(record["av_id"]) != pre.vin:
                continue
            overlap = (min(pre.timestamp[-1], record["timestamp"][-1])
                       - max(pre.timestamp[0], record["timestamp"][0]))
            if overlap > best_overlap:
                best, best_overlap = record, overlap
        if best is None:
            continue
        grid = round_decimals(pre.timestamp, 6)
        released_time = np.asarray(best["timestamp"])
        start = int(np.searchsorted(grid, released_time[0] - 1e-7))
        if start + len(released_time) > len(grid):
            continue
        window = slice(start, start + len(released_time))
        if not np.allclose(grid[window], released_time, atol=1e-6):
            continue
        our_x = round_decimals(FT_TO_METER * pre.x_position, 6)[window]
        target[pre.index] = float(np.median(our_x - np.asarray(best["x_position"])))

    bias = gpsmatch.matching_bias(preprocessed, MOTION_2022_11_16, 16)
    return bias, target


@pytest.fixture(scope="session")
def av_runs():
    """Control-vehicle runs from the assembled GPS file for the same day."""
    from mvtpy import avdist

    gps_file = (WORKSPACE / "results" / "gps" / "CIRCLES_GPS_10Hz_2022-11-16.json")
    if not gps_file.is_file():
        pytest.skip("assembled GPS results not available")
    return avdist.load_av_runs(gps_file)


@pytest.fixture(scope="session")
def released_segments():
    """Released slim segments, keyed by their full trajectory id."""
    if not SLIM_SEGMENT.is_file():
        pytest.skip("released results/ not available")
    return {segment["trajectory_id"]["x_oid"]: segment
            for segment in _read_records(SLIM_SEGMENT, SLIM_RECORDS)}


@pytest.fixture(scope="session")
def produced_segments(raw_segment_path, driving_line):
    """Segments this port produces from the raw data, keyed by trajectory id."""
    from itertools import islice

    from mvtpy import lanes
    from mvtpy.rawio import iter_trajectories

    westbound = (record for record in iter_trajectories(raw_segment_path)
                 if record["direction"] < 0)

    produced = {}
    for record in islice(westbound, LANE_TRAJECTORIES):
        y_corr, lane = lanes.assign_lanes(record, driving_line)
        for segment in lanes.clip_lane_changes(record, y_corr, lane):
            produced[segment["trajectory_id"]] = segment
    return produced


@pytest.fixture(scope="session")
def grade_map():
    from mvtpy.kinematics import GradeMap

    if not GRADE_FIT.is_file():
        pytest.skip(f"missing {GRADE_FIT}")
    return GradeMap.from_csv(GRADE_FIT)


@pytest.fixture(scope="session")
def segment_pairs():
    """Pairs of (raw window, released segment) with identical sample times.

    A released segment is a lane-change-clipped piece of a raw trajectory, so
    the raw window is located by matching the segment's first timestamp and
    verifying that the whole rounded time vector agrees.
    """
    if not RAW_SEGMENT.is_file() or not SLIM_SEGMENT.is_file():
        pytest.skip("raw data/ or released results/ not available")

    raw_records = _read_records(RAW_SEGMENT, RAW_RECORDS)
    released = _read_records(SLIM_SEGMENT, SLIM_RECORDS)

    by_oid = {}
    for record in raw_records:
        identifier = record["_id"]
        by_oid[identifier["$oid"] if isinstance(identifier, dict) else identifier] = record

    pairs = []
    for segment in released:
        oid = segment["trajectory_id"]["x_oid"].rsplit("-", 1)[0]
        record = by_oid.get(oid)
        if record is None:
            continue

        raw_time = np.asarray(record["timestamp"], dtype=float)
        raw_x = np.asarray(record["x_position"], dtype=float)
        segment_time = np.asarray(segment["timestamp"], dtype=float)

        start = int(np.argmin(np.abs(raw_time - segment_time[0])))
        stop = start + len(segment_time)
        if stop > len(raw_time):
            continue
        if not np.array_equal(round_decimals(raw_time[start:stop], 4), segment_time):
            continue

        pairs.append({
            "segment": segment,
            "time": raw_time[start:stop],
            "x_feet": raw_x[start:stop],
            "direction": record["direction"],
        })

    if not pairs:
        pytest.skip("no raw/released segment pairs could be matched")
    return pairs
