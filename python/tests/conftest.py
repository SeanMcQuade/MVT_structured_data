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
